import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import mitekmer as mk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_canonical_is_strand_independent():
    s = "ACGTTGCAAC"
    assert mk.canonical(s) == mk.canonical(mk.revcomp(s))


def test_kmer_table_multiplicity(tmp_path):
    fa = tmp_path / "te.fa"
    # second copy of the first 8 bases (as reverse complement) -> shared canonical k-mers
    seq = "ACGGTCAGTTCAGGCATCGA"
    fa.write_text(f">a\n{seq}\n>b\n{mk.revcomp(seq[:8])}TTTTTTT\n")
    df = mk.te_kmer_table(fa, 5)
    shared = df[df.te == "a"].head(4)
    assert (shared.multiplicity == 2).all()
    assert (df[df.te == "a"].tail(3).multiplicity == 1).all()


def test_baseline_peak_skips_error_peak():
    counts = np.arange(1, 61)
    freqs = np.round(1e6 * np.exp(-counts / 1.5)) + np.round(
        2e5 * np.exp(-0.5 * ((counts - 30) / 5.0) ** 2))
    assert abs(mk.baseline_peak(counts, freqs.astype(int)) - 30) <= 1


def rand_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def simulate_reads(genome_seqs, path, cov, rng, rl=100, err=0.003):
    total = sum(len(g) for g in genome_seqs)
    n = int(cov * total / (2 * rl))
    with gzip.open(path[0], "wt") as f1, gzip.open(path[1], "wt") as f2:
        for i in range(n):
            g = genome_seqs[rng.choice(len(genome_seqs), p=[len(x) / total for x in genome_seqs])]
            p = rng.integers(0, len(g) - 300)
            frag = g[p:p + 300]
            r1, r2 = frag[:rl], mk.revcomp(frag)[:rl]
            for rd, fh in ((r1, f1), (r2, f2)):
                rd = list(rd)
                for j in np.nonzero(rng.random(rl) < err)[0]:
                    rd[j] = rng.choice(list("ACGT"))
                fh.write(f"@r{i}\n{''.join(rd)}\n+\n{'I' * rl}\n")


@pytest.mark.skipif(shutil.which("jellyfish") is None, reason="jellyfish not installed")
def test_synthetic_copy_number_recovered(tmp_path):
    """Strain REF has 1 copy of each TE; strain AMP has TE_B at 10x. Both 30x backbone coverage."""
    rng = np.random.default_rng(7)
    backbone = rand_seq(rng, 200_000)
    te_a, te_b = rand_seq(rng, 800), rand_seq(rng, 800)
    (tmp_path / "te.fa").write_text(f">TE_A\n{te_a}\n>TE_B\n{te_b}\n")
    sheet = ["strain\tr1\tr2"]
    for name, copies_b in (("REF", 1), ("AMP", 10)):
        # one genome "molecule" per copy so fragments never span copies
        seqs = [backbone] + [te_a + rand_seq(rng, 400)] + [te_b + rand_seq(rng, 400)] * copies_b
        r1, r2 = tmp_path / f"{name}_1.fq.gz", tmp_path / f"{name}_2.fq.gz"
        simulate_reads(seqs, (r1, r2), 30, rng)
        sheet.append(f"{name}\t{r1}\t{r2}")
    (tmp_path / "samples.tsv").write_text("\n".join(sheet) + "\n")
    env_work = tmp_path / "work"
    for _, row in pd.read_csv(tmp_path / "samples.tsv", sep="\t").iterrows():
        subprocess.run([str(ROOT / "scripts/count.sh"), row.strain, row.r1, row.r2, "21", "2", "64M", "20"],
                       check=True, env={"WORKDIR": str(env_work), "TE_FASTA": str(tmp_path / "te.fa"),
                                        "PATH": __import__("os").environ["PATH"]})
    subprocess.run([sys.executable, str(ROOT / "scripts/mitekmer.py"), "analyze", "--fasta", str(tmp_path / "te.fa"),
                    "--samples", str(tmp_path / "samples.tsv"), "--workdir", str(env_work), "--k", "21",
                    "--reference", "REF", "--out", str(tmp_path / "out")], check=True)
    res = pd.read_csv(tmp_path / "out/te_copy_number.k21.tsv", sep="\t").set_index("te")
    assert res.loc["TE_A", "fc_AMP"] == pytest.approx(1.0, rel=0.25)
    assert res.loc["TE_B", "fc_AMP"] == pytest.approx(10.0, rel=0.25)
    assert not res.loc["TE_A", "spike_AMP"]
    assert res.loc["TE_B", "spike_AMP"]


def test_window_profile_flags_partial_copies():
    n = 200
    depths = pd.DataFrame({"te": "x", "pos": np.arange(n), "multiplicity": 1,
                           "depth_full": 100.0, "depth_partial": 100.0})
    depths.loc[80:119, "depth_partial"] = 10.0  # 40 k-mers missing from most copies
    # a k-mer shared with another TE must be ignored
    depths.loc[0, "multiplicity"] = 2
    depths.loc[0, "depth_full"] = 9999.0
    win, summ = mk.window_profile(depths, {"full": 10, "partial": 10}, window=40)
    s = summ.set_index("strain")
    assert s.loc["full", "median_cn"] == 10 and s.loc["full", "breadth"] == 1
    assert s.loc["partial", "breadth"] < 1
    assert s.loc["partial", "min_window_ratio"] == pytest.approx(0.1)
    assert win[(win.strain == "full")].cn.max() == 10  # multiplicity>1 k-mer dropped


def test_window_profile_zero_depth_is_nan_not_error():
    depths = pd.DataFrame({"te": "x", "pos": np.arange(50), "multiplicity": 1, "depth_s": 0.0})
    _, summ = mk.window_profile(depths, {"s": 10}, window=25)
    assert np.isnan(summ.breadth.iloc[0])


def test_gc_count():
    assert mk.gc_count("AAGCTTGC") == 4


def test_correct_gc_divides_by_bin_bias_and_skips_sparse_bins():
    depths = pd.DataFrame({"kmer": ["AAAAA", "GGGGG", "ACGTA"],
                           "depth_s": [8.0, 12.0, 10.0]})
    tab = pd.DataFrame({"gc": [0, 2, 5], "n_kmers": [5000, 10, 5000], "bias": [0.8, 5.0, 1.2]})
    out = mk.correct_gc(depths, {"s": tab})
    assert out.depth_s.tolist() == pytest.approx([10.0, 10.0, 10.0])  # 8/.8, 12/1.2, gc=2 sparse -> /1


@pytest.mark.skipif(shutil.which("jellyfish") is None, reason="needs jellyfish")
def test_gc_bias_table_from_index(tmp_path):
    rng = np.random.default_rng(3)
    # low-GC reads at depth ~2x the GC-rich ones
    lo = rand_seq(rng, 4000).replace("G", "A").replace("C", "T")
    hi = rand_seq(rng, 4000)
    fq = tmp_path / "r.fa"
    with open(fq, "w") as fh:
        for name, g, cov in (("lo", lo, 40), ("hi", hi, 20)):
            for i in range(cov * len(g) // 60):
                st = int(rng.integers(0, len(g) - 60))
                fh.write(f">{name}{i}\n{g[st:st + 60]}\n")
    idx = tmp_path / "x.jf"
    subprocess.run(["jellyfish", "count", "-m", "15", "-s", "1M", "-C", "-o", str(idx), str(fq)], check=True)
    tab = mk.gc_bias_table(idx, 15, baseline=20, lo=0.3, hi=3.0)
    assert tab.loc[tab.gc <= 1, "mean_depth"].mean() > 1.3 * tab.loc[tab.gc >= 6, "mean_depth"].mean()


def _hist(freq_fn, overflow=0, top=1000):
    counts = np.arange(1, top + 2)
    freqs = np.array([freq_fn(c) for c in counts[:-1]] + [overflow], dtype=float)
    return counts, freqs


def test_baseline_peak_ignores_overflow_bin():
    # jellyfish histo -h 1000 puts every k-mer with a higher count in bin 1001
    counts, freqs = _hist(lambda c: 1e6 * np.exp(-c / 1.5) + 2e5 * np.exp(-0.5 * ((c - 30) / 5.0) ** 2),
                          overflow=5e6)
    assert abs(mk.baseline_peak(counts, freqs) - 30) <= 1


def test_baseline_peak_nan_when_no_peak():
    # shallow strain: errors and single-copy peak merge into one falling curve with a shoulder
    head = [40162480, 34617967, 37328309, 36711675, 33117262, 27784742, 21910363,
            16417589, 11758585, 8162998, 5542059, 3704232, 2488393, 1693122]
    counts = np.arange(1, 1002)
    freqs = np.array(head + [head[-1] * 0.68 ** i for i in range(1, 987)] + [155638], dtype=float)
    assert np.isnan(mk.baseline_peak(counts, freqs))


def test_baseline_peak_nan_for_monotone_decline():
    counts, freqs = _hist(lambda c: 1e7 * np.exp(-c / 20.0), overflow=1e3)
    assert np.isnan(mk.baseline_peak(counts, freqs))


def test_baseline_peak_deep_strain_still_found():
    counts, freqs = _hist(lambda c: 1e6 * np.exp(-c / 1.5) + 8e5 * np.exp(-0.5 * ((c - 133) / 15.0) ** 2),
                          overflow=1e5)
    assert abs(mk.baseline_peak(counts, freqs) - 133) <= 2
