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
