import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import discover_bursts as db  # noqa: E402
import mitekmer as mk  # noqa: E402


def rand_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def kmer_counts(seq, k, n):
    return {mk.canonical(seq[i:i + k]): n for i in range(len(seq) - k + 1)}


def test_unitig_recovers_one_sequence_either_strand():
    rng = np.random.default_rng(1)
    s = rand_seq(rng, 200)
    for query in (s, mk.revcomp(s)):
        u = db.build_unitigs(kmer_counts(query, 21, 100))
        assert len(u) == 1
        assert u[0][0] in (s, mk.revcomp(s))


def test_error_variants_are_pruned_not_split():
    rng = np.random.default_rng(2)
    s = rand_seq(rng, 300)
    counts = kmer_counts(s, 21, 1000)
    # a read with one error at position 150 adds 21 rare k-mers
    bad = s[:150] + ("A" if s[150] != "A" else "C") + s[151:]
    for km in kmer_counts(bad, 21, 40):
        counts.setdefault(km, 40)
    unitigs = db.build_unitigs(counts, ratio=0.1)
    longest = max(unitigs, key=lambda u: len(u[0]))[0]
    assert longest in (s, mk.revcomp(s))


def test_two_unrelated_sequences_give_two_unitigs():
    rng = np.random.default_rng(3)
    a, b = rand_seq(rng, 150), rand_seq(rng, 180)
    counts = {**kmer_counts(a, 21, 50), **kmer_counts(b, 21, 50)}
    u = db.build_unitigs(counts)
    assert sorted(len(x[0]) for x in u) == [150, 180]


def test_cycle_does_not_loop_forever():
    # a circular sequence: its k-mers form a cycle
    rng = np.random.default_rng(4)
    s = rand_seq(rng, 80)
    circ = s + s[:20]
    u = db.build_unitigs(kmer_counts(circ, 21, 30))
    assert len(u) >= 1


@pytest.mark.skipif(shutil.which("jellyfish") is None, reason="needs jellyfish")
def test_end_to_end_finds_burst_and_drops_common_family(tmp_path):
    rng = np.random.default_rng(5)
    k = 21
    genome = rand_seq(rng, 20000)
    burst = rand_seq(rng, 300)     # high only in the test strain
    adapter = rand_seq(rng, 70)    # high only in the test strain, like adapter read-through
    common = rand_seq(rng, 300)    # high in every strain

    def make(name, burst_copies, depth=30, seed=0, adapter_copies=0):
        r = np.random.default_rng(seed)
        fa = tmp_path / f"{name}.fa"
        with open(fa, "w") as fh:
            n = 0
            for seq, copies in ((genome, 1), (burst, burst_copies), (common, 40), (adapter, adapter_copies)):
                if copies == 0:
                    continue
                reads = int(len(seq) * copies * depth / 100)
                for _ in range(reads):
                    if copies == 1 and seq is genome:
                        st = int(r.integers(0, len(seq) - 100))
                        read = list(seq[st:st + 100])
                    else:
                        # an insertion sits in flanking DNA: reads may overhang either end
                        st = int(r.integers(-99, len(seq)))
                        pad_l = "".join(r.choice(list("ACGT"), max(0, -st)))
                        pad_r = "".join(r.choice(list("ACGT"), 100))
                        read = list((pad_l + seq[max(0, st):] + pad_r)[:100])
                    if r.random() < 0.3:      # one error in 30% of reads
                        p = int(r.integers(0, 100))
                        read[p] = "ACGT"[(("ACGT".index(read[p])) + 1) % 4]
                    fh.write(f">r{n}\n{''.join(read)}\n")
                    n += 1
        jf = tmp_path / f"{name}.k{k}.jf"
        subprocess.run(["jellyfish", "count", "-m", str(k), "-s", "20M", "-C", "-t", "2", "-o", str(jf), str(fa)],
                       check=True)
        histo = subprocess.run(["jellyfish", "histo", "-h", "1000", str(jf)], check=True,
                               capture_output=True, text=True).stdout
        (tmp_path / f"{name}.k{k}.histo").write_text(histo)

    make("pos", burst_copies=60, seed=1, adapter_copies=200)
    make("c1", burst_copies=2, seed=2)
    make("c2", burst_copies=1, seed=3)
    lib = tmp_path / "lib.fa"
    lib.write_text(f">burst\n{burst}\n>common\n{common}\n")
    out = tmp_path / "out"
    mask = tmp_path / "mask.fa"
    mask.write_text(f">adapter\n{adapter}\n")
    tab = db.run(type("A", (), dict(positive="pos", cohort=["c1", "c2"], workdir=str(tmp_path), k=k,
                                    min_cn=20.0, min_fold=10.0, extend_fold=0.0, min_seed_fraction=0.25,
                                    min_depth=15.0, allow_low_depth=False, threads=2, 
                                    edge_ratio=0.1, min_unitig=60, mask=[str(mask)],
                                    fasta=str(lib), tmpdir=str(tmp_path), out=str(out)))())
    assert len(tab) >= 1
    top = tab.sort_values("length", ascending=False).iloc[0]
    assert top.best_library_match == "burst"
    assert top.fraction_in_library > 0.95
    assert abs(top.length - len(burst)) <= 2
    assert not (tab.best_library_match == "common").any()
    # the adapter is reported, but flagged, and the burst is not
    assert bool(tab[tab.artifact].shape[0] >= 1)
    assert not bool(top.artifact)
    assert top.n_carriers == 0
    assert list(tab.score) == sorted(tab.score, reverse=True)
    assert (out / "kmers.tsv.zst").exists() or (out / "kmers.tsv.gz").exists()


def test_bubble_from_common_second_allele_is_popped():
    rng = np.random.default_rng(7)
    s = rand_seq(rng, 300)
    counts = kmer_counts(s, 21, 1000)
    alt = s[:150] + ("A" if s[150] != "A" else "C") + s[151:]
    for km in kmer_counts(alt, 21, 200):         # 17% second allele: above the 0.1 edge ratio
        counts.setdefault(km, 200)
    without = db.build_unitigs(counts, ratio=0.1, pop=False)
    with_pop = db.build_unitigs(counts, ratio=0.1, pop=True)
    assert max(len(u[0]) for u in without) < 300       # split at the SNP
    assert max(len(u[0]) for u in with_pop) == 300
    assert max(with_pop, key=lambda u: len(u[0]))[0] in (s, mk.revcomp(s))


def test_two_real_paralogs_far_apart_are_not_popped_as_a_bubble():
    rng = np.random.default_rng(8)
    shared_l, shared_r = rand_seq(rng, 100), rand_seq(rng, 100)
    a, b = rand_seq(rng, 100), rand_seq(rng, 100)       # arms of 100 bp: longer than 2k = 42
    counts = {}
    for arm in (a, b):
        counts.update(kmer_counts(shared_l + arm + shared_r, 21, 500))
    g = db.Graph(dict(counts), 0.1)
    assert g.pop_bubbles() == 0


def test_even_k_is_rejected():
    with pytest.raises(ValueError):
        db.Graph({"ACGT": 5}, 0.1)


def test_low_depth_test_strain_is_refused(tmp_path):
    (tmp_path / "t.k21.histo").write_text("".join(f"{c} {int(1e6 * (0.5 ** c) + 5e4 * np.exp(-0.5 * ((c - 8) / 2) ** 2))}\n"
                                                 for c in range(1, 1002)))
    args = type("A", (), dict(positive="t", cohort=["c"], workdir=str(tmp_path), k=21, min_cn=20.0,
                              min_fold=10.0, extend_fold=0.0, min_seed_fraction=0.25, min_depth=15.0,
                              allow_low_depth=False, threads=2, edge_ratio=0.1,
                              min_unitig=60, mask=[], fasta=None, tmpdir=str(tmp_path), out=str(tmp_path / "o")))()
    with pytest.raises(SystemExit):
        db.run(args)


@pytest.mark.skipif(shutil.which("jellyfish") is None, reason="needs jellyfish")
def test_hysteresis_extends_a_seed_region(tmp_path):
    rng = np.random.default_rng(11)
    k = 21
    genome = rand_seq(rng, 20000)
    elem = rand_seq(rng, 800)
    depth = 30

    def make(name, copies_by_region, seed):
        """One element in random flanks. The number of reads that start in each region sets its copy level."""
        r = np.random.default_rng(seed)
        fa = tmp_path / f"{name}.fa"
        with open(fa, "w") as fh:
            n = 0
            for _ in range(int(len(genome) * depth / 100)):
                st = int(r.integers(0, len(genome) - 100))
                fh.write(f">r{n}\n{genome[st:st + 100]}\n")
                n += 1
            for (lo, hi), copies in copies_by_region:
                for _ in range(int((hi - lo) * copies * depth / 100)):
                    st = int(r.integers(lo, hi))
                    pad_l = "".join(r.choice(list("ACGT"), max(0, -st)))
                    pad_r = "".join(r.choice(list("ACGT"), 100))
                    read = (pad_l + elem[max(0, st):st + 100] + pad_r)[:100]
                    fh.write(f">r{n}\n{read}\n")
                    n += 1
        jf = tmp_path / f"{name}.k{k}.jf"
        subprocess.run(["jellyfish", "count", "-m", str(k), "-s", "20M", "-C", "-t", "2", "-o", str(jf), str(fa)],
                       check=True)
        (tmp_path / f"{name}.k{k}.histo").write_text(subprocess.run(
            ["jellyfish", "histo", "-h", "1000", str(jf)], check=True, capture_output=True, text=True).stdout)

    # Reads are 100 bp, so a change in copy level blurs over about 80 bp. The 500 bp middle is long
    # enough to reach a plateau. Reads may start 99 bp before the element.
    regions = lambda ends, mid: [((-99, 150), ends), ((150, 650), mid), ((650, 800), ends)]
    make("pos", regions(80, 80), 1)
    make("c1", regions(14, 0.5), 2)
    make("c2", regions(14, 0.5), 3)

    def go(min_fold, extend_fold, tag):
        a = type("A", (), dict(positive="pos", cohort=["c1", "c2"], workdir=str(tmp_path), k=k, min_cn=20.0,
                               min_fold=min_fold, extend_fold=extend_fold, min_seed_fraction=0.25, min_depth=15.0,
                               allow_low_depth=False, threads=2, edge_ratio=0.1, min_unitig=60, mask=[], fasta=None,
                               tmpdir=str(tmp_path), out=str(tmp_path / tag)))()
        return db.run(a)

    strict = go(30.0, 0.0, "strict")     # the ends have fold ~5, the middle ~40: only the middle is a seed
    hyst = go(30.0, 4.0, "hyst")         # relaxed neighbours rejoin the ends to the seed region
    assert len(strict) >= 1 and strict.length.max() < 500
    assert hyst.length.max() >= 780
