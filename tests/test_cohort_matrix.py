import shutil
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import cohort_matrix as cm  # noqa: E402
import mitekmer as mk  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("jellyfish") is None or shutil.which("zstd") is None,
                                reason="needs jellyfish and zstd")


def rand_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def test_standout_mask_and_carrier_counts():
    cn = np.array([100.0, 100.0, 10.0, 100.0])
    med = np.array([1.0, 30.0, 0.0, 5.0])
    m = cm.standout_mask(cn, med, min_cn=20, min_fold=5)
    assert m.tolist() == [True, False, False, True]    # 100 vs 31 is fold 3.2, 10 is below min_cn
    M = np.array([[100.0, 60.0, 0.0, 1.0],   # one other strain at 60 >= max(3, 25)
                  [100.0, 2.0, 1.0, 0.0]])   # nobody else carries it
    assert cm.carrier_counts(M, np.array([100.0, 100.0]), 0).tolist() == [1, 0]


def test_end_to_end_private_burst_versus_lineage_repeat(tmp_path):
    rng = np.random.default_rng(21)
    k, depth = 21, 30
    genome = rand_seq(rng, 12000)
    burst, lineage, common = rand_seq(rng, 300), rand_seq(rng, 300), rand_seq(rng, 300)
    strains = {            # copies of (burst, lineage, common)
        "T": (60, 50, 40), "L": (1, 50, 40), "c1": (1, 0, 40), "c2": (2, 0, 40),
        "c3": (1, 0, 40), "c4": (0, 0, 40)}
    work = tmp_path / "work"
    work.mkdir()

    def element_reads(r, seq, copies, fh, n0):
        n = n0
        for _ in range(int(len(seq) * copies * depth / 100)):
            st = int(r.integers(-99, len(seq)))
            pad_l = "".join(r.choice(list("ACGT"), max(0, -st)))
            pad_r = "".join(r.choice(list("ACGT"), 100))
            fh.write(f">r{n}\n{(pad_l + seq[max(0, st):] + pad_r)[:100]}\n")
            n += 1
        return n

    for i, (name, (cb, cl, cc)) in enumerate(strains.items()):
        r = np.random.default_rng(100 + i)
        fa = work / f"{name}.fa"
        with open(fa, "w") as fh:
            n = 0
            for _ in range(int(len(genome) * depth / 100)):
                st = int(r.integers(0, len(genome) - 100))
                fh.write(f">r{n}\n{genome[st:st + 100]}\n")
                n += 1
            for seq, c in ((burst, cb), (lineage, cl), (common, cc)):
                if c:
                    n = element_reads(r, seq, c, fh, n)
        jf = work / f"{name}.k{k}.jf"
        subprocess.run(["jellyfish", "count", "-m", str(k), "-s", "20M", "-C", "-t", "2", "-o", str(jf), str(fa)],
                       check=True)
        (work / f"{name}.k{k}.histo").write_text(subprocess.run(
            ["jellyfish", "histo", "-h", "1000", str(jf)], check=True, capture_output=True, text=True).stdout)
        fa.unlink()

    out = tmp_path / "matrix"
    base = dict(workdir=str(work), k=k, out=str(out), part=0, nparts=1, tmpdir=str(tmp_path), min_cn=20.0,
                min_fold=5.0, max_carriers=None, strain="T", min_depth=15.0, allow_low_depth=False,
                edge_ratio=0.1, min_unitig=60, mask=[], fasta=None, element="burst", tag="")
    A = lambda **kw: Namespace(**{**base, **kw})
    cm.stage_strains(A())
    assert len(pd.read_csv(out / "strains.tsv", sep="\t")) == 6
    # two parts, as a SLURM array would run them
    for part in (0, 1):
        cm.stage_dump(A(part=part, nparts=2))
    cm.stage_union(A())
    for part in (0, 1):
        cm.stage_query(A(part=part, nparts=2))
    cm.stage_build(A())
    M, med = cm.load_matrix(out)
    assert M.shape[1] == 6 and M.shape[0] == len(med) > 500

    lib = tmp_path / "lib.fa"
    lib.write_text(f">burst\n{burst}\n>lineage\n{lineage}\n>common\n{common}\n")
    tab = cm.stage_call(A(fasta=str(lib)))
    matches = set(tab.best_library_match)
    assert "burst" in matches and "lineage" in matches and "common" not in matches
    b = tab[tab.best_library_match == "burst"].iloc[0]
    assert b.n_carriers == 0 and abs(b.length - 300) <= 2
    assert tab[tab.best_library_match == "lineage"].iloc[0].n_carriers == 1   # strain L also carries it
    # asking for strain-private sequences only drops the lineage repeat
    private = cm.stage_call(A(fasta=str(lib), max_carriers=0))
    assert set(private.best_library_match) == {"burst"}

    ev = cm.stage_evaluate(A(fasta=str(lib), element="burst")).set_index("strain")
    assert ev.loc["T", "fraction_called"] > 0.95
    assert ev.drop(index="T").fraction_called.max() < 0.05
