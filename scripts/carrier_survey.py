#!/usr/bin/env python3
"""Copy number of given sequences in every available k-mer index, to see who carries them.

For each sequence in a FASTA, take its canonical k-mers, look each up in every strain's index, divide by that
strain's single-copy depth, and report the median over the sequence's k-mers. A repeat shared by related
strains shows many carriers. Contamination or a library artifact shows one or two unrelated carriers.

usage: carrier_survey.py --fasta unitigs.fa --ids u1 u2 --workdir work --k 21 --out table.tsv
Needs jellyfish on PATH.
"""
import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mitekmer as mk  # noqa: E402
import discover_bursts as db  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--ids", nargs="*", help="sequence ids to use (default: all)")
    ap.add_argument("--workdir", default="work")
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--strains", nargs="*", help="default: every strain with an index and a baseline")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seqs = {n: q for n, q in mk.read_fasta(a.fasta)}
    ids = a.ids or list(seqs)
    k = a.k
    per_seq = {i: sorted({mk.canonical(seqs[i][j:j + k]) for j in range(len(seqs[i]) - k + 1)
                          if set(seqs[i][j:j + k]) <= set("ACGT")}) for i in ids}
    allk = sorted({x for v in per_seq.values() for x in v})
    pos = {x: n for n, x in enumerate(allk)}
    work = Path(a.workdir)
    strains = a.strains or sorted(p.name[:-len(f".k{k}.jf")] for p in work.glob(f"*.k{k}.jf"))
    rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        for s in strains:
            if not (work / f"{s}.k{k}.histo").exists():
                continue
            b = db.single_copy_depth(work, s, k)
            if not np.isfinite(b):
                continue
            cn = db.query_counts(work / f"{s}.k{k}.jf", allk, tmp) / b
            rows[s] = {i: float(np.median(cn[[pos[x] for x in per_seq[i]]])) for i in ids}
    tab = pd.DataFrame(rows).T
    tab.index.name = "strain"
    tab.round(2).to_csv(a.out, sep="\t")
    print(f"{len(tab)} strains x {len(ids)} sequences -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
