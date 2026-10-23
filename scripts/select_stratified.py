#!/usr/bin/env python3
"""Pick a stratified set of Chen et al. 2019 accessions across the mPing range.

Joins Supplementary Data 1 to the ENA PRJEB6180 run list. Keeps accessions with
approximate depth >= --min-depth whose runs all have exactly two FASTQ files, and
that are not already in --exclude. Samples a fixed number per mPing bin (fixed seed).
Writes TSV: strain, mping, subgroup, approx_depth, gb, r1_urls, r2_urls (comma-separated).
"""
import argparse
import numpy as np
import pandas as pd
from select_panel import load_table, GENOME_MB

BINS = [("0", 0, 0, 12), ("1-5", 1, 5, 14), ("6-10", 6, 10, 14), ("11-20", 11, 20, 14),
        ("21-50", 21, 50, 20), ("51-100", 51, 100, 24), (">100", 101, 10**6, 3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--ena-runs", required=True)
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--min-depth", type=float, default=10.0)
    ap.add_argument("--max-gb", type=float, default=15.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    t = load_table(a.xlsx).dropna(subset=["mPing"])
    e = pd.read_csv(a.ena_runs, sep="\t").dropna(subset=["fastq_ftp"])
    e["files"] = e.fastq_ftp.str.split(";")
    e["sizes"] = e.fastq_bytes.astype(str).str.split(";").map(lambda v: sum(int(x) for x in v if x.isdigit()))
    e["paired"] = e.files.map(len) == 2
    g = e.groupby("sample_alias")
    info = pd.DataFrame({"all_paired": g["paired"].all(), "gb": g["sizes"].sum() / 1e9,
                         "depth": g["base_count"].sum() / (GENOME_MB * 1e6)})
    r1 = e[e.paired].groupby("sample_alias")["files"].apply(lambda s: ",".join(f[0] for f in s))
    r2 = e[e.paired].groupby("sample_alias")["files"].apply(lambda s: ",".join(f[1] for f in s))
    info["r1"], info["r2"] = r1, r2
    t = t.join(info, on="Taxa")
    ok = t[(t.all_paired == True) & (t.depth >= a.min_depth) & (t.gb <= a.max_gb)  # noqa: E712
           & ~t.Taxa.isin(a.exclude)]
    rng = np.random.default_rng(a.seed)
    picks = []
    for name, lo, hi, n in BINS:
        pool = ok[(ok.mPing >= lo) & (ok.mPing <= hi)]
        take = pool.iloc[rng.permutation(len(pool))[:n]]
        print(f"bin {name:7} pool={len(pool):4d} picked={len(take)}")
        picks.append(take)
    sel = pd.concat(picks)
    out = sel[["Taxa", "mPing", "Subgroups", "depth", "gb", "r1", "r2"]].rename(
        columns={"Taxa": "strain", "mPing": "paper_mping", "Subgroups": "subgroup",
                 "depth": "approx_depth"})
    out.to_csv(a.out, sep="\t", index=False)
    print(f"{len(out)} accessions, {out.gb.sum():.0f} GB to stream, "
          f"subgroups: {out.subgroup.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
