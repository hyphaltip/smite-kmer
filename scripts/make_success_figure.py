#!/usr/bin/env python3
"""Figure of mPing copy number for the accessions whose k-mer analysis succeeded.

Reads an `analyze` output directory and the Chen 2019 paper values, and writes:
  <out>/success_k<k>.png   two panels: paper vs k-mer estimate, and estimates by paper bin
  <out>/success_k<k>.tsv   per-accession table (paper value, estimate, group, bin)

Accessions with no baseline (no clear single-copy peak) are listed in the table with an empty
estimate. They are not plotted.

usage: make_success_figure.py --results results/success/k21 --k 21 --out docs/report/success
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
BINS = [(0, 0, "0"), (1, 5, "1-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 50, "21-50"),
        (51, 100, "51-100"), (101, 10**6, ">100")]
COL = {"burst": "#c0392b", "reference": "#34495e", "low_panel": "#2a7fb8", "dose": "#7f8c8d",
       "landrace_unverified": "#e08e0b"}


def bin_label(v):
    for lo, hi, lab in BINS:
        if lo <= v <= hi:
            return lab
    return None


def load_paper():
    grp = pd.read_csv(ROOT / "resources/strain_groups.tsv", sep="\t").set_index("strain")
    stream = pd.read_csv(ROOT / "resources/Chen2019/stream100.tsv", sep="\t").set_index("strain")
    rows = []
    for s in grp.index:
        rows.append((s, grp.loc[s, "group"], grp.loc[s, "paper_mping"]))
    for s in stream.index:
        if s not in grp.index:
            rows.append((s, "dose", stream.loc[s, "paper_mping"]))
    return pd.DataFrame(rows, columns=["strain", "group", "paper"]).drop_duplicates("strain").set_index("strain")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res, out = Path(a.results), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cn = pd.read_csv(res / f"te_copy_number.k{a.k}.tsv", sep="\t").set_index("te")
    base = pd.read_csv(res / f"baseline.k{a.k}.tsv", sep="\t", index_col=0)["single_copy_depth"]
    paper = load_paper()
    strains = [s for s in base.index if s in paper.index]
    rows = []
    for s in strains:
        est = cn.loc["mPing", f"cn_{s}"] if f"cn_{s}" in cn.columns else np.nan
        rows.append({"strain": s, "group": paper.loc[s, "group"], "paper_mPing": paper.loc[s, "paper"],
                     "single_copy_depth": base[s] if pd.notna(base[s]) else np.nan,
                     f"kmer_mPing_k{a.k}": est})
    tab = pd.DataFrame(rows)
    tab.to_csv(out / f"success_k{a.k}.tsv", sep="\t", index=False)

    ok = tab.dropna(subset=[f"kmer_mPing_k{a.k}", "paper_mPing"]).copy()
    x = ok.paper_mPing.to_numpy(float)
    y = ok[f"kmer_mPing_k{a.k}"].to_numpy(float)
    r = pearsonr(x, y)[0]
    rho = spearmanr(x, y)[0]
    n = len(ok)
    n_missing = int(tab[f"kmer_mPing_k{a.k}"].isna().sum())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.9), gridspec_kw={"width_ratios": [1.1, 1]})
    plt.rcParams.update({"font.size": 10})
    ax1.axhspan(231, 503, color="#c0392b", alpha=0.08, lw=0)
    lim = 220
    ax1.plot([0, lim], [0, lim], color="#999", lw=1, ls="--")
    for g in ["dose", "low_panel", "burst", "reference", "landrace_unverified"]:
        m = ok.group == g
        if m.any():
            ax1.scatter(ok.paper_mPing[m] + np.random.default_rng(1).uniform(-0.25, 0.25, m.sum()) * (g == "dose"),
                        ok[f"kmer_mPing_k{a.k}"][m], s=22 if g == "dose" else 46,
                        color=COL[g], alpha=0.8, label={"dose": "3K accessions", "low_panel": "low-mPing panel",
                                                        "burst": "burst (Chen 2019)", "reference": "Nipponbare / NB",
                                                        "landrace_unverified": "landrace (not in paper)"}[g],
                        zorder=3)
    ax1.set_xlabel("mPing copies, Chen 2019 read-depth estimate")
    ax1.set_ylabel(f"mPing copies, k-mer estimate (k={a.k})")
    ax1.set_xlim(-5, 190)
    ax1.set_ylim(-5, 560)
    ax1.set_title(f"{n} accessions with a baseline and a paper value", fontsize=10)
    ax1.text(0.03, 0.97, f"Pearson r = {r:.3f}\nSpearman = {rho:.3f}", transform=ax1.transAxes,
             va="top", fontsize=9)
    ax1.legend(frameon=False, fontsize=8, loc="lower right")
    ax1.spines[["top", "right"]].set_visible(False)

    labels, data, cols = [], [], []
    for lo, hi, lab in BINS:
        m = (ok.paper_mPing >= lo) & (ok.paper_mPing <= hi)
        if m.any():
            labels.append(f"{lab}\n(n={int(m.sum())})")
            data.append(ok[f"kmer_mPing_k{a.k}"][m].to_numpy(float))
    pos = np.arange(len(data))
    ax2.axhspan(231, 503, color="#c0392b", alpha=0.08, lw=0)
    for i, d in enumerate(data):
        j = np.random.default_rng(2).uniform(-0.18, 0.18, len(d))
        ax2.scatter(np.full(len(d), i) + j, d, s=20, color="#2a7fb8", alpha=0.7, zorder=3)
        ax2.plot([i - 0.25, i + 0.25], [np.median(d)] * 2, color="#222", lw=2, zorder=4)
    ax2.set_xticks(pos)
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_xlabel("paper mPing copies (bin)")
    ax2.set_ylabel(f"k-mer mPing copies (k={a.k})")
    ax2.set_yscale("symlog", linthresh=10)
    ax2.set_ylim(0, 600)
    ax2.set_title("k-mer estimate by paper bin (bar = median)", fontsize=10)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"success_k{a.k}.png", dpi=150)
    plt.close(fig)

    print(f"k={a.k}: plotted {n} accessions; {n_missing} with no estimate (no baseline)")
    print(f"Pearson r={r:.3f} Spearman={rho:.3f}")
    print(f"wrote {out / f'success_k{a.k}.png'} and {out / f'success_k{a.k}.tsv'}")


if __name__ == "__main__":
    sys.exit(main())
