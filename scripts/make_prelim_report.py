#!/usr/bin/env python3
"""Build the preliminary report files from an `analyze` output directory.

Writes to docs/report/ (tracked by git, unlike results/):
  REPORT.md            report with relative links to the figures and the PDF
  REPORT.pdf           same content, built with pandoc + xelatex
  report_preview.png   2x2 montage of the figures, for embedding with a link to the PDF
  figures/*.png        the four figures
  mping_by_strain.tsv  per-strain table used in the report

usage: make_prelim_report.py [--results results/full] [--k 21] [--out docs/report] [--no-pdf]
Needs pandoc and xelatex on PATH for the PDF.
"""
import argparse
import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import mitekmer as mk  # noqa: E402

COL = {"burst": "#c0392b", "reference": "#34495e", "low_panel": "#2a7fb8", "landrace_unverified": "#e08e0b"}
LAB = {"burst": "burst (Chen 2019)", "reference": "Nipponbare reference", "low_panel": "low-mPing panel"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def load(results, k):
    cn = pd.read_csv(results / f"te_copy_number.k{k}.tsv", sep="\t").set_index("te")
    depths = pd.read_csv(results / f"te_kmer_depths.k{k}.tsv", sep="\t")
    base = pd.read_csv(results / f"baseline.k{k}.tsv", sep="\t", index_col=0)["single_copy_depth"]
    grp = pd.read_csv(ROOT / "resources/strain_groups.tsv", sep="\t").set_index("strain")
    strains = [s for s in base.index if s in grp.index]
    return cn, depths, base, grp, strains


def fig_recovery(mp, grp, k, path):
    fig, ax = plt.subplots(figsize=(5.6, 5))
    pp = grp.loc[mp.index, "paper_mping"].dropna()
    pp = pp[mp[pp.index].notna()]
    x, y = pp.values, mp[pp.index].values
    lim = 60
    ax.plot([0, lim], [0, lim], color="#999", lw=1, ls="--", label="y = x")
    for s in pp.index:
        ax.scatter(pp[s], mp[s], s=48, color=COL[grp.loc[s, "group"]], zorder=3)
    ax.set_xlabel("mPing copies, Chen 2019 read-depth estimate")
    ax.set_ylabel(f"mPing copies, k-mer estimate (k={k})")
    low = (grp.loc[pp.index, "group"] == "low_panel").values
    r, rho = pearsonr(x[low], y[low])[0], spearmanr(x[low], y[low])[0]
    ax.text(0.04, 0.95, f"{low.sum()} low-mPing strains: Pearson r = {r:.3f}\nSpearman = {rho:.3f}",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.scatter([], [], color=COL["low_panel"], label="low-mPing panel")
    ax.scatter([], [], color=COL["reference"], label="NB")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return r, rho, x[low], y[low]


def fig_by_strain(mp, grp, path):
    order = mp.dropna().sort_values().index
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axhspan(231, 503, color="#c0392b", alpha=0.10, lw=0)
    ax.text(len(mp) - 3.6, 340, "Chen 2019 burst range\n(231-503)", color="#a93226",
            fontsize=8.5, va="center", ha="right")
    for i, s in enumerate(order):
        h = max(mp[s], 0.3)
        ax.bar(i, h, color=COL[grp.loc[s, "group"]], width=0.7)
        ax.text(i, h * 1.12, f"{mp[s]:.0f}" if mp[s] >= 10 else f"{mp[s]:.1f}", ha="center", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 1500)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("mPing copy number (k-mer estimate)")
    for g in ["burst", "reference", "low_panel"]:
        ax.bar(np.nan, np.nan, color=COL[g], label=LAB[g])
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    missing = list(mp[mp.isna()].index)
    if missing:
        fig.text(0.01, 0.005, "No baseline, not shown: " + ", ".join(missing), ha="left", va="bottom",
                 fontsize=8, color="#555")
    fig.tight_layout(rect=(0, 0.03 if missing else 0, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return order


def fig_profile(depths, base, order, fasta, path):
    win, _ = mk.window_profile(depths, base.to_dict(), 40)
    w = win[win.te == "mPing"]
    M = w.pivot(index="strain", columns="start", values="cn")
    M = M.loc[[s for s in order if s in M.index]]
    L = np.log2(M.div(M.median(axis=1), axis=0).replace(0, np.nan))
    seq = dict(mk.read_fasta(fasta))["mPing"]
    starts = list(M.columns)
    ends = w[w.strain == M.index[0]].sort_values("start")["end"].tolist()
    gc = [mk.gc_count(seq[a:b + 21]) / len(seq[a:b + 21]) for a, b in zip(starts, ends)]
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(7.4, 5.4),
                                 gridspec_kw={"height_ratios": [1, 6], "hspace": 0.08}, sharex=True)
    a0.plot(range(len(gc)), gc, marker="o", color="#555")
    a0.set_ylabel("GC", fontsize=9)
    a0.set_ylim(0.2, 0.8)
    im = a1.imshow(L.values, aspect="auto", cmap="RdBu", vmin=-2.5, vmax=2.5)
    a1.set_yticks(range(len(L)))
    a1.set_yticklabels(L.index, fontsize=8)
    a1.set_xticks(range(len(starts)))
    a1.set_xticklabels([f"{a}-{b}" for a, b in zip(starts, ends)], rotation=45, ha="right", fontsize=8)
    a1.set_xlabel("mPing position (bp)")
    cb = fig.colorbar(im, ax=[a0, a1], fraction=0.03, pad=0.02)
    cb.set_label("log2 (window copy number / strain median)", fontsize=8)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame({"start": starts, "end": ends, "gc": gc})


def fig_all_tes(cn, strains, order, path):
    tes = list(cn.index)
    H = pd.DataFrame({s: [cn.loc[t, f"cn_{s}"] for t in tes] for s in strains}, index=tes)
    Hl = np.log10(H.replace(0, np.nan).clip(lower=0.1))
    cols = [s for s in order if s in Hl.columns]
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    im = ax.imshow(Hl[cols].values, aspect="auto", cmap="viridis", vmin=-1, vmax=2.7)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(tes)))
    ax.set_yticklabels([f"{t}  (n={int(cn.loc[t, 'n_kmers_unique'])})" for t in tes], fontsize=8)
    for i, t in enumerate(tes):
        for j, s in enumerate(cols):
            v = H.loc[t, s]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:.0f}" if v >= 10 else f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                    color="white" if Hl.loc[t, s] < 1.2 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("log10 copy number", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return H


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(v) for v in r) + " |")
    return "\n".join(lines)


def montage(paths, out):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.6))
    for ax, p in zip(axes.flat, paths):
        ax.imshow(plt.imread(p))
        ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.03, hspace=0.03)
    fig.savefig(out, dpi=110, facecolor="white")
    plt.close(fig)


def git_state():
    def run(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = " (uncommitted changes present)" if run("status", "--porcelain") else ""
    return run("rev-parse", "--short", "HEAD") + dirty


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=str(ROOT / "results/full"))
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--out", default=str(ROOT / "docs/report"))
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--tex-bin", default=os.environ.get("TEX_BIN", ""),
                    help="directory with a TeX Live xelatex to put first on PATH (env TEX_BIN)")
    a = ap.parse_args()
    out, figs = Path(a.out), Path(a.out) / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    cn, depths, base, grp, strains = load(Path(a.results), a.k)
    mp = pd.Series({s: cn.loc["mPing", f"cn_{s}"] for s in strains})

    r, rho, xl, yl = fig_recovery(mp, grp, a.k, figs / "fig1_recovery.png")
    order = fig_by_strain(mp, grp, figs / "fig2_mping_by_strain.png")
    no_base = [s_ for s_ in mp.index if s_ not in order]
    order_all = pd.Index(list(order[::-1]) + no_base)  # high to low, strains without a baseline last
    gcw = fig_profile(depths, base, order_all, ROOT / "RiceTE.fa", figs / "fig3_mping_profile.png")
    H = fig_all_tes(cn, strains, order_all, figs / "fig4_all_tes.png")
    montage([figs / f for f in ("fig1_recovery.png", "fig2_mping_by_strain.png",
                                "fig3_mping_profile.png", "fig4_all_tes.png")], out / "report_preview.png")

    rows = []
    for s in list(order[::-1]) + no_base:
        paper = grp.loc[s, "paper_mping"]
        rows.append({"strain": s, "group": grp.loc[s, "group"], "single_copy_depth": "" if pd.isna(base[s]) else int(base[s]),
                     "paper_mPing": "" if pd.isna(paper) else f"{paper:.0f}",
                     f"kmer_mPing_k{a.k}": "n/a (no baseline)" if pd.isna(mp[s]) else f"{mp[s]:.1f}",
                     "ratio": "" if pd.isna(paper) or paper == 0 or pd.isna(mp[s]) else f"{mp[s] / paper:.2f}"})
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "mping_by_strain.tsv", sep="\t", index=False)
    ratios = np.array([float(x) for x in tab[tab.group == "low_panel"].ratio if x != ""])
    others = [t for t in H.index if t not in ("mPing", "nDart")]
    spread = (H.loc[others].max(axis=1) / H.loc[others].min(axis=1)).round(1)
    mp_pos = mp[mp > 0]
    hi_gc = gcw.loc[gcw.gc.idxmax()]
    lowest = []
    _, summ = mk.window_profile(depths, base.to_dict(), 40)
    sm = summ[(summ.te == "mPing")].set_index("strain")["min_window_ratio"]
    lowest = sm[(sm < 0.4) & (mp > 1)].round(2)

    md = f"""# mPing k-mer screen: preliminary results

[![Report preview](report_preview.png)](REPORT.pdf)

*Click the image for the [PDF](REPORT.pdf). Generated {datetime.date.today()} from code at commit {git_state()}.*

Copy numbers of rice transposable elements (TEs) estimated from k-mer depth in short reads, compared with
published values. The algorithm is described in [`docs/ALGORITHM.md`](../ALGORITHM.md).

**Scope.** k={a.k}; {len(strains)} strains{(" (" + ", ".join(no_base) + " has no baseline and no copy numbers)") if no_base else ""}; {len(cn)} TEs in `RiceTE.fa`; no GC correction applied. NB_DRR028132 uses only the first 2.27 GB of R1 (single-end, about 7.6x), because its download was incomplete.

## Summary

- For the {len(xl)} low-mPing strains, the k-mer mPing estimate and the Chen et al. 2019 read-depth estimate
  have Pearson r = {r:.3f} and Spearman {rho:.3f}. The k-mer value is a median {np.median(ratios):.2f} times the paper's
  (range {ratios.min():.2f} to {ratios.max():.2f}).
- Nipponbare (NB) gives {mp['NB']:.1f} mPing copies. The paper lists 52 for its Nipponbare entry.
- HEG4 gives {mp['HEG4']:.0f}. The paper reports 231 to 503 for its four bursting accessions and no value for HEG4 alone.

## Figure 1. mPing copies against the paper

![mPing k-mer estimate against the Chen 2019 estimate](figures/fig1_recovery.png){{width=60%}}

With {len(xl)} points, the correlation is a small-sample result.

## Figure 2. mPing copies by strain

![mPing copy number by strain](figures/fig2_mping_by_strain.png){{width=85%}}

The low-mPing panel spans {mp[[s for s in strains if grp.loc[s,'group']=='low_panel']].min():.0f} to {mp[[s for s in strains if grp.loc[s,'group']=='low_panel']].max():.0f} copies, NB has {mp['NB']:.0f} and HEG4 has {mp['HEG4']:.0f}.
Fold change against NB understates how low the panel is, because NB carries {mp['NB']:.0f} copies.
The HEG4 value has no independent check.

{md_table(tab)}

## Figure 3. Position profile along mPing

![mPing window copy number along the element](figures/fig3_mping_profile.png){{width=85%}}

Windows are 40 consecutive unique k-mers. The window at positions {int(hi_gc.start)} to {int(hi_gc.end)} is the highest-GC stretch
of mPing (GC {hi_gc.gc:.2f}, against {sum(mk.gc_count(x) for x in [dict(mk.read_fasta(ROOT / 'RiceTE.fa'))['mPing']])/len(dict(mk.read_fasta(ROOT / 'RiceTE.fa'))['mPing']):.2f} for the whole element).
Several strains show a dip in this window. The lowest window is below 0.4 of the strain median in
{', '.join(f'{s} ({v})' for s, v in lowest.items())}. A dip shared across strains at the highest-GC window suggests a
sequencing bias. This has not been tested. The GC correction is implemented in `scripts/mitekmer.py` and has not yet
been applied to these results. A strain with no mPing k-mers has an empty row.

## Figure 4. All TEs in the library

![Copy number of every TE in every strain](figures/fig4_all_tes.png){{width=90%}}

For the TEs other than mPing and nDart, the largest-to-smallest ratio across strains ranges from {spread.min()} to {spread.max()}.
For mPing the ratio is about {mp_pos.max() / mp_pos.min():.0f} ({mp_pos.max():.0f} against {mp_pos.min():.1f}, leaving out strains with 0 copies).
nDart: {H.loc['nDart'].max():.1f} in {H.loc['nDart'].idxmax()} against {H.loc['nDart'].drop(H.loc['nDart'].idxmax()).min():.1f} to {H.loc['nDart'].drop(H.loc['nDart'].idxmax()).max():.1f} elsewhere.
There are no published copy numbers for these TEs in the data used here, so the remaining spread cannot be assigned to biology or noise.
`n` is the number of unique k-mers per TE.

## Limits

- Each strain is scaled by its own single-copy depth, taken as the highest peak above the error valley of its k-mer histogram.
  Low-depth strains (single-copy depth 5 to 9) sit close to that valley, so their scaling is less certain.
- A k-mer must match exactly, so diverged copies are undercounted.
- A strain whose k-mer histogram has no clear single-copy peak gets no copy numbers instead of a guess
  (the peak must be at least 1.1 times the error valley, with the valley at count 200 or below).
- See [`docs/ALGORITHM.md`](../ALGORITHM.md) for the full list.

## Not yet included

- GC-corrected copy numbers. The correction is in `scripts/mitekmer.py`, but it is not applied to the figures here.
- Per-line comparison for the 100-line stratified set. See `docs/report/success/REPORT.md`.
- A119, A123 and A157 have no per-line value in Chen et al. 2019, so they appear without a paper comparison.
"""
    (out / "REPORT.md").write_text(md)
    print("wrote", out / "REPORT.md")
    if a.no_pdf:
        return
    env = dict(os.environ)
    if a.tex_bin:
        env["PATH"] = a.tex_bin + os.pathsep + env["PATH"]
    for tool in ("pandoc", "xelatex"):
        if shutil.which(tool, path=env["PATH"]) is None:
            sys.exit(f"{tool} not on PATH; PDF not built (Markdown and figures are written)")
    # PDF variant: no preview montage (the figures follow below) and no dead relative links
    pdf_md = "\n".join(l for l in md.splitlines()
                       if not l.startswith("[![Report preview]") and not l.startswith("*Click the image"))
    pdf_md = pdf_md.replace("[`docs/ALGORITHM.md`](../ALGORITHM.md)", "`docs/ALGORITHM.md`")
    pdf_md = pdf_md.replace("# mPing k-mer screen: preliminary results",
                            f"# mPing k-mer screen: preliminary results\n\n*Generated {datetime.date.today()} from code at commit {git_state()}.*", 1)
    (out / ".report_pdf_src.md").write_text(pdf_md)
    # keep each figure next to its text instead of floating to a later page
    (out / ".report_pdf_header.tex").write_text("\\usepackage{float}\n\\floatplacement{figure}{H}\n")
    subprocess.run(["pandoc", ".report_pdf_src.md", "-o", "REPORT.pdf", "--pdf-engine=xelatex",
                    "-V", "geometry:margin=2cm", "-V", "colorlinks=true", "-V", "fontsize=10pt",
                    "-H", ".report_pdf_header.tex", "--resource-path", "."], cwd=out, check=True, env=env)
    (out / ".report_pdf_src.md").unlink()
    (out / ".report_pdf_header.tex").unlink()
    print("wrote", out / "REPORT.pdf")


if __name__ == "__main__":
    main()
