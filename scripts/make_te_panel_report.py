#!/usr/bin/env python3
"""Report: GC correction, k=21 against k=31, and all 14 TEs in the library.

Reads two analyze output directories (k=21 and k=31) and the GC-bias tables in work/gcbias/.
For each k it computes median copy numbers with and without GC correction
(`mitekmer.copy_number_table`). It then writes:
  <out>/REPORT.md, <out>/REPORT.pdf (when pandoc and xelatex are available)
  <out>/figures/*.png
  <out>/mping_metrics.tsv      mPing against the paper, per k and correction
  <out>/te_panel.tsv           per TE and k: median, spread, depth dependence, k21 vs k31 agreement
  <out>/te_copy_numbers_k*.tsv per-strain copy numbers, uncorrected and corrected

usage: make_te_panel_report.py --k21 results/success/k21 --k31 results/success/k31 \
         --gc-dir work/gcbias --out docs/report/te_panel
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

BINS = [(0, 0, "0"), (1, 5, "1-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 50, "21-50"),
        (51, 100, "51-100"), (101, 10**6, ">100")]


def load_paper():
    """Paper mPing copy number per strain: the strain-group table plus the stratified 3K list."""
    grp = pd.read_csv(ROOT / "resources/strain_groups.tsv", sep="\t").set_index("strain")
    stream = pd.read_csv(ROOT / "resources/Chen2019/stream100.tsv", sep="\t").set_index("strain")
    rows = [(s, grp.loc[s, "group"], grp.loc[s, "paper_mping"]) for s in grp.index]
    rows += [(s, "3K", stream.loc[s, "paper_mping"]) for s in stream.index if s not in grp.index]
    return pd.DataFrame(rows, columns=["strain", "group", "paper"]).drop_duplicates("strain").set_index("strain")

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def load_k(results, k, gc_dir):
    depths = pd.read_csv(Path(results) / f"te_kmer_depths.k{k}.tsv", sep="\t")
    base = pd.read_csv(Path(results) / f"baseline.k{k}.tsv", sep="\t", index_col=0)["single_copy_depth"]
    base = base.dropna()
    gc = mk.load_gc_tables(gc_dir, k, base.index)
    raw = mk.copy_number_table(depths, base.to_dict())
    corr = mk.copy_number_table(depths, base.to_dict(), gc)
    return dict(depths=depths, base=base, raw=raw, corr=corr, gc=gc)


def mping_metrics(cn, paper, label):
    s = [x for x in cn.columns if x in paper.index and pd.notna(paper.loc[x, "paper"])]
    x = paper.loc[s, "paper"].astype(float).to_numpy()
    y = cn.loc["mPing", s].astype(float).to_numpy()
    ok = np.isfinite(y)
    x, y = x[ok], y[ok]
    pos = x > 0
    rec = {"set": label, "n": len(x),
           "pearson": pearsonr(x, y)[0] if len(x) > 2 else np.nan,
           "spearman": spearmanr(x, y)[0] if len(x) > 2 else np.nan,
           "median_ratio_pos": np.median(y[pos] / x[pos]) if pos.any() else np.nan,
           "median_abs_err": np.median(np.abs(y - x)),
           "max_abs_err": np.max(np.abs(y - x)) if len(x) else np.nan}
    return rec


def te_summary(raw, corr, base, k):
    rows = []
    for te in corr.index:
        v = corr.loc[te].dropna()
        v = v[v > 0]
        common = corr.loc[te].index.intersection(base.index)
        d = pd.DataFrame({"cn": corr.loc[te, common], "b": base[common]}).dropna()
        rho_depth = spearmanr(d.cn, d.b)[0] if len(d) > 2 else np.nan
        rows.append({"te": te, "k": k, "n_strains": int(len(v)),
                     "median_cn": float(v.median()) if len(v) else np.nan,
                     "q10": float(v.quantile(0.1)) if len(v) else np.nan,
                     "q90": float(v.quantile(0.9)) if len(v) else np.nan,
                     "log10_sd": float(np.log10(v).std()) if len(v) > 1 else np.nan,
                     "rho_cn_vs_depth": rho_depth,
                     "gc_corrected_strains": int(sum(1 for s in common if s in raw.columns)),
                     })
    return pd.DataFrame(rows)


def fig_gc_curves(gc_tables, path):
    fig, ax = plt.subplots(figsize=(6.4, 4))
    for name, tabs in gc_tables.items():
        for s, t in tabs.items():
            ax.plot(t.gc, t.bias, lw=0.8, alpha=0.5, color="#2a7fb8" if name == "k21" else "#c0392b")
    ax.axhline(1, color="#666", lw=0.8, ls="--")
    ax.set_xlabel("G+C bases in the k-mer")
    ax.set_ylabel("relative depth (bias)")
    ax.set_title("Depth by G+C content, one line per strain", fontsize=9)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color="#2a7fb8", label="k=21"), Line2D([], [], color="#c0392b", label="k=31")],
              frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_mping(paper, panels, path):
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 4.4), squeeze=False)
    for ax, (title, cn) in zip(axes[0], panels):
        s = [x for x in cn.columns if x in paper.index and pd.notna(paper.loc[x, "paper"])]
        x = paper.loc[s, "paper"].astype(float)
        y = cn.loc["mPing", s].astype(float)
        m = y.notna()
        ax.scatter(x[m], y[m], s=18, color="#2a7fb8", alpha=0.8)
        ax.plot([0, 190], [0, 190], color="#999", lw=0.8, ls="--")
        ax.axhspan(231, 503, color="#c0392b", alpha=0.07, lw=0)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("paper mPing (Chen 2019)")
        ax.set_ylabel("k-mer mPing")
        ax.set_xlim(-3, 190)
        ax.set_ylim(-3, 260)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_heatmap(cn, path):
    order = cn.loc["mPing"].sort_values(ascending=False).index
    M = np.log10(cn[order].clip(lower=0.1))
    fig, ax = plt.subplots(figsize=(9, 4.8))
    im = ax.imshow(M.values, aspect="auto", cmap="viridis", vmin=-1, vmax=2.7)
    ax.set_yticks(range(len(M.index)))
    ax.set_yticklabels(M.index, fontsize=7)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=70, ha="right", fontsize=5.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("log10 copy number (GC-corrected)", fontsize=8)
    ax.set_title("Copy number of each TE in each strain", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_te_spread(corr, path):
    tes = list(corr.index)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for i, te in enumerate(tes):
        v = corr.loc[te].dropna()
        v = v[v > 0]
        jit = np.random.default_rng(i).uniform(-0.18, 0.18, len(v))
        ax.scatter(np.full(len(v), i) + jit, v, s=8, alpha=0.6, color="#2a7fb8" if te != "mPing" else "#c0392b")
    ax.set_yscale("log")
    ax.set_xticks(range(len(tes)))
    ax.set_xticklabels(tes, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("copy number (GC-corrected, k=21)")
    ax.set_title("Spread of copy number across strains, per TE", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def mping_window_stats(d, min_cn=8, window=40):
    """Relative copy number per mPing window, without and with GC correction.
    Only strains with a median mPing copy number of at least min_cn are used, so counting noise is small."""
    out = {}
    for name, depths in (("raw", d["depths"]),
                         ("corrected", mk.correct_gc(d["depths"], d["gc"]))):
        win, _ = mk.window_profile(depths, d["base"].to_dict(), window)
        w = win[win.te == "mPing"].copy()
        w["rel"] = w.cn / w.groupby("strain").cn.transform("median")
        keep = w.groupby("strain").cn.median()
        keep = keep[keep >= min_cn].index
        out[name] = w[w.strain.isin(keep)]
    return out


def fig_mping_windows(stats, seq_gc, path):
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(6.6, 5), sharex=True, gridspec_kw={"height_ratios": [1, 3]})
    starts = sorted(stats["raw"].start.unique())
    a0.plot(range(len(starts)), seq_gc, marker="o", color="#555", ms=3)
    a0.set_ylabel("GC")
    for name, color in (("raw", "#7f8c8d"), ("corrected", "#2a7fb8")):
        w = stats[name]
        g = w.groupby("start").rel
        med, q10, q90 = g.median(), g.quantile(0.1), g.quantile(0.9)
        x = range(len(med))
        a1.plot(x, med.values, color=color, label=f"{name} (median)", marker="o", ms=3)
        a1.fill_between(x, q10.values, q90.values, color=color, alpha=0.15, lw=0)
    a1.axhline(1, color="#666", lw=0.8, ls="--")
    a1.set_xticks(range(len(starts)))
    a1.set_xticklabels([str(x) for x in starts], rotation=45, fontsize=7)
    a1.set_xlabel("mPing window start (bp)")
    a1.set_ylabel("copy number / strain median")
    a1.legend(frameon=False, fontsize=8)
    a1.set_title("Band: 10th to 90th percentile across strains", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_bins(paper, cn, path):
    s = [x for x in cn.columns if x in paper.index and pd.notna(paper.loc[x, "paper"])]
    d = pd.DataFrame({"paper": paper.loc[s, "paper"].astype(float), "est": cn.loc["mPing", s].astype(float)}).dropna()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.axhspan(231, 503, color="#c0392b", alpha=0.08, lw=0)
    rows, labels = [], []
    for i, (lo, hi, lab) in enumerate([b for b in BINS if ((d.paper >= b[0]) & (d.paper <= b[1])).any()]):
        v = d[(d.paper >= lo) & (d.paper <= hi)]
        ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.18, 0.18, len(v)), v.est,
                   s=14, color="#2a7fb8", alpha=0.7, zorder=3)
        ax.plot([i - 0.25, i + 0.25], [v.est.median()] * 2, color="#222", lw=2, zorder=4)
        pos = v[v.paper > 0]
        rows.append({"paper bin": lab, "strains": len(v), "median estimate": v.est.median(),
                     "median estimate / paper": (pos.est / pos.paper).median() if len(pos) else np.nan})
        labels.append(f"{lab}\n(n={len(v)})")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylim(0, 600)
    ax.set_xlabel("paper mPing copies (bin)")
    ax.set_ylabel("k-mer mPing copies (k=21, GC-corrected)")
    ax.set_title("Bar: median. Shaded band: paper's burst range (231 to 503)", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return pd.DataFrame(rows)


def data_status(k_results, strain_list, paper):
    """Which strains of the stratified list have data, no baseline, or no output."""
    base = pd.read_csv(Path(k_results) / "baseline.k21.tsv", sep="\t", index_col=0).iloc[:, 0]
    ok = [x for x in strain_list if x in base.index and pd.notna(base[x])]
    nobase = [x for x in strain_list if x in base.index and pd.isna(base[x])]
    missing = [x for x in strain_list if x not in base.index]
    return ok, nobase, missing


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join("" if pd.isna(v) else (f"{v:.3g}" if isinstance(v, float) else str(v))
                                       for v in r) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k21", required=True)
    ap.add_argument("--k31", required=True)
    ap.add_argument("--gc-dir", default=str(ROOT / "work/gcbias"))
    ap.add_argument("--out", default=str(ROOT / "docs/report"))
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--tex-bin", default=os.environ.get("TEX_BIN", ""),
                    help="directory with a TeX Live xelatex to put first on PATH (env TEX_BIN)")
    a = ap.parse_args()
    out, figs = Path(a.out), Path(a.out) / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    tabs_dir = out / "tables"
    tabs_dir.mkdir(exist_ok=True)
    paper = load_paper()

    d21 = load_k(a.k21, 21, a.gc_dir)
    d31 = load_k(a.k31, 31, a.gc_dir)

    metrics = []
    for k, d in ((21, d21), (31, d31)):
        metrics.append(mping_metrics(d["raw"], paper, f"k{k} uncorrected, all strains"))
        metrics.append(mping_metrics(d["corr"], paper, f"k{k} GC-corrected, all strains"))
    common = [s for s in d21["raw"].columns if s in d31["raw"].columns]
    for k, d in ((21, d21), (31, d31)):
        metrics.append(mping_metrics(d["raw"][common], paper, f"k{k} uncorrected, shared strains"))
        metrics.append(mping_metrics(d["corr"][common], paper, f"k{k} GC-corrected, shared strains"))
    met = pd.DataFrame(metrics)
    met.to_csv(tabs_dir / "mping_metrics.tsv", sep="\t", index=False)

    t21 = te_summary(d21["raw"], d21["corr"], d21["base"], 21)
    t31 = te_summary(d31["raw"], d31["corr"], d31["base"], 31)
    panel = t21.merge(t31, on="te", suffixes=("_k21", "_k31"))
    rho_k = []
    for te in d21["corr"].index:
        s = [x for x in common if pd.notna(d21["corr"].loc[te, x]) and pd.notna(d31["corr"].loc[te, x])]
        rho_k.append(spearmanr(d21["corr"].loc[te, s], d31["corr"].loc[te, s])[0] if len(s) > 2 else np.nan)
    panel["rho_k21_vs_k31"] = rho_k
    panel["q90_over_q10_k21"] = panel.q90_k21 / panel.q10_k21
    panel["median_ratio_k31_over_k21"] = panel.median_cn_k31 / panel.median_cn_k21
    panel.to_csv(tabs_dir / "te_panel.tsv", sep="\t", index=False)

    d21["corr"].round(3).to_csv(tabs_dir / "te_copy_numbers_k21.tsv", sep="\t")
    d31["corr"].round(3).to_csv(tabs_dir / "te_copy_numbers_k31.tsv", sep="\t")

    wstats = mping_window_stats(d21)
    mp_seq = dict(mk.read_fasta(ROOT / "RiceTE.fa"))["mPing"]
    st = sorted(wstats["raw"].start.unique())
    en = wstats["raw"].groupby("start").end.first()
    seq_gc = [mk.gc_count(mp_seq[a:en[a] + 21]) / len(mp_seq[a:en[a] + 21]) for a in st]
    fig_mping_windows(wstats, seq_gc, figs / "mping_windows.png")
    low_win = min(st, key=lambda a: wstats["raw"][wstats["raw"].start == a].rel.median())
    def wq(name):
        v = wstats[name][wstats[name].start == low_win].rel
        return v.median(), v.quantile(0.1), v.nunique() and len(v)
    r_med, r_q10, n_strains = wq("raw")
    c_med, c_q10, _ = wq("corrected")
    dip_text = (
        f"The lowest window is the one starting at base {low_win}. Across the {n_strains} strains with at least 8 mPing copies, "
        f"its copy number is a median {r_med:.2f} of the strain median without correction and {c_med:.2f} with correction. "
        f"The 10th percentile strain is at {r_q10:.2f} without and {c_q10:.2f} with correction. "
        f"GC correction therefore removes part of the median dip, but the lowest-copy strains stay low. "
        f"The cause of that remaining dip is not known from these data."
    )
    fig_gc_curves({"k21": d21["gc"], "k31": d31["gc"]}, figs / "gc_curves.png")
    fig_mping(paper, [("k=21, uncorrected", d21["raw"]), ("k=21, GC-corrected", d21["corr"]),
                      ("k=31, uncorrected", d31["raw"]), ("k=31, GC-corrected", d31["corr"])],
              figs / "mping_vs_paper.png")
    fig_heatmap(d21["corr"], figs / "te_heatmap_k21.png")
    fig_te_spread(d21["corr"], figs / "te_spread_k21.png")

    bins_tab = fig_bins(paper, d21["corr"], figs / "mping_by_paper_bin.png")
    bins_tab.to_csv(tabs_dir / "mping_by_paper_bin.tsv", sep="\t", index=False)
    grp = pd.read_csv(ROOT / "resources/strain_groups.tsv", sep="\t").set_index("strain")
    ref_rows = []
    for st_ in grp.index[grp.group.isin(["burst", "reference", "landrace_unverified"])]:
        ref_rows.append({"strain": st_, "group": grp.loc[st_, "group"],
                         "paper": grp.loc[st_, "paper_mping"],
                         "k=21 corrected": d21["corr"].loc["mPing", st_] if st_ in d21["corr"].columns else np.nan,
                         "k=31 corrected": d31["corr"].loc["mPing", st_] if st_ in d31["corr"].columns else np.nan})
    ref_tab = pd.DataFrame(ref_rows)
    ref_tab.to_csv(tabs_dir / "reference_and_burst_strains.tsv", sep="\t", index=False)
    strat = pd.read_csv(ROOT / "resources/Chen2019/stream100.tsv", sep="\t").strain.tolist()
    ok21, nobase21, miss21 = data_status(a.k21, strat, paper)
    s31 = pd.read_csv(Path(a.k31) / "baseline.k31.tsv", sep="\t", index_col=0).iloc[:, 0]
    ok31 = [x for x in strat if x in s31.index and pd.notna(s31[x])]
    nobase31 = [x for x in strat if x in s31.index and pd.isna(s31[x])]
    miss31 = [x for x in strat if x not in s31.index]
    b191_note = (" B191: NCBI has no full-format .sra for run ERR624821 (only an SRA Lite file), so the job stopped "
                 "before counting. SRA Lite files have simplified quality scores." if "B191" in miss21 else "")
    b191_line = ("\n  - " + b191_note.strip()) if b191_note else ""
    nd21, nd31 = len(d21["gc"]), len(d31["gc"])
    n21, n31 = d21["raw"].shape[1], d31["raw"].shape[1]
    today = datetime.date.today().isoformat()
    m = met.copy()
    mp_tab = m.rename(columns={"set": "setting", "median_ratio_pos": "median ratio to paper (paper > 0)",
                               "median_abs_err": "median abs error", "max_abs_err": "max abs error"})
    te_tab = panel[["te", "median_cn_k21", "q90_over_q10_k21", "log10_sd_k21", "rho_cn_vs_depth_k21",
                    "median_cn_k31", "median_ratio_k31_over_k21", "rho_k21_vs_k31"]].round(2)
    te_tab.columns = ["TE", "median k=21", "q90/q10 k=21", "log10 SD k=21", "rho vs depth k=21",
                      "median k=31", "median k=31 / k=21", "rho k=21 vs k=31"]
    others = panel[~panel.te.isin(["mPing"])]
    varying = panel.sort_values("log10_sd_k21", ascending=False).head(2)
    steady = panel[~panel.te.isin(varying.te)]
    disagree = panel[panel.rho_k21_vs_k31 < 0.6].sort_values("rho_k21_vs_k31")
    depth_neg = panel[panel.rho_cn_vs_depth_k21 <= -0.3].sort_values("rho_cn_vs_depth_k21")
    findings = (
        f"- **Which TEs vary between strains.** {varying.te.iloc[0]} and {varying.te.iloc[1]} vary most "
        f"(log10 SD {varying.log10_sd_k21.iloc[0]:.2f} and {varying.log10_sd_k21.iloc[1]:.2f}). The other "
        f"{len(steady)} TEs have log10 SD from {steady.log10_sd_k21.min():.2f} to {steady.log10_sd_k21.max():.2f}. "
        f"For them, the 90th percentile strain is {steady.q90_over_q10_k21.min():.1f} to "
        f"{steady.q90_over_q10_k21.max():.1f} times the 10th percentile strain. "
        f"These TEs are recovered at similar copy numbers across strains. There is no published per-strain "
        f"value for them, so this shows consistency across strains, not accuracy.\n"
        f"- **k=21 against k=31.** For mPing the two agree (rank correlation {panel.set_index('te').loc['mPing','rho_k21_vs_k31']:.2f}, "
        f"median ratio {panel.set_index('te').loc['mPing','median_ratio_k31_over_k21']:.2f}). For {len(disagree)} TEs the rank "
        f"correlation is below 0.6: " + ", ".join(f"{r.te} ({r.rho_k21_vs_k31:.2f})" for r in disagree.itertuples()) + ". "
        f"For most TEs the k=31 median is lower than the k=21 median (ratio {panel.median_ratio_k31_over_k21.min():.2f} to "
        f"{panel.median_ratio_k31_over_k21.max():.2f}). Because there is no reference value for these TEs, these data do not say "
        f"which k is right for them. One explanation to test is that diverged copies lose exact 31-base matches more often "
        f"than 21-base matches. mPing copies are near-identical, which fits its agreement.\n"
        f"- **Dependence on depth.** At k=21 the copy number falls as single-copy depth rises for "
        + ", ".join(f"{r.te} ({r.rho_cn_vs_depth_k21:.2f})" for r in depth_neg.itertuples())
        + ". A rank correlation of this size can come from the baseline estimate or from real differences between strains. "
        "These data do not separate the two."
    )

    sh = met.set.str.contains("shared")
    s21 = met[sh & met.set.str.startswith("k21 GC")].iloc[0]
    s31 = met[sh & met.set.str.startswith("k31 GC")].iloc[0]
    n_uniq21 = int(((d21["depths"].te == "mPing") & (d21["depths"].multiplicity == 1)).sum())
    n_uniq31 = int(((d31["depths"].te == "mPing") & (d31["depths"].multiplicity == 1)).sum())
    better = "k=21" if (s21.median_abs_err < s31.median_abs_err) else "k=31"
    gap = abs(s21.median_abs_err - s31.median_abs_err)
    size = "small" if gap < 0.25 else "clear"
    k_choice = (
        f"On the {int(s21.n)} strains that have both k values and a paper value, the GC-corrected "
        f"k=21 estimate has Spearman {s21.spearman:.3f} and median absolute error {s21.median_abs_err:.2f} copies. "
        f"The GC-corrected k=31 estimate has Spearman {s31.spearman:.3f} and median absolute error "
        f"{s31.median_abs_err:.2f} copies. {better} has the smaller error, and the difference is {size} "
        f"({gap:.2f} copies). mPing has {n_uniq21} unique k-mers at k=21 and {n_uniq31} at k=31, "
        f"so specificity does not separate the two values of k for mPing."
    )
    s31_nobase = int(pd.read_csv(Path(a.k31) / "baseline.k31.tsv", sep="\t", index_col=0).iloc[:, 0].isna().sum())
    s21_nobase = int(pd.read_csv(Path(a.k21) / "baseline.k21.tsv", sep="\t", index_col=0).iloc[:, 0].isna().sum())
    k_choice_short = (f"On the {int(s21.n)} strains with both k values and a paper value, k=21 has median absolute error "
                      f"{s21.median_abs_err:.2f} copies and k=31 has {s31.median_abs_err:.2f}. Spearman is {s21.spearman:.3f} and "
                      f"{s31.spearman:.3f}. The difference is {size}. k=31 gives no baseline for more strains "
                      f"({s31_nobase} against {s21_nobase} at k=21).")
    r21c = met[met.set == "k21 GC-corrected, all strains"].iloc[0]
    r21u = met[met.set == "k21 uncorrected, all strains"].iloc[0]
    summary = (
        f"- mPing copy numbers from k-mer depth agree with the paper's read-depth values. For {int(r21c.n)} strains "
        f"at k=21 with GC correction, the rank correlation is {r21c.spearman:.3f} and the median ratio to the paper is "
        f"{r21c.median_ratio_pos:.2f}.\n"
        f"- GC correction changes the median absolute error from {r21u.median_abs_err:.2f} to {r21c.median_abs_err:.2f} copies "
        f"at k=21. The rank correlation changes from {r21u.spearman:.3f} to {r21c.spearman:.3f}.\n"
        f"- {k_choice_short}\n"
        f"- {varying.te.iloc[0]} and {varying.te.iloc[1]} vary most across strains. The other {len(steady)} TEs are recovered at "
        f"similar copy numbers in all strains. For {len(disagree)} TEs ({', '.join(disagree.te)}) k=21 and k=31 rank the "
        f"strains differently.\n"
        f"- The dip near mPing base {low_win} is only partly explained by GC bias.")
    md = f"""# TE copy numbers from k-mer depth: results

Generated {today}. This is the single current report. Earlier versions are in the git history.

## Summary

{summary}

## Data

- Stratified 3K set: {len(strat)} accessions.
  - k=21: {len(ok21)} have a copy number. {len(nobase21)} have no clear single-copy peak. No output: {len(miss21)} ({', '.join(miss21) if miss21 else 'none'}).
  - k=31: {len(ok31)} have a copy number. {len(nobase31)} have no clear single-copy peak. No output: {len(miss31)} ({', '.join(miss31) if miss31 else 'none'}).
{b191_line}
- All strains in the analysis: k=21 {n21}, k=31 {n31}, counting only those with a baseline. {nd21} (k=21) and {nd31} (k=31) have a GC-bias curve. Strains without a curve are not corrected.
- A strain with no clear single-copy peak (no baseline) has no values. It is missing from the figures, not set to zero.

## mPing across the paper's range

![k-mer mPing estimate by paper bin](figures/mping_by_paper_bin.png)

{md_table(bins_tab.round(2))}

The 3K accessions were chosen at random within paper bins, not at random from the 3,000 lines, so overall rates need weighting by bin.

## Reference and burst strains

{md_table(ref_tab.round(1))}

The paper gives a value only for Nipponbare. For the burst lines it gives only a range for the four bursting accessions (231 to 503). HEG4, A119 and A123 are three of them. A157 is a landrace expected to have high mPing, from the project notes. The paper does not list it.

## GC bias

![Depth against G+C content](figures/gc_curves.png)

Each line is one strain. A flat line at 1 means no bias. Lines that slope down at high G+C show lower depth
for G+C-rich k-mers. The correction divides each k-mer's depth by its strain's value for the k-mer's G+C count.
Count bins with fewer than 1,000 genome k-mers are not corrected.

## mPing along its sequence

![mPing copy number by window, with and without GC correction](figures/mping_windows.png)

Windows are 40 consecutive mPing k-mers (k=21). Only strains with at least 8 mPing copies are used. {dip_text}

## mPing against the paper

![mPing copy number against the paper, four settings](figures/mping_vs_paper.png)

{md_table(mp_tab)}

Each row uses the paper's read-depth value (Chen et al. 2019) as the reference.
"Shared strains" restricts both k values to the strains that have both.

## Which k is better for mPing

{k_choice}

## All TEs in the library

The library has 14 TEs. The paper gives values only for mPing (and for Ping and Pong, which are not in
`RiceTE.fa`). For the other 13 TEs there is no reference value, so their copy numbers can only be checked
for internal consistency:

- **Agreement between k values** (`rho k=21 vs k=31`). A TE with a stable copy number should agree
  between k=21 and k=31.
- **Dependence on depth** (`rho vs depth`). A copy number that rises or falls with the strain's single-copy
  depth may reflect an artifact.

{md_table(te_tab)}

### Findings

{findings}

![Copy number of each TE in each strain](figures/te_heatmap_k21.png)

![Spread of copy numbers across strains](figures/te_spread_k21.png)

## Limits

- The 100-accession set was sampled within paper bins, so these results are not a population estimate.
- Strains in the figures are the ones with both a baseline and a GC curve, unless the table says otherwise.
- The paper's values come from read depth. They are independent of k-mer depth, but they are only available for mPing.
- The correction assumes the bias is a function of G+C content alone. Other sequence features are not corrected.
"""
    (out / "REPORT.md").write_text(md)
    print("wrote", out / "REPORT.md")
    if a.no_pdf:
        return
    env = dict(os.environ)
    if a.tex_bin:
        env["PATH"] = a.tex_bin + os.pathsep + env["PATH"]
    for tool in ("pandoc", "xelatex"):
        if shutil.which(tool, path=env.get("PATH")) is None:
            sys.exit(f"{tool} not on PATH; PDF not built")
    subprocess.run(["pandoc", "REPORT.md", "-o", "REPORT.pdf", "--pdf-engine=xelatex",
                    "-V", "geometry:margin=2cm", "-V", "colorlinks=true", "-V", "fontsize=10pt",
                    "--resource-path", "."], cwd=out, check=True, env=env)
    print("wrote", out / "REPORT.pdf")


if __name__ == "__main__":
    main()
