#!/usr/bin/env python3
"""TE k-mer copy-number screen against per-strain jellyfish indexes.

Subcommands:
  kmers     canonical k-mers of the TE fasta, with in-file multiplicity
  query     depth of those k-mers in one jellyfish index
  baseline  single-copy depth (modal peak) from a jellyfish histogram
  analyze   normalise, compare to a reference strain, flag, plot
  gcbias    depth-vs-GC curve of one index, from single-copy-like k-mers
  profile   windowed copy number along each TE from an analyze output, with a breadth summary
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def canonical(s):
    r = revcomp(s)
    return s if s <= r else r


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def te_kmer_table(fasta, k):
    """One row per k-mer position; multiplicity counts canonical k-mers file-wide."""
    rows = []
    for te, seq in read_fasta(fasta):
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if set(kmer) <= set("ACGT"):
                rows.append((te, i, canonical(kmer)))
    df = pd.DataFrame(rows, columns=["te", "pos", "kmer"])
    df["multiplicity"] = df["kmer"].map(Counter(df["kmer"]))
    return df


def query_depths(index, kmers, jellyfish="jellyfish"):
    """Depth of each canonical k-mer in a jellyfish index (0 if absent)."""
    uniq = sorted(set(kmers))
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for i, km in enumerate(uniq):
            fh.write(f">{i}\n{km}\n")
        qpath = fh.name
    try:
        out = subprocess.run([jellyfish, "query", str(index), "-s", qpath],
                             check=True, capture_output=True, text=True).stdout
    finally:
        Path(qpath).unlink()
    depth = {}
    for line in out.splitlines():
        km, n = line.split()
        depth[canonical(km)] = int(n)
    return depth


def read_histo(path):
    arr = np.loadtxt(path, dtype=np.int64, ndmin=2)
    return arr[:, 0], arr[:, 1]


def baseline_peak(counts, freqs, max_valley=200, min_rise=1.1):
    """Single-copy depth: modal peak above the error valley of a k-mer histogram.

    The last bin of `jellyfish histo -h N` holds every k-mer with a higher count, so it is
    ignored. The error valley is searched only at counts <= max_valley. Returns NaN when there
    is no valley in that range, or when the peak above it is less than min_rise times the valley
    height (shallow data, where errors and the single-copy peak merge).
    """
    counts = np.asarray(counts)
    freqs = np.asarray(freqs, dtype=float)
    keep = counts < counts.max()
    top = int(counts[keep].max())
    h = np.zeros(top + 2)
    h[counts[keep]] = freqs[keep]
    sm = np.convolve(h, np.ones(3) / 3, mode="same")[:top + 1]
    limit = min(max_valley, top - 1)
    valley = 2
    while valley < limit and sm[valley] > sm[valley + 1]:
        valley += 1
    if valley >= limit:
        print(f"WARNING: no error valley at counts <= {max_valley}; baseline not estimated",
              file=sys.stderr)
        return float("nan")
    peak = valley + int(np.argmax(sm[valley:]))
    if sm[peak] < min_rise * sm[valley]:
        print(f"WARNING: no clear single-copy peak (valley={valley}, peak={peak}); baseline not estimated",
              file=sys.stderr)
        return float("nan")
    if peak <= valley + 1:
        print(f"WARNING: single-copy peak is next to the error valley (valley={valley}, peak={peak})",
              file=sys.stderr)
    return peak


def block_bootstrap_fc(a, b, block, n_boot=2000, seed=1):
    """CI for median(a)/median(b) over k-mer positions, resampling blocks of
    adjacent k-mers because overlapping k-mers are not independent."""
    rng = np.random.default_rng(seed)
    n = len(a)
    nb = max(1, int(np.ceil(n / block)))
    starts_max = max(1, n - block + 1)
    fcs = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate([np.arange(s, min(s + block, n))
                              for s in rng.integers(0, starts_max, nb)])[:n]
        mb = np.median(b[idx])
        fcs[i] = np.median(a[idx]) / mb if mb > 0 else np.nan
    return np.nanpercentile(fcs, [2.5, 97.5])


def gc_count(kmer):
    return kmer.count("G") + kmer.count("C")


def gc_bias_table(index, k, baseline, lo=0.5, hi=1.5, jellyfish="jellyfish"):
    """Relative depth by G+C count, from k-mers whose count is within [lo, hi] x the
    single-copy depth (mostly single-copy genome, so depth differences between GC
    bins are sequencing bias). bias = bin mean depth / mean depth over all bins.
    Streams jellyfish dump through awk because the index can hold ~4e8 k-mers."""
    if not np.isfinite(baseline):
        raise ValueError("no baseline for this index; cannot select single-copy-like k-mers")
    lo_c, hi_c = max(1, int(np.floor(lo * baseline))), int(np.ceil(hi * baseline))
    awk = ('{g=gsub(/[GC]/,"",$1); n[g]++; s[g]+=$2} '
           'END{for(g in n) print g"\t"n[g]"\t"s[g]}')
    dump = subprocess.Popen([jellyfish, "dump", "-c", "-L", str(lo_c), "-U", str(hi_c), str(index)],
                            stdout=subprocess.PIPE, text=True)
    agg = subprocess.run(["awk", awk], stdin=dump.stdout, capture_output=True, text=True, check=True)
    dump.stdout.close()
    if dump.wait() != 0:
        raise RuntimeError("jellyfish dump failed")
    rows = [(int(g), int(n), float(t)) for g, n, t in
            (line.split("\t") for line in agg.stdout.splitlines())]
    df = pd.DataFrame(sorted(rows), columns=["gc", "n_kmers", "sum_depth"])
    df["mean_depth"] = df["sum_depth"] / df["n_kmers"]
    overall = df["sum_depth"].sum() / df["n_kmers"].sum()
    df["bias"] = df["mean_depth"] / overall
    return df[["gc", "n_kmers", "mean_depth", "bias"]]


def correct_gc(depths, bias_by_strain, min_kmers=1000):
    """Divide each TE k-mer depth by its strain's bias for that k-mer's G+C count.
    GC bins with fewer than min_kmers genome k-mers get bias 1 (no correction)."""
    out = depths.copy()
    gc = out["kmer"].map(gc_count)
    for strain, tab in bias_by_strain.items():
        col = f"depth_{strain}"
        if col not in out:
            continue
        ok = tab[tab["n_kmers"] >= min_kmers].set_index("gc")["bias"]
        out[col] = out[col] / gc.map(ok).fillna(1.0)
    return out


def window_profile(depths, base, window, min_frac=0.5):
    """Windowed copy number along each TE from a te_kmer_depths table.

    Only k-mers unique within the TE file are used. Windows are runs of `window`
    consecutive unique k-mers; the last short run is kept if it has >= window/2 k-mers.
    Returns (windows, summary). `breadth` is the fraction of windows whose copy
    number is >= min_frac * the strain's median window copy number, so partial-length
    copies pull it below 1. `min_window_ratio` is the lowest window / median.
    """
    uniq = depths[depths["multiplicity"] == 1]
    win_rows, sum_rows = [], []
    for te, g in uniq.groupby("te", sort=False):
        g = g.sort_values("pos").reset_index(drop=True)
        edges = list(range(0, len(g), window))
        if len(edges) > 1 and len(g) - edges[-1] < window / 2:
            edges.pop()
        bounds = list(zip(edges, edges[1:] + [len(g)]))
        for s, b in base.items():
            col = f"depth_{s}"
            if col not in g:
                continue
            cn = np.array([np.median(g[col].iloc[lo:hi]) / b for lo, hi in bounds])
            for (lo, hi), c in zip(bounds, cn):
                win_rows.append((te, s, int(g["pos"].iloc[lo]), int(g["pos"].iloc[hi - 1]), hi - lo, c))
            med = float(np.median(cn))
            if med > 0:
                breadth = float(np.mean(cn >= min_frac * med))
                low = float(cn.min() / med)
            else:
                breadth = low = np.nan
            sum_rows.append((te, s, len(g), len(bounds), med, breadth, low))
    windows = pd.DataFrame(win_rows, columns=["te", "strain", "start", "end", "n_kmers", "cn"])
    summary = pd.DataFrame(sum_rows, columns=["te", "strain", "n_kmers_unique", "n_windows",
                                              "median_cn", "breadth", "min_window_ratio"])
    return windows, summary


def profile(args):
    res = Path(args.results)
    depths = pd.read_csv(res / f"te_kmer_depths.k{args.k}.tsv", sep="\t")
    base = pd.read_csv(res / f"baseline.k{args.k}.tsv", sep="\t", index_col=0)["single_copy_depth"]
    if args.gc_bias_dir:
        tabs = {s_: pd.read_csv(Path(args.gc_bias_dir) / f"{s_}.k{args.k}.gcbias.tsv", sep="\t")
                for s_ in base.index
                if (Path(args.gc_bias_dir) / f"{s_}.k{args.k}.gcbias.tsv").exists()}
        depths = correct_gc(depths, tabs)
        print(f"GC-corrected {len(tabs)} of {len(base)} strains", file=sys.stderr)
    windows, summary = window_profile(depths, base.to_dict(), args.window, args.min_frac)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    windows.to_csv(out / f"te_windows.k{args.k}.tsv", sep="\t", index=False)
    summary.to_csv(out / f"te_breadth.k{args.k}.tsv", sep="\t", index=False)
    if args.te:
        print(summary[summary.te == args.te].round(3).to_string(index=False))


def analyze(args):
    from scipy.stats import wilcoxon
    samples = pd.read_csv(args.samples, sep="\t")
    strains = list(samples["strain"])
    if args.reference not in strains:
        sys.exit(f"reference {args.reference} not in sample sheet")
    work = Path(args.workdir)
    kt = te_kmer_table(args.fasta, args.k)
    uniq = kt[kt["multiplicity"] == 1].copy()
    base, gb = {}, {}
    for s in strains:
        c, f = read_histo(work / f"{s}.k{args.k}.histo")
        base[s] = baseline_peak(c, f)
        if np.isnan(base[s]):
            print(f"WARNING: {s}: no baseline, copy numbers left empty", file=sys.stderr)
        depth = pd.read_csv(work / f"{s}.k{args.k}.depth.tsv", sep="\t",
                            index_col="kmer")["depth"]
        kt[f"depth_{s}"] = kt["kmer"].map(depth)
        uniq[f"depth_{s}"] = uniq["kmer"].map(depth)
        uniq[f"cn_{s}"] = uniq[f"depth_{s}"] / base[s]
    if np.isnan(base[args.reference]):
        sys.exit(f"reference {args.reference} has no baseline; choose another reference")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    pd.Series(base, name="single_copy_depth").to_csv(out / f"baseline.k{args.k}.tsv", sep="\t")
    kt.to_csv(out / f"te_kmer_depths.k{args.k}.tsv", sep="\t", index=False)

    rows = []
    for te, g in uniq.groupby("te", sort=False):
        g = g.sort_values("pos")
        row = {"te": te, "n_kmers_total": int((kt.te == te).sum()),
               "n_kmers_unique": len(g)}
        for s in strains:
            row[f"cn_{s}"] = float(np.median(g[f"cn_{s}"])) if len(g) else np.nan
        ref = g[f"cn_{args.reference}"].to_numpy()
        for s in strains:
            if s == args.reference:
                continue
            x = g[f"cn_{s}"].to_numpy()
            if np.isnan(base[s]):
                for key in ("fc", "wilcoxon_p", "fc_lo", "fc_hi"):
                    row[f"{key}_{s}"] = np.nan
                row[f"spike_{s}"] = False
                continue
            mref = np.median(ref) if len(ref) else np.nan
            fc = np.median(x) / mref if mref > 0 else np.nan
            row[f"fc_{s}"] = fc
            if len(g) >= 2 and not np.all(x == ref):
                row[f"wilcoxon_p_{s}"] = wilcoxon(x, ref).pvalue
            else:
                row[f"wilcoxon_p_{s}"] = np.nan
            if len(g) >= 2:
                lo, hi = block_bootstrap_fc(x, ref, args.k)
            else:
                lo = hi = np.nan
            row[f"fc_lo_{s}"], row[f"fc_hi_{s}"] = lo, hi
            row[f"spike_{s}"] = bool(fc >= args.threshold and lo > 1)
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(out / f"te_copy_number.k{args.k}.tsv", sep="\t", index=False)
    print(res.round(3).to_string(index=False))
    plot_tracks(uniq, strains, args.k, out / f"te_tracks.k{args.k}.png")


def plot_tracks(uniq, strains, k, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tes = list(dict.fromkeys(uniq["te"]))
    ncol = 3
    nrow = int(np.ceil(len(tes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 2.6 * nrow), squeeze=False)
    for ax, te in zip(axes.flat, tes):
        g = uniq[uniq.te == te].sort_values("pos")
        for s in strains:
            y = g[f"cn_{s}"].rolling(15, center=True, min_periods=1).median()
            ax.plot(g["pos"], y, lw=1, label=s)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_title(te, fontsize=9)
        ax.set_xlabel("position (bp)", fontsize=8)
        ax.set_ylabel("copy number", fontsize=8)
    for ax in list(axes.flat)[len(tes):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=7)
    fig.suptitle(f"Unique TE k-mer depth / single-copy depth (k={k}, 15-mer rolling median)")
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)

    a = sp.add_parser("kmers")
    a.add_argument("--fasta", required=True)
    a.add_argument("--k", type=int, required=True)
    a.add_argument("--out", required=True)

    a = sp.add_parser("query")
    a.add_argument("--index", required=True)
    a.add_argument("--fasta", required=True)
    a.add_argument("--k", type=int, required=True)
    a.add_argument("--out", required=True)

    a = sp.add_parser("baseline")
    a.add_argument("--histo", required=True)

    a = sp.add_parser("gcbias")
    a.add_argument("--index", required=True)
    a.add_argument("--k", type=int, required=True)
    a.add_argument("--histo", required=True, help="histogram used to get the single-copy depth")
    a.add_argument("--out", required=True)

    a = sp.add_parser("analyze")
    a.add_argument("--fasta", required=True)
    a.add_argument("--samples", required=True)
    a.add_argument("--workdir", required=True)
    a.add_argument("--k", type=int, required=True)
    a.add_argument("--reference", default="NB")
    a.add_argument("--threshold", type=float, default=2.0)
    a.add_argument("--out", required=True)

    a = sp.add_parser("profile")
    a.add_argument("--results", required=True, help="analyze output directory")
    a.add_argument("--k", type=int, required=True)
    a.add_argument("--window", type=int, default=40, help="unique k-mers per window")
    a.add_argument("--min-frac", type=float, default=0.5)
    a.add_argument("--te", help="print the summary for this TE")
    a.add_argument("--gc-bias-dir", help="directory of <strain>.k<k>.gcbias.tsv; corrects depths if given")
    a.add_argument("--out", required=True)

    args = p.parse_args()
    if args.cmd == "kmers":
        te_kmer_table(args.fasta, args.k).to_csv(args.out, sep="\t", index=False)
    elif args.cmd == "query":
        kt = te_kmer_table(args.fasta, args.k)
        d = query_depths(args.index, kt["kmer"])
        pd.DataFrame({"kmer": list(d), "depth": list(d.values())}).to_csv(
            args.out, sep="\t", index=False)
    elif args.cmd == "baseline":
        print(baseline_peak(*read_histo(args.histo)))
    elif args.cmd == "gcbias":
        b = baseline_peak(*read_histo(args.histo))
        gc_bias_table(args.index, args.k, b).to_csv(args.out, sep="\t", index=False)
    elif args.cmd == "profile":
        profile(args)
    else:
        analyze(args)


if __name__ == "__main__":
    main()
