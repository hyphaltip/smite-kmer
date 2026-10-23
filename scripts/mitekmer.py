#!/usr/bin/env python3
"""TE k-mer copy-number screen against per-strain jellyfish indexes.

Subcommands:
  kmers     canonical k-mers of the TE fasta, with in-file multiplicity
  query     depth of those k-mers in one jellyfish index
  baseline  single-copy depth (modal peak) from a jellyfish histogram
  analyze   normalise, compare to a reference strain, flag, plot
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


def baseline_peak(counts, freqs):
    """Modal peak above the error valley of a k-mer histogram."""
    h = np.zeros(int(counts.max()) + 2)
    h[counts] = freqs
    s = np.convolve(h, np.ones(3) / 3, mode="same")
    valley = 2
    while valley + 1 < len(s) and s[valley] > s[valley + 1]:
        valley += 1
    peak = valley + int(np.argmax(s[valley:]))
    if peak <= valley + 1:
        print(f"WARNING: no clear single-copy peak (valley={valley}, peak={peak})",
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
        depth = pd.read_csv(work / f"{s}.k{args.k}.depth.tsv", sep="\t",
                            index_col="kmer")["depth"]
        kt[f"depth_{s}"] = kt["kmer"].map(depth)
        uniq[f"depth_{s}"] = uniq["kmer"].map(depth)
        uniq[f"cn_{s}"] = uniq[f"depth_{s}"] / base[s]
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

    a = sp.add_parser("analyze")
    a.add_argument("--fasta", required=True)
    a.add_argument("--samples", required=True)
    a.add_argument("--workdir", required=True)
    a.add_argument("--k", type=int, required=True)
    a.add_argument("--reference", default="NB")
    a.add_argument("--threshold", type=float, default=2.0)
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
    else:
        analyze(args)


if __name__ == "__main__":
    main()
