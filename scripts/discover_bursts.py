#!/usr/bin/env python3
"""Find sequences that are high copy in one strain but not in a comparison set, without a TE library.

Steps:
  1. Dump the k-mers of the test strain whose count is at least MIN_CN times its single-copy depth.
  2. Look each one up in every comparison strain. Keep k-mers that are at least MIN_FOLD times higher
     than the comparison median (+1), so a family that is high everywhere is dropped.
  3. Build a de Bruijn graph from the kept k-mers. Drop edges to k-mers that are much rarer than their
     siblings (sequencing-error variants), then read off the unitigs.
  4. Write the unitigs with their copy number in the test strain and in each comparison strain, ranked by
     fold x length. Unitigs that are mostly k-mers of the --mask files (adapters and so on) are flagged.
  Bubbles (two paths that rejoin within 2k steps, as from a common second allele) are popped before the unitigs
  are read off, so an element is not split at a SNP. Every step-1 k-mer is saved in kmers.tsv.zst.
If --fasta is given, each unitig is also compared with that library by shared k-mers. This is a check
and does not affect which sequences are found.

usage: discover_bursts.py --positive HEG4 --cohort B031 CX55 ... --workdir work --k 21 --out DIR
Needs jellyfish on PATH. Indexes are work/<strain>.k<k>.jf, histograms work/<strain>.k<k>.histo.
"""
import argparse
import os
import subprocess
import time
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mitekmer as mk  # noqa: E402

BASES = "ACGT"


# ---------- graph ----------

class Graph:
    """De Bruijn graph over canonical k-mers (odd k, so no k-mer is its own reverse complement).

    counts: {canonical k-mer: count}. An edge x->y is dropped when y is rarer than ratio times the most
    abundant successor of x, or x is rarer than ratio times the most abundant predecessor of y.
    These are sequencing-error variants.
    """

    def __init__(self, counts, ratio=0.1):
        self.counts = counts
        self.ratio = ratio
        self.k = len(next(iter(counts))) if counts else 0
        if self.k % 2 == 0 and self.k:
            raise ValueError("k must be odd")

    @staticmethod
    def canon(x):
        r = mk.revcomp(x)
        return x if x <= r else r

    def cnt(self, x):
        return self.counts.get(self.canon(x), 0)

    def raw_succ(self, x):
        base = x[1:]
        return [base + b for b in BASES if self.cnt(base + b)]

    def succ(self, x):
        s = self.raw_succ(x)
        if len(s) < 2:
            return s
        top = max(self.cnt(y) for y in s)
        return [y for y in s if self.cnt(y) >= self.ratio * top]

    def pred(self, y):
        s = self.raw_succ(mk.revcomp(y))
        if len(s) >= 2:
            top = max(self.cnt(z) for z in s)
            s = [z for z in s if self.cnt(z) >= self.ratio * top]
        return [mk.revcomp(z) for z in s]

    def out_edges(self, x):
        return [y for y in self.succ(x) if x in self.pred(y)]

    def in_edges(self, y):
        return [x for x in self.pred(y) if y in self.succ(x)]

    def pop_bubbles(self, max_len=None):
        """Remove the weaker arm of every bubble: two paths from one k-mer that rejoin within max_len steps.

        A common second allele (for example 15% of reads) survives the edge-ratio filter and would
        otherwise split a unitig in two. The arm with the larger total count is kept. Returns the
        number of arms removed.
        """
        max_len = max_len or 2 * self.k
        removed = 0
        changed = True
        while changed:
            changed = False
            for c in list(self.counts):
                if c not in self.counts:
                    continue
                for x in (c, mk.revcomp(c)):
                    outs = self.out_edges(x)
                    if len(outs) < 2:
                        continue
                    arms = {}
                    for y in outs:
                        path, cur, merge = [y], y, None
                        while len(path) <= max_len:
                            if len(self.in_edges(cur)) > 1:
                                merge = cur
                                path.pop()
                                break
                            o = self.out_edges(cur)
                            if len(o) != 1:
                                break
                            cur = o[0]
                            path.append(cur)
                        if merge is not None:
                            arms.setdefault(merge, []).append(path)
                    for merge, paths in arms.items():
                        if len(paths) < 2:
                            continue
                        paths.sort(key=lambda p: sum(self.cnt(n) for n in p), reverse=True)
                        for loser in paths[1:]:
                            for n in loser:
                                if self.counts.pop(self.canon(n), None) is not None:
                                    changed = True
                            removed += 1
                    if changed:
                        break
        return removed

    def unitigs(self):
        """Maximal unitigs as (sequence, [k-mer counts along it])."""
        visited, out = set(), []
        for start in list(self.counts):
            if start in visited:
                continue
            path = [start]
            visited.add(start)
            x = start
            while True:
                o = self.out_edges(x)
                if len(o) != 1:
                    break
                y = o[0]
                if len(self.in_edges(y)) != 1 or self.canon(y) in visited:
                    break
                path.append(y)
                visited.add(self.canon(y))
                x = y
            x, back = mk.revcomp(start), []
            while True:
                o = self.out_edges(x)
                if len(o) != 1:
                    break
                y = o[0]
                if len(self.in_edges(y)) != 1 or self.canon(y) in visited:
                    break
                back.append(y)
                visited.add(self.canon(y))
                x = y
            if back:
                path = [mk.revcomp(b) for b in reversed(back)] + path
            seq = path[0] + "".join(p[-1] for p in path[1:])
            out.append((seq, [self.cnt(p) for p in path]))
        return out


def build_unitigs(counts, ratio=0.1, pop=True):
    """Unitigs of the graph of `counts` (dict is copied). With pop=True, bubbles are popped first."""
    if not counts:
        return []
    g = Graph(dict(counts), ratio)
    if pop:
        g.pop_bubbles()
    return g.unitigs()


# ---------- jellyfish ----------

def dump_high(index, min_count):
    """k-mers (canonical) and counts at or above min_count, streamed from jellyfish dump."""
    p = subprocess.Popen(["jellyfish", "dump", "-c", "-L", str(int(min_count)), str(index)],
                         stdout=subprocess.PIPE, text=True)
    kmers, counts = [], []
    for line in p.stdout:
        km, n = line.split()
        kmers.append(km)
        counts.append(int(n))
    if p.wait() != 0:
        raise RuntimeError(f"jellyfish dump failed for {index}")
    return kmers, np.array(counts, dtype=np.int64)


def query_counts(index, kmers, tmpdir):
    """Count of each k-mer in an index, in the order given. Checks that each returned k-mer is the one asked for."""
    q = Path(tempfile.mkstemp(suffix=".fa", dir=tmpdir)[1])
    with open(q, "w") as fh:
        for i, km in enumerate(kmers):
            fh.write(f">{i}\n{km}\n")
    try:
        out = subprocess.run(["jellyfish", "query", "-s", str(q), str(index)],
                             check=True, capture_output=True, text=True).stdout.split("\n")
    finally:
        q.unlink(missing_ok=True)
    out = [line for line in out if line]
    if len(out) != len(kmers):
        raise RuntimeError(f"query returned {len(out)} lines for {len(kmers)} k-mers ({index})")
    res = np.empty(len(kmers), dtype=np.int64)
    for i, line in enumerate(out):
        km, n = line.split()
        if km != kmers[i]:
            raise RuntimeError(f"query line {i} is {km}, expected {kmers[i]} ({index})")
        res[i] = int(n)
    return res


def single_copy_depth(workdir, strain, k):
    c, f = mk.read_histo(Path(workdir) / f"{strain}.k{k}.histo")
    return mk.baseline_peak(c, f)


def read_mask_kmers(paths, k):
    """Canonical k-mers of FASTA files (.gz allowed): adapters, PhiX, vectors, organelles, rDNA."""
    import gzip
    mask = set()
    for path in paths:
        opener = gzip.open if str(path).endswith(".gz") else open
        seq = []

        def flush():
            q = "".join(seq).upper()
            for i in range(len(q) - k + 1):
                km = q[i:i + k]
                if set(km) <= set("ACGT"):
                    mask.add(mk.canonical(km))
        with opener(path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    flush()
                    seq = []
                elif line:
                    seq.append(line)
        flush()
    return mask


def write_zstd_table(df, path):
    """Write a TSV, compressed with zstd when it is installed, otherwise gzip."""
    import shutil
    tmp = str(path) + ".tmp"
    df.to_csv(tmp, sep="\t", index=False)
    if shutil.which("zstd"):
        subprocess.run(["zstd", "-q", "-T0", "-f", "--rm", "-o", str(path) + ".zst", tmp], check=True)
    else:
        subprocess.run(["gzip", "-f", tmp], check=True)
        os.replace(tmp + ".gz", str(path) + ".gz")


# ---------- main ----------

def run(args):
    import pandas as pd
    if args.k % 2 == 0:
        sys.exit("k must be odd")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    work = Path(args.workdir)
    k = args.k
    b_pos = single_copy_depth(work, args.positive, k)
    if not np.isfinite(b_pos):
        sys.exit(f"{args.positive} has no baseline")
    if b_pos < args.min_depth and not args.allow_low_depth:
        sys.exit(f"{args.positive}: single-copy depth {b_pos:.0f} is below --min-depth {args.min_depth:g}. "
                 f"Copy numbers are unreliable at this depth. Use --allow-low-depth to run anyway.")
    t0 = time.perf_counter()

    def lap(msg):
        print(f"[{time.perf_counter() - t0:6.1f}s] {msg}", file=sys.stderr)

    print(f"{args.positive}: single-copy depth {b_pos:.0f}", file=sys.stderr)
    kmers, cnt = dump_high(work / f"{args.positive}.k{k}.jf", args.min_cn * b_pos)
    cn_pos = cnt / b_pos
    lap(f"step 1: {len(kmers):,} k-mers with copy number >= {args.min_cn}")

    cohort, depths_c = [], []
    for s_ in args.cohort:
        if s_ == args.positive:
            continue
        b = single_copy_depth(work, s_, k)
        if not np.isfinite(b):
            print(f"WARNING: {s_} has no baseline, left out of the comparison set", file=sys.stderr)
            continue
        cohort.append(s_)
        depths_c.append(b)
    cohort_cn = []
    if cohort:
        from concurrent.futures import ThreadPoolExecutor
        with tempfile.TemporaryDirectory(dir=args.tmpdir or None) as tmp:
            def look(s_):
                return query_counts(work / f"{s_}.k{k}.jf", kmers, tmp)
            with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
                for s_, b, counts_c in zip(cohort, depths_c, pool.map(look, cohort)):
                    cohort_cn.append(counts_c / b)
                    lap(f"  looked up {len(kmers):,} k-mers in {s_}")
    if not cohort:
        sys.exit("no usable comparison strains")
    M = np.vstack(cohort_cn)
    med = np.median(M, axis=0)
    extend_fold = args.extend_fold if args.extend_fold else args.min_fold
    if extend_fold > args.min_fold:
        sys.exit("--extend-fold must not be larger than --min-fold")
    fold_k = cn_pos / (med + 1.0)
    seed = fold_k >= args.min_fold
    keep = fold_k >= extend_fold          # k-mers allowed in the graph; seeds are the strict subset
    lap(f"step 2: {int(seed.sum()):,} seed k-mers (fold >= {args.min_fold}), {int(keep.sum()):,} in the graph "
        f"(fold >= {extend_fold}), of {len(kmers):,}; {len(cohort)} comparison strains")

    # every step-1 k-mer with its copy numbers, so thresholds can be changed without new lookups
    tab_k = pd.DataFrame({"kmer": kmers, "count": cnt, "cn_test": np.round(cn_pos, 2),
                          "cn_comparison_median": np.round(med, 2), "fold": np.round(fold_k, 2),
                          "seed": seed.astype(int), "kept": keep.astype(int)})
    for s, row in zip(cohort, M):
        tab_k[f"cn_{s}"] = np.round(row, 2)
    write_zstd_table(tab_k, out / "kmers.tsv")

    idx = np.flatnonzero(keep)
    pos_in = {kmers[i]: j for j, i in enumerate(idx)}
    kept_cn, kept_M, kept_fold, kept_seed = cn_pos[idx], M[:, idx], fold_k[idx], seed[idx]
    counts = {kmers[i]: int(cnt[i]) for i in idx}
    g = Graph(counts, args.edge_ratio) if counts else None
    popped = g.pop_bubbles() if g else 0
    unitigs = g.unitigs() if g else []
    lap(f"step 3: {len(unitigs):,} unitigs after popping {popped} bubble arms")

    masked = read_mask_kmers(args.mask, k) if args.mask else set()
    known = None
    if args.fasta:
        known = {}
        for name, seq in mk.read_fasta(args.fasta):
            for i in range(len(seq) - k + 1):
                km = seq[i:i + k]
                if set(km) <= set("ACGT"):
                    known.setdefault(mk.canonical(km), set()).add(name)

    rows, long_fa, short_fa = [], [], []
    for n, (seq, _) in enumerate(unitigs):
        ks = [seq[i:i + k] for i in range(len(seq) - k + 1)]
        can = [mk.canonical(x) for x in ks]
        j = [pos_in[c] for c in can if c in pos_in]
        cn_test = float(np.median(kept_cn[j])) if j else np.nan
        per_strain = np.median(kept_M[:, j], axis=1) if j else np.full(len(cohort), np.nan)
        comp = float(np.median(per_strain)) if j else np.nan
        seed_frac = float(kept_seed[j].mean()) if j else 0.0
        if seed_frac < args.min_seed_fraction:
            continue           # reached only by relaxed k-mers: no seed support
        # a comparison strain carries the sequence when it has at least 3 copies and a quarter of the test copy number
        carrier_cut = max(3.0, 0.25 * cn_test)
        row = {"id": f"u{n + 1}", "length": len(seq), "n_kmers": len(ks), "median_cn_test": round(cn_test, 2),
               "median_cn_comparison": round(comp, 2),
               "fold": round(float(np.median(kept_fold[j])), 2) if j else np.nan,
               "seed_fraction": round(seed_frac, 3),
               "n_carriers": int((per_strain >= carrier_cut).sum()),
               "fraction_masked": round(sum(c in masked for c in can) / len(can), 3) if masked else 0.0}
        for s, v in zip(cohort, per_strain):
            row[f"cn_{s}"] = round(float(v), 2)
        if known is not None:
            hits = Counter(t for c in can for t in known.get(c, ()))
            top = hits.most_common(1)
            row["best_library_match"] = top[0][0] if top else ""
            row["fraction_in_library"] = round(top[0][1] / len(ks), 3) if top else 0.0
        row["artifact"] = bool(row["fraction_masked"] >= 0.5)
        record = f">{row['id']} len={len(seq)} cn_test={cn_test:.0f}\n{seq}\n"
        if len(seq) >= args.min_unitig:
            rows.append(row)
            long_fa.append((row["fold"] * len(seq), record))
        else:
            short_fa.append(record)

    tab = pd.DataFrame(rows)
    if len(tab):
        tab["score"] = tab.fold * tab.length
        tab = tab.sort_values("score", ascending=False)
        order = {r.id: i for i, r in enumerate(tab.itertuples())}
        long_fa.sort(key=lambda t: order[t[1].split()[0][1:]])
    lap(f"done: {len(rows)} unitigs reported")
    tab.to_csv(out / "unitigs.tsv", sep="\t", index=False)
    (out / "unitigs.fa").write_text("".join(r for _, r in long_fa))
    (out / "unitigs_short.fa").write_text("".join(short_fa))
    summary = {"positive": args.positive, "comparison_strains": len(cohort), "k": k,
               "single_copy_depth": b_pos, "step1_kmers": len(kmers), "seed_kmers": int(seed.sum()),
               "graph_kmers": int(keep.sum()),
               "bubble_arms_popped": popped, "unitigs_total": len(unitigs), "unitigs_reported": len(rows),
               "unitigs_short_not_reported": len(short_fa),
               "unitigs_flagged_artifact": int(tab.artifact.sum()) if len(tab) else 0,
               "bp_reported": int(tab.length.sum()) if len(tab) else 0}
    pd.Series(summary).to_csv(out / "summary.tsv", sep="\t", header=False)
    print(pd.Series(summary).to_string(), file=sys.stderr)
    return tab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positive", required=True, help="strain suspected to have a burst")
    ap.add_argument("--cohort", nargs="+", required=True, help="comparison strains")
    ap.add_argument("--workdir", default="work")
    ap.add_argument("--k", type=int, default=21, help="odd k")
    ap.add_argument("--threads", type=int, default=1, help="comparison strains looked up at the same time")
    ap.add_argument("--min-cn", type=float, default=20.0, help="minimum copy number in the test strain")
    ap.add_argument("--min-fold", type=float, default=10.0,
                    help="test copy number must be >= this x (comparison median + 1)")
    ap.add_argument("--extend-fold", type=float, default=0.0,
                    help="relaxed fold for graph neighbours of seed k-mers (hysteresis). 0 = same as --min-fold")
    ap.add_argument("--min-seed-fraction", type=float, default=0.25,
                    help="keep a unitig only if at least this fraction of its k-mers are seeds")
    ap.add_argument("--min-depth", type=float, default=15.0,
                    help="smallest single-copy depth accepted for the test strain")
    ap.add_argument("--allow-low-depth", action="store_true", help="run even below --min-depth")
    ap.add_argument("--edge-ratio", type=float, default=0.1,
                    help="drop graph edges to k-mers rarer than this fraction of their siblings")
    ap.add_argument("--min-unitig", type=int, default=60, help="shortest unitig to report, in bp")
    ap.add_argument("--mask", nargs="*", default=[], help="FASTA files (.gz ok) of artifacts to tag, "
                    "such as adapters, PhiX, vectors, organelles, rDNA")
    ap.add_argument("--fasta", help="optional library (for a check only)")
    ap.add_argument("--tmpdir", help="directory for temporary files (default: system)")
    ap.add_argument("--out", required=True)
    run(ap.parse_args())


if __name__ == "__main__":
    main()
