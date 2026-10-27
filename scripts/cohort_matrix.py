#!/usr/bin/env python3
"""Cohort-wide k-mer copy-number matrix, and blind burst calls from it.

A k-mer that is high copy in one strain and low in the rest of the cohort marks a burst. This keeps one row
per candidate k-mer and one column per strain, so a strain is judged against the whole cohort, and a repeat shared
by a lineage can be told from a repeat private to one strain.

Stages (each is a subcommand; the heavy ones take --part/--nparts so a SLURM array can split the strains):
  strains   list the strains that have an index and a single-copy depth        -> strains.tsv
  dump      k-mers at or above MIN_CN x single-copy depth, per strain          -> cand/<strain>.txt.zst
  union     every k-mer that is high in at least one strain                    -> candidates.txt.zst
  query     copy number of every candidate in each strain's index              -> cn/<strain>.npy
  build     stack the columns, compute the cohort median                       -> matrix.npy, median.npy
  call      unitigs of the k-mers that stand out in one strain                 -> <out>/calls/<strain>/
  evaluate  how much of a known element each strain's calls recover            -> evaluate.tsv

usage: cohort_matrix.py STAGE --workdir work --k 21 --out work/matrix [stage options]
Needs jellyfish on PATH (dump and query stages) and zstd.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mitekmer as mk  # noqa: E402
import discover_bursts as db  # noqa: E402


# ---------- helpers ----------

def zstd_write_lines(path, lines):
    p = subprocess.Popen(["zstd", "-q", "-T0", "-f", "-o", str(path)], stdin=subprocess.PIPE, text=True)
    for x in lines:
        p.stdin.write(x + "\n")
    p.stdin.close()
    if p.wait() != 0:
        raise RuntimeError(f"zstd failed for {path}")


def zstd_read_lines(path):
    p = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE, text=True)
    out = [line.rstrip("\n") for line in p.stdout]
    if p.wait() != 0:
        raise RuntimeError(f"zstd failed for {path}")
    return out


def my_part(items, part, nparts):
    return [x for i, x in enumerate(items) if i % nparts == part]


def load_strains(out):
    t = pd.read_csv(Path(out) / "strains.tsv", sep="\t")
    return list(t.strain), dict(zip(t.strain, t.baseline))


# ---------- stages ----------

def stage_strains(a):
    work, k = Path(a.workdir), a.k
    rows = []
    for p in sorted(work.glob(f"*.k{k}.jf")):
        s = p.name[:-len(f".k{k}.jf")]
        if not (work / f"{s}.k{k}.histo").exists():
            continue
        b = db.single_copy_depth(work, s, k)
        if np.isfinite(b):
            rows.append((s, float(b)))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["strain", "baseline"]).to_csv(out / "strains.tsv", sep="\t", index=False)
    print(f"{len(rows)} strains with an index and a baseline", file=sys.stderr)


def stage_dump(a):
    strains, base = load_strains(a.out)
    (Path(a.out) / "cand").mkdir(exist_ok=True)
    for s in my_part(strains, a.part, a.nparts):
        f = Path(a.out) / "cand" / f"{s}.txt.zst"
        if f.exists():
            continue
        kmers, _ = db.dump_high(Path(a.workdir) / f"{s}.k{a.k}.jf", a.min_cn * base[s])
        zstd_write_lines(f, kmers)
        print(f"dump {s}: {len(kmers):,} k-mers", file=sys.stderr)


def stage_union(a):
    strains, _ = load_strains(a.out)
    seen = set()
    for s in strains:
        f = Path(a.out) / "cand" / f"{s}.txt.zst"
        if not f.exists():
            sys.exit(f"missing {f}; run the dump stage first")
        seen.update(zstd_read_lines(f))
    zstd_write_lines(Path(a.out) / "candidates.txt.zst", sorted(seen))
    print(f"union: {len(seen):,} candidate k-mers from {len(strains)} strains", file=sys.stderr)


def load_candidates(out):
    return zstd_read_lines(Path(out) / "candidates.txt.zst")


def query_fasta(index, fasta, kmers):
    """Counts of the k-mers in `fasta` (one per record, in the order of `kmers`)."""
    out = subprocess.run(["jellyfish", "query", "-s", str(fasta), str(index)],
                         check=True, capture_output=True, text=True).stdout.split("\n")
    if out and not out[-1]:
        out.pop()
    if len(out) != len(kmers):
        raise RuntimeError(f"query returned {len(out)} lines for {len(kmers)} k-mers ({index})")
    pairs = [x.split() for x in out]
    if [p[0] for p in pairs] != kmers:
        raise RuntimeError(f"query output order differs from the k-mers asked for ({index})")
    return np.array([int(p[1]) for p in pairs], dtype=np.int64)


def stage_query(a):
    strains, base = load_strains(a.out)
    (Path(a.out) / "cn").mkdir(exist_ok=True)
    todo = [s for s in my_part(strains, a.part, a.nparts) if not (Path(a.out) / "cn" / f"{s}.npy").exists()]
    if not todo:
        return
    kmers = load_candidates(a.out)
    with tempfile.TemporaryDirectory(dir=a.tmpdir or None) as tmp:
        fa = Path(tmp) / "candidates.fa"
        with open(fa, "w") as fh:
            for i, km in enumerate(kmers):
                fh.write(f">{i}\n{km}\n")
        for s in todo:
            cnt = query_fasta(Path(a.workdir) / f"{s}.k{a.k}.jf", fa, kmers)
            np.save(Path(a.out) / "cn" / f"{s}.npy", (cnt / base[s]).astype(np.float32))
            print(f"query {s}: {len(kmers):,} k-mers", file=sys.stderr)


def stage_build(a):
    strains, _ = load_strains(a.out)
    cols = []
    for s in strains:
        f = Path(a.out) / "cn" / f"{s}.npy"
        if not f.exists():
            sys.exit(f"missing {f}; run the query stage first")
        cols.append(np.load(f))
    M = np.column_stack(cols).astype(np.float32)          # k-mers x strains
    np.save(Path(a.out) / "matrix.npy", M)
    np.save(Path(a.out) / "median.npy", np.median(M, axis=1).astype(np.float32))
    print(f"matrix {M.shape[0]:,} k-mers x {M.shape[1]} strains", file=sys.stderr)


def load_matrix(out):
    # about 0.75 GB for 2.3 million k-mers x 80 strains; reading a column from disk is slow
    return np.load(Path(out) / "matrix.npy"), np.load(Path(out) / "median.npy")


# ---------- calling ----------

def standout_mask(cn, med, min_cn, min_fold):
    """k-mers whose copy number in a strain is at least min_cn and at least min_fold x (cohort median + 1)."""
    return (cn >= min_cn) & (cn >= min_fold * (med + 1.0))


def carrier_counts(M_rows, cn, strain_col):
    """For each k-mer, the number of OTHER strains that carry it: at least 3 copies and a quarter of this strain's."""
    cut = np.maximum(3.0, 0.25 * cn)[:, None]
    mine = np.zeros(M_rows.shape[1], dtype=bool)
    mine[strain_col] = True
    return ((M_rows >= cut) & ~mine[None, :]).sum(axis=1)


def stage_call(a):
    strains, base = load_strains(a.out)
    if a.strain not in strains:
        sys.exit(f"{a.strain} is not in the cohort")
    j = strains.index(a.strain)
    if base[a.strain] < a.min_depth and not a.allow_low_depth:
        sys.exit(f"{a.strain}: single-copy depth {base[a.strain]:.0f} is below --min-depth {a.min_depth:g}")
    M, med = load_matrix(a.out)
    cn = np.asarray(M[:, j])
    sel = np.flatnonzero(standout_mask(cn, med, a.min_cn, a.min_fold))
    kmers_all = load_candidates(a.out)
    kmers = [kmers_all[i] for i in sel]
    rows_sel = np.asarray(M[sel, :])
    carriers = carrier_counts(rows_sel, cn[sel], j)
    if a.max_carriers is not None:
        keep = carriers <= a.max_carriers
        sel, kmers, rows_sel, carriers = sel[keep], [x for x, k_ in zip(kmers, keep) if k_], rows_sel[keep], carriers[keep]
    print(f"{a.strain}: {len(sel):,} k-mers stand out (cn >= {a.min_cn}, fold >= {a.min_fold}"
          f"{'' if a.max_carriers is None else f', carriers <= {a.max_carriers}'})", file=sys.stderr)

    out = Path(a.out) / "calls" / (a.strain + (f"_{a.tag}" if a.tag else ""))
    out.mkdir(parents=True, exist_ok=True)
    counts = {km: int(round(float(c) * 100)) + 1 for km, c in zip(kmers, cn[sel])}
    pos_in = {km: i for i, km in enumerate(kmers)}
    g = db.Graph(counts, a.edge_ratio) if counts else None
    popped = g.pop_bubbles() if g else 0
    unitigs = g.unitigs() if g else []

    masked = db.read_mask_kmers(a.mask, a.k) if a.mask else set()
    known = None
    if a.fasta:
        known = {}
        for name, seq in mk.read_fasta(a.fasta):
            for i in range(len(seq) - a.k + 1):
                km = seq[i:i + a.k]
                if set(km) <= set("ACGT"):
                    known.setdefault(mk.canonical(km), set()).add(name)

    cn_sel, med_sel = cn[sel], med[sel]
    rows, long_fa, short_fa = [], [], []
    for n, (seq, _) in enumerate(unitigs):
        ks = [seq[i:i + a.k] for i in range(len(seq) - a.k + 1)]
        can = [mk.canonical(x) for x in ks]
        idx = [pos_in[c] for c in can if c in pos_in]
        cn_test = float(np.median(cn_sel[idx])) if idx else np.nan
        cohort_med = float(np.median(med_sel[idx])) if idx else np.nan
        row = {"id": f"u{n + 1}", "length": len(seq), "n_kmers": len(ks),
               "median_cn_test": round(cn_test, 2), "median_cn_cohort": round(cohort_med, 2),
               "fold": round(cn_test / (cohort_med + 1.0), 2),
               "n_carriers": int(np.median(carriers[idx])) if idx else 0,
               "fraction_masked": round(sum(c in masked for c in can) / len(can), 3) if masked else 0.0}
        if known is not None:
            hits = Counter(t for c in can for t in known.get(c, ()))
            top = hits.most_common(1)
            row["best_library_match"] = top[0][0] if top else ""
            row["fraction_in_library"] = round(top[0][1] / len(ks), 3) if top else 0.0
        row["artifact"] = bool(row["fraction_masked"] >= 0.5)
        rec = f">{row['id']} len={len(seq)} cn_test={cn_test:.0f}\n{seq}\n"
        if len(seq) >= a.min_unitig:
            rows.append(row)
            long_fa.append((row["id"], rec))
        else:
            short_fa.append(rec)
    tab = pd.DataFrame(rows)
    if len(tab):
        tab["score"] = tab.fold * tab.length
        tab = tab.sort_values("score", ascending=False)
        order = {r.id: i for i, r in enumerate(tab.itertuples())}
        long_fa.sort(key=lambda t: order[t[0]])
    tab.to_csv(out / "unitigs.tsv", sep="\t", index=False)
    (out / "unitigs.fa").write_text("".join(r for _, r in long_fa))
    (out / "unitigs_short.fa").write_text("".join(short_fa))
    summary = {"strain": a.strain, "k": a.k, "single_copy_depth": base[a.strain], "cohort_strains": len(strains),
               "standout_kmers": len(sel), "bubble_arms_popped": popped, "unitigs_total": len(unitigs),
               "unitigs_reported": len(rows), "unitigs_short_not_reported": len(short_fa),
               "unitigs_flagged_artifact": int(tab.artifact.sum()) if len(tab) else 0}
    pd.Series(summary).to_csv(out / "summary.tsv", sep="\t", header=False)
    print(pd.Series(summary).to_string(), file=sys.stderr)
    return tab


def stage_evaluate(a):
    """Recovery of a known element's unique k-mers in each strain, using the same call rule as `call`.
    The library is used only to score the calls."""
    strains, base = load_strains(a.out)
    M, med = load_matrix(a.out)
    kt = mk.te_kmer_table(a.fasta, a.k)
    target = kt[(kt.te == a.element) & (kt.multiplicity == 1)].kmer.tolist()
    kmers_all = load_candidates(a.out)
    index = {km: i for i, km in enumerate(kmers_all)}
    rows_t = [index[x] for x in target if x in index]
    print(f"{a.element}: {len(target)} unique k-mers, {len(rows_t)} are candidates", file=sys.stderr)
    sub = np.asarray(M[rows_t, :])
    rows = []
    for j, s in enumerate(strains):
        cn = np.asarray(M[:, j])
        called = standout_mask(sub[:, j], med[rows_t], a.min_cn, a.min_fold)
        rows.append({"strain": s, "baseline": base[s],
                     f"median_cn_{a.element}": float(np.median(sub[:, j])),
                     "kmers_called": int(called.sum()), "fraction_called": float(called.mean()),
                     "genome_kmers_standing_out": int(standout_mask(cn, med, a.min_cn, a.min_fold).sum())})
    tab = pd.DataFrame(rows)
    tab.to_csv(Path(a.out) / f"evaluate_{a.element}.tsv", sep="\t", index=False)
    print(tab.sort_values(f"median_cn_{a.element}", ascending=False).head(12).round(2).to_string(index=False),
          file=sys.stderr)
    return tab


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["strains", "dump", "union", "query", "build", "call", "evaluate"])
    ap.add_argument("--workdir", default="work")
    ap.add_argument("--k", type=int, default=21, help="odd k")
    ap.add_argument("--out", default="work/matrix")
    ap.add_argument("--part", type=int, default=0, help="this task's part (0-based) for dump and query")
    ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--tmpdir")
    ap.add_argument("--min-cn", type=float, default=20.0, help="minimum copy number in the strain")
    ap.add_argument("--min-fold", type=float, default=5.0, help="strain must be >= this x (cohort median + 1)")
    ap.add_argument("--max-carriers", type=int, default=None,
                    help="call: keep k-mers carried by at most this many other strains")
    ap.add_argument("--strain", help="call: the test strain")
    ap.add_argument("--tag", default="", help="call: suffix for the output directory, to keep settings apart")
    ap.add_argument("--min-depth", type=float, default=15.0)
    ap.add_argument("--allow-low-depth", action="store_true")
    ap.add_argument("--edge-ratio", type=float, default=0.1)
    ap.add_argument("--min-unitig", type=int, default=60)
    ap.add_argument("--mask", nargs="*", default=[])
    ap.add_argument("--fasta", help="library, used only to label or score calls")
    ap.add_argument("--element", default="mPing", help="evaluate: element of --fasta to score")
    a = ap.parse_args()
    if a.k % 2 == 0:
        sys.exit("k must be odd")
    for tool in ("zstd",):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not on PATH")
    {"strains": stage_strains, "dump": stage_dump, "union": stage_union, "query": stage_query,
     "build": stage_build, "call": stage_call, "evaluate": stage_evaluate}[a.stage](a)


if __name__ == "__main__":
    main()
