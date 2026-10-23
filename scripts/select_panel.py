#!/usr/bin/env python3
"""Pick a low-mPing accession panel from Chen et al. 2019 Supplementary Data 1
joined to the ENA PRJEB6180 run list. Writes a panel table and a run-level download list.

usage: select_panel.py --xlsx X --ena-runs T.tsv[.gz] --out-prefix P [--n-per-group ...]
"""
import argparse
import numpy as np
import pandas as pd

GENOME_MB = 373  # approximate Nipponbare genome size, used only for rough depth

# accessions per subgroup; temperate japonica is low because that group carries the bursters
QUOTA = {"indica": 3, "tropical japonica": 2, "aus/boro": 2, "aromatic": 1,
         "temperate japonica": 2}


def load_table(xlsx):
    import openpyxl
    rows = list(openpyxl.load_workbook(xlsx, read_only=True).worksheets[0].iter_rows(values_only=True))
    t = pd.DataFrame(rows[2:], columns=rows[1]).dropna(subset=["Taxa"])
    t["Taxa"] = t["Taxa"].astype(str).str.strip()
    for c in ("mPing", "Ping", "Pong"):  # one accession has 'NA'
        t[c] = pd.to_numeric(t[c], errors="coerce")
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--ena-runs", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--max-mping", type=int, default=12)
    ap.add_argument("--min-depth", type=float, default=10.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    t = load_table(a.xlsx)
    e = pd.read_csv(a.ena_runs, sep="\t")
    cov = e.groupby("sample_alias")["base_count"].sum().div(GENOME_MB * 1e6).rename("approx_depth")
    nruns = e.groupby("sample_alias").size().rename("n_runs")
    t = t.join(cov, on="Taxa").join(nruns, on="Taxa")
    ok = t[(t.mPing <= a.max_mping) & (t.Ping <= 3) & (t.approx_depth >= a.min_depth)]
    rng = np.random.default_rng(a.seed)
    pick = []
    for grp, n in QUOTA.items():
        pool = ok[ok.Subgroups == grp]
        pick.append(pool.iloc[rng.permutation(len(pool))[:n]])
    panel = pd.concat(pick)
    panel.to_csv(a.out_prefix + ".panel.tsv", sep="\t", index=False)
    runs = e[e.sample_alias.isin(panel.Taxa)][["sample_alias", "run_accession", "fastq_ftp", "fastq_bytes"]]
    runs.to_csv(a.out_prefix + ".runs.tsv", sep="\t", index=False)
    print(panel.drop(columns=["Origin"]).to_string(index=False))
    print(f"\n{len(panel)} accessions, {len(runs)} runs, {runs.fastq_bytes.map(lambda s: sum(int(x) for x in str(s).split(';'))).sum()/1e9:.1f} GB")


if __name__ == "__main__":
    main()
