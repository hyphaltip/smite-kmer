#!/usr/bin/env python3
"""Write the run list for NCBI SRA download: one row per run of each accession.

Runs come from the ENA run table for the Chen 2019 project (PRJEB6180). The read_count there is
the expected number of read pairs, used to check the converted FASTQ.

usage: make_sra_runs.py --strains resources/Chen2019/stream100.tsv --out resources/Chen2019/sra_runs.tsv
Output columns: strain  run  ena_read_count
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strains", required=True, help="TSV with a 'strain' column")
    ap.add_argument("--ena-runs", default=str(ROOT / "resources/Chen2019/PRJEB6180.ena_runs.tsv.gz"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    want = pd.read_csv(a.strains, sep="\t")["strain"].tolist()
    ena = pd.read_csv(a.ena_runs, sep="\t")
    ena = ena[ena.sample_alias.isin(want)][["sample_alias", "run_accession", "read_count"]]
    missing = sorted(set(want) - set(ena.sample_alias))
    if missing:
        sys.exit(f"no ENA runs for: {missing}")
    ena = ena.rename(columns={"sample_alias": "strain", "run_accession": "run",
                              "read_count": "ena_read_count"})
    ena.sort_values(["strain", "run"]).to_csv(a.out, sep="\t", index=False)
    print(f"{ena.strain.nunique()} accessions, {len(ena)} runs -> {a.out}")


if __name__ == "__main__":
    main()
