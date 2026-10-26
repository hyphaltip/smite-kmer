#!/usr/bin/env bash
# Login-node step: prefetch the .sra file of every run of one accession from NCBI SRA.
# usage: sra_fetch.sh RUNS_TSV STRAIN SRA_ROOT
# RUNS_TSV columns: strain run ena_read_count (from make_sra_runs.py)
# Output: SRA_ROOT/<run>/<run>.sra. Prefetch validates each file and skips complete ones.
set -euo pipefail
runs=$1 strain=$2 root=$3
source /etc/profile.d/modules.sh
module load sratoolkit/3.2.0
mkdir -p "$root"
awk -F'\t' -v s="$strain" 'NR>1 && $1==s {print $2}' "$runs" | while read -r run; do
  [[ -n $run ]] || continue
  prefetch "$run" -O "$root" --max-size 30G
done
echo "sra_fetch: $strain prefetched into $root"
