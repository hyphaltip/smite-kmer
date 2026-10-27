#!/usr/bin/env bash
# SLURM array wrapper: runs sra_node_task.sh for accession number SLURM_ARRAY_TASK_ID in STRAINS.
# usage (sbatch): sbatch --array=1-N%4 scripts/sra_node_array_task.sh RUNS_TSV STRAINS OUTDIR KS
set -euo pipefail
runs=$1 strains=$2 outdir=$3 ks=$4
root=${PROJ_ROOT:?set PROJ_ROOT}
strain=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$strains" | cut -f1)  # no header line in the strains file
[[ -n $strain ]] || { echo "no accession for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
echo "task $SLURM_ARRAY_TASK_ID: $strain"
"$root/scripts/sra_node_task.sh" "$runs" "$strain" "$outdir" "$ks" "${SLURM_CPUS_PER_TASK:-8}" 2G
