#!/usr/bin/env bash
# SLURM array task: rebuild one GC-bias table from an existing jellyfish index.
# usage (sbatch): sbatch --array=1-N scripts/gcbias_rebuild_task.sh LIST_TSV WORKDIR
# LIST_TSV columns: strain k (header line first)
set -euo pipefail
list=$1 work=$2
root=${PROJ_ROOT:?set PROJ_ROOT}
read -r s k < <(tail -n +2 "$list" | sed -n "${SLURM_ARRAY_TASK_ID}p")
source /etc/profile.d/modules.sh
module load jellyfish/2.3.0
cd "$root"
python3 scripts/mitekmer.py gcbias --index "$work/$s.k$k.jf" --k "$k" \
  --histo "$work/$s.k$k.histo" --out "$work/gcbias/$s.k$k.gcbias.tsv"
echo "gcbias_rebuild: $s k=$k done"
