#!/usr/bin/env bash
# SLURM task for one cohort_matrix.py stage.
# usage (sbatch): sbatch [--array=0-N] scripts/cohort_task.sh STAGE [extra cohort_matrix.py args]
# For array tasks, SLURM_ARRAY_TASK_ID is the part and NPARTS (env) is the number of parts.
set -euo pipefail
stage=$1
shift
root=${PROJ_ROOT:?set PROJ_ROOT}
source /etc/profile.d/modules.sh
module load jellyfish/2.3.0
cd "$root"
extra=()
if [[ -n ${SLURM_ARRAY_TASK_ID:-} ]]; then extra=(--part "$SLURM_ARRAY_TASK_ID" --nparts "${NPARTS:?set NPARTS}"); fi
python3 scripts/cohort_matrix.py "$stage" --workdir work --k 21 --out "${MATRIX_DIR:-$root/work/matrix}" \
  --tmpdir "${SCRATCH:?}" "${extra[@]}" "$@"
