#!/usr/bin/env bash
# SLURM array task: run discover_bursts.py for one test strain.
# usage (sbatch): sbatch --array=1-N -c 8 scripts/discover_task.sh JOBS_FILE OUTROOT [extra discover_bursts.py args]
# JOBS_FILE: one job per line, "strain<TAB>min_fold<TAB>extend_fold" (extend_fold 0 = no hysteresis).
# Output goes to OUTROOT/f<min_fold>e<extend_fold>/<strain>.
# env: COHORT_FILE  file with the comparison strains, space separated (default: the 10 strains below).
# The test strain is dropped from the comparison set if present.
set -euo pipefail
jobs=$1 outroot=$2
shift 2
root=${PROJ_ROOT:?set PROJ_ROOT}
read -r strain fold ext < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$jobs")
source /etc/profile.d/modules.sh
module load jellyfish/2.3.0
cd "$root"
res=/opt/linux/rocky/8.x/x86_64/pkgs/BBMap/39.84/bin/resources
if [[ -n ${COHORT_FILE:-} ]]; then cohort=$(cat "$COHORT_FILE")
else cohort="B031 CX269 IRIS_313-8171 IRIS_313-11652 IRIS_313-8321 IRIS_313-10871 IRIS_313-11418 IRIS_313-8385 IRIS_313-10761 NB"; fi
python3 scripts/discover_bursts.py --positive "$strain" --cohort $cohort --workdir work --k 21 \
  --min-fold "$fold" --extend-fold "$ext" --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --mask "$res/adapters.fa" "$res/phix174_ill.ref.fa.gz" "$root/resources/masks/NC_001320.1.fa" "$root/resources/masks/NC_011033.1.fa" \
  --fasta RiceTE.fa --tmpdir "${SCRATCH:?}" --out "$outroot/f${fold}e${ext}/$strain" "$@"
