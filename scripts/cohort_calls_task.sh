#!/usr/bin/env bash
# One SLURM job: evaluate mPing recovery in every strain, then call bursts for a list of test strains,
# with and without a limit on how many other strains may carry the sequence.
# usage (sbatch): sbatch scripts/cohort_calls_task.sh "STRAIN1 STRAIN2 ..."
set -euo pipefail
strains=$1
root=${PROJ_ROOT:?set PROJ_ROOT}
source /etc/profile.d/modules.sh
module load jellyfish/2.3.0
cd "$root"
res=/opt/linux/rocky/8.x/x86_64/pkgs/BBMap/39.84/bin/resources
mask=("$res/adapters.fa" "$res/phix174_ill.ref.fa.gz" "$root/resources/masks/NC_001320.1.fa" "$root/resources/masks/NC_011033.1.fa")
common=(--workdir work --k 21 --out work/matrix --mask "${mask[@]}" --fasta RiceTE.fa)
python3 scripts/cohort_matrix.py evaluate "${common[@]}" --element mPing
for s in $strains; do
  python3 scripts/cohort_matrix.py call "${common[@]}" --strain "$s" --tag all
  python3 scripts/cohort_matrix.py call "${common[@]}" --strain "$s" --max-carriers 5 --tag private
done
