#!/usr/bin/env bash
# One SLURM array task: stream the accession on line SLURM_ARRAY_TASK_ID+1 of LIST.
# usage (from submit_stream.sh): stream_task.sh LIST_TSV KS OUTDIR
# env: PROJ_ROOT (project directory), SLURM_ARRAY_TASK_ID, SLURM_CPUS_PER_TASK
# LIST_TSV: strain in column 1, R1 and R2 URL lists in the last two columns.
set -euo pipefail
list=$1 ks=$2 outdir=$3
root=${PROJ_ROOT:?set PROJ_ROOT to the project directory}
if ! command -v module > /dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
fi
if command -v module > /dev/null 2>&1; then module load jellyfish/2.3.0; fi
cd "$root"
line=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n{print $1"\t"$(NF-1)"\t"$NF}' "$list")
IFS=$'\t' read -r s r1 r2 <<< "$line"
[[ -n $s && -n $r1 && -n $r2 ]] || { echo "stream_task: bad list line $((SLURM_ARRAY_TASK_ID + 1)) in $list" >&2; exit 1; }
lastk=${ks##*,}
if [[ -s $outdir/gcbias/$s.k$lastk.gcbias.tsv ]]; then echo "$s already done"; exit 0; fi
scripts/profile_stream.sh "$s" "$outdir" "$ks" "$r1" "$r2" $(( ${SLURM_CPUS_PER_TASK:-2} / 2 )) 2G 20
