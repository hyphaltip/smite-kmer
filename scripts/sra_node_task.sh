#!/usr/bin/env bash
# One accession, all on one compute node: prefetch the runs from NCBI SRA into $SCRATCH,
# convert them, check each pair count against ENA, and count k-mers. Only the count tables
# are written to shared storage. $SCRATCH is removed by the scheduler when the job ends.
# usage: sra_node_task.sh RUNS_TSV STRAIN OUTDIR KS [THREADS] [HASH_SIZE]
# env: PROJ_ROOT (project directory), SCRATCH (set by the job)
set -euo pipefail
runs=$1 strain=$2 outdir=$3 ks=$4 threads=${5:-8} hash=${6:-2G}
root=${PROJ_ROOT:?set PROJ_ROOT to the project directory}
scratch=${SCRATCH:?SCRATCH is not set (run inside a SLURM job)}
lastk=${ks##*,}
if [[ -s $outdir/gcbias/$strain.k$lastk.gcbias.tsv ]]; then echo "$strain already done"; exit 0; fi

source /etc/profile.d/modules.sh
module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7 jellyfish/2.3.0
work=$scratch/$strain
sra=$work/sra fq=$work/fq
mkdir -p "$sra" "$fq"

r1s=() r2s=()
while IFS=$'\t' read -r s run expected; do
  [[ $s == "$strain" ]] || continue
  echo "sra_node_task: $(date +%T) prefetch $run"
  prefetch "$run" -O "$sra" --max-size 30G
  dir="$fq/$run"; mkdir -p "$dir"
  echo "sra_node_task: $(date +%T) convert $run"
  parallel-fastq-dump --sra-id "$sra/$run/$run.sra" --threads "$threads" --outdir "$dir" \
    --split-files --gzip --tmpdir "$scratch"
  n1=$(( $(gzip -dc "$dir/${run}_1.fastq.gz" | wc -l) / 4 ))
  n2=$(( $(gzip -dc "$dir/${run}_2.fastq.gz" | wc -l) / 4 ))
  if [[ $n1 != "$expected" || $n2 != "$expected" ]]; then
    echo "sra_node_task: $run mismatch: converted $n1/$n2 pairs, ENA expects $expected" >&2
    exit 1
  fi
  echo "sra_node_task: $run ok, $expected pairs"
  r1s+=("$dir/${run}_1.fastq.gz"); r2s+=("$dir/${run}_2.fastq.gz")
done < <(tail -n +2 "$runs")
[[ ${#r1s[@]} -gt 0 ]] || { echo "sra_node_task: no runs for $strain" >&2; exit 1; }

join() { local IFS=,; echo "$*"; }
cd "$root"
echo "sra_node_task: $(date +%T) count $strain"
PROJ_ROOT="$root" scripts/profile_stream.sh "$strain" "$outdir" "$ks" "$(join "${r1s[@]}")" "$(join "${r2s[@]}")" "$threads" "$hash" 20
echo "sra_node_task: $(date +%T) $strain done"
