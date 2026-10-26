#!/usr/bin/env bash
# SLURM step for one accession: convert the prefetched .sra files, check read counts against
# ENA, count k-mers, then remove the FASTQ and .sra files that have been used.
# usage: sra_count_task.sh RUNS_TSV STRAIN SRA_ROOT FASTQ_ROOT OUTDIR KS THREADS HASH
# env: PROJ_ROOT (project directory), SCRATCH (node-local temporary files)
set -euo pipefail
runs=$1 strain=$2 sra=$3 fq=$4 outdir=$5 ks=$6 threads=${7:-8} hash=${8:-2G}
root=${PROJ_ROOT:?set PROJ_ROOT}
scratch=${SCRATCH:?SCRATCH is not set}
lastk=${ks##*,}
if [[ -s $outdir/gcbias/$strain.k$lastk.gcbias.tsv ]]; then echo "$strain already done"; exit 0; fi

source /etc/profile.d/modules.sh
module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7 jellyfish/2.3.0

r1s=() r2s=() sras=()
while IFS=$'\t' read -r s run expected; do
  [[ $s == "$strain" ]] || continue
  dir="$fq/$strain/$run"
  mkdir -p "$dir"
  parallel-fastq-dump --sra-id "$sra/$run/$run.sra" --threads "$threads" --outdir "$dir" \
    --split-files --gzip --tmpdir "$scratch"
  for mate in 1 2; do
    f="$dir/${run}_$mate.fastq.gz"
    [[ -s $f ]] || { echo "sra_count_task: missing $f" >&2; exit 1; }
  done
  # Count read pairs: a record is 4 lines, so the count is lines / 4.
  n1=$(( $(gzip -dc "$dir/${run}_1.fastq.gz" | wc -l) / 4 ))
  n2=$(( $(gzip -dc "$dir/${run}_2.fastq.gz" | wc -l) / 4 ))
  if [[ $n1 != "$expected" || $n2 != "$expected" ]]; then
    echo "sra_count_task: $run mismatch: converted $n1/$n2 pairs, ENA expects $expected" >&2
    exit 1
  fi
  echo "sra_count_task: $run ok, $expected pairs"
  r1s+=("$dir/${run}_1.fastq.gz"); r2s+=("$dir/${run}_2.fastq.gz"); sras+=("$sra/$run/$run.sra")
done < <(tail -n +2 "$runs")

[[ ${#r1s[@]} -gt 0 ]] || { echo "sra_count_task: no runs for $strain" >&2; exit 1; }
join() { local IFS=,; echo "$*"; }
cd "$root"
PROJ_ROOT="$root" scripts/profile_stream.sh "$strain" "$outdir" "$ks" "$(join "${r1s[@]}")" "$(join "${r2s[@]}")" "$threads" "$hash" 20

# All counts are written and verified, so the reads are no longer needed here.
# NCBI keeps the originals. Remove the FASTQ and the .sra files.
rm -rf "$fq/$strain"
rm -f "${sras[@]}"
echo "sra_count_task: $strain done; FASTQ and .sra removed"
