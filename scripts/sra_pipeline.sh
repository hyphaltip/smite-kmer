#!/usr/bin/env bash
# Driver for the NCBI SRA route. Run on the login node (prefetch needs the network).
# For each accession in STRAINS_TSV: prefetch its runs, then submit one SLURM count job.
# Accessions that already have a final GC-bias table are skipped, so a rerun resumes.
# usage: sra_pipeline.sh STRAINS_TSV RUNS_TSV [KS]
# env: SRA_ROOT (default /bigdata/.../sra_cache), FASTQ_ROOT (default /bigdata/.../fastq_tmp),
#      OUTDIR (default work), ACCOUNT (default stajichlab), PARTITION (default batch),
#      CPUS (default 8), MEM (default 32G), TIME (default 06:00:00), HASH_SIZE (default 2G)
set -euo pipefail
strains=$1 runs=$2 ks=${3:-21,31}
here=$(cd "$(dirname "$0")/.." && pwd)
sra_root=${SRA_ROOT:-$here/sra_cache}
fq_root=${FASTQ_ROOT:-$here/fastq_tmp}
outdir=${OUTDIR:-$here/work}
acct=${ACCOUNT:-stajichlab} part=${PARTITION:-stajichlab}
cpus=${CPUS:-8} mem=${MEM:-32G} time=${TIME:-06:00:00} hash=${HASH_SIZE:-2G}
lastk=${ks##*,}
mkdir -p "$here/logs" "$sra_root" "$fq_root" "$outdir"
for strain in $(tail -n +2 "$strains" | cut -f1); do
  if [[ -s $outdir/gcbias/$strain.k$lastk.gcbias.tsv ]]; then echo "skip $strain (done)"; continue; fi
  "$here/scripts/sra_fetch.sh" "$runs" "$strain" "$sra_root"
  jid=$(sbatch --parsable -A "$acct" -p "$part" -J "sra.$strain" -c "$cpus" --mem="$mem" -t "$time" \
    -o "$here/logs/sra.$strain.%j.out" \
    --export=ALL,PROJ_ROOT="$here" \
    "$here/scripts/sra_count_task.sh" "$runs" "$strain" "$sra_root" "$fq_root" "$outdir" "$ks" "$cpus" "$hash")
  echo "submitted $strain: job $jid"
done
