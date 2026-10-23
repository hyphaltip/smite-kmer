#!/usr/bin/env bash
# Submit a SLURM array that streams ENA reads into profile_stream.sh, one accession per task.
# usage: submit_stream.sh LIST_TSV [KS] [MAX_CONCURRENT]
# LIST_TSV columns: strain ... r1 r2 (r1/r2 are the last two columns; comma-separated URLs).
# env: PARTITION ACCOUNT (default exfab) CPUS MEM TIME KEEP_INDEX (default 1) OUTDIR (default work)
# Tasks whose final gcbias table exists are skipped, so a rerun does only the missing lines.
set -euo pipefail
list=$1 ks=${2:-21,31} maxc=${3:-4}
here=$(cd "$(dirname "$0")/.." && pwd)
part=${PARTITION:-exfab} acct=${ACCOUNT:-exfab}
cpus=${CPUS:-8} mem=${MEM:-48G} time=${TIME:-08:00:00}
outdir=${OUTDIR:-$here/work} keep=${KEEP_INDEX-1}
n=$(( $(wc -l < "$list") - 1 ))
lastk=${ks##*,}
mkdir -p "$here/logs" "$outdir"
sbatch --parsable -p "$part" -A "$acct" -J "stream.k${ks//,/_}" -c "$cpus" --mem="$mem" -t "$time" \
  --array="1-$n%$maxc" -o "$here/logs/stream.%A_%a.out" \
  --export=ALL,PROJ_ROOT="$here",KEEP_INDEX="$keep" --wrap \
  "module load jellyfish/2.3.0; cd '$here'; \
   IFS=\$'\\t' read -r s r1 r2 < <(awk -F'\\t' -v n=\$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print \$1\"\\t\"\$(NF-1)\"\\t\"\$NF}' '$list'); \
   if [[ -s '$outdir/gcbias/'\$s.k$lastk.gcbias.tsv ]]; then echo \"\$s already done\"; exit 0; fi; \
   scripts/profile_stream.sh \"\$s\" '$outdir' '$ks' \"\$r1\" \"\$r2\" \$((SLURM_CPUS_PER_TASK/2)) 2G 20"
