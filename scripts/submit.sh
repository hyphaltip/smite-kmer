#!/usr/bin/env bash
# Submit one SLURM count job per strain in samples.tsv, then one analysis job
# that starts after all count jobs succeed.
# usage: submit.sh K [REFERENCE]
# env overrides: PARTITION ACCOUNT THREADS MEM TIME HASH_SIZE ANALYZE_TIME
#   EXTRA_DEP (job ids, colon-separated; count jobs wait on them)  ANALYZE_DEP (extra ids the analysis waits on)
#   ANALYZE_SAMPLES (sheet for the analysis job, default = SAMPLES)
#   SAMPLES (sheet, default samples.tsv)  NO_ANALYZE=1 (submit count jobs only)  OUT_NAME (results subdir, default full)
set -euo pipefail
k=$1 ref=${2:-NB}
here=$(cd "$(dirname "$0")/.." && pwd)
part=${PARTITION:-batch} acct=${ACCOUNT:-stajichlab}
threads=${THREADS:-16} mem=${MEM:-64G} time=${TIME:-24:00:00}
hash=${HASH_SIZE:-2G} atime=${ANALYZE_TIME:-04:00:00}
samples=${SAMPLES:-samples.tsv} out=${OUT_NAME:-full}
mkdir -p "$here/logs" "$here/work" "$here/results"

common=(--parsable -p "$part" -A "$acct" --export=ALL,PROJ_ROOT="$here")
ids=()
while IFS=$'\t' read -r s r1 r2 _; do
  [[ -n $s ]] || continue
  id=$(sbatch "${common[@]}" ${EXTRA_DEP:+--dependency=afterok:$EXTRA_DEP} -J "kmer.$s.k$k" -c "$threads" --mem="$mem" -t "$time" \
    -o "$here/logs/$s.k$k.%j.out" --wrap \
    "module load jellyfish/2.3.0; cd '$here'; scripts/count.sh '$s' '$r1' '${r2:--}' $k $threads $hash 20")
  echo "count $s: job $id"
  ids+=("$id")
done < <(tail -n +2 "$here/$samples")

if [[ -n ${NO_ANALYZE:-} ]]; then echo "count jobs only (NO_ANALYZE set)"; exit 0; fi
dep=$(IFS=:; echo "${ids[*]}")${ANALYZE_DEP:+:$ANALYZE_DEP}
aid=$(sbatch "${common[@]}" -J "kmer.analyze.k$k" -c 4 --mem=16G -t "$atime" \
  --dependency="afterok:$dep" -o "$here/logs/analyze.k$k.%j.out" --wrap \
  "cd '$here'; python3 scripts/mitekmer.py analyze --fasta RiceTE.fa --samples '${ANALYZE_SAMPLES:-$samples}' --workdir work --k $k --reference '$ref' --out results/$out")
echo "analyze: job $aid (after $dep)"
