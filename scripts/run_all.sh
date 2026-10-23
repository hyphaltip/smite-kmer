#!/usr/bin/env bash
# Count every strain in the sample sheet for one k, then analyze.
# usage: run_all.sh K [THREADS] [HASH_SIZE] [REFERENCE]
set -euo pipefail
k=$1 threads=${2:-8} hash=${3:-2G} ref=${4:-NB}
here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here"
mkdir -p logs work results
tail -n +2 samples.tsv | while IFS=$'\t' read -r s r1 r2 _; do
  /usr/bin/time -v scripts/count.sh "$s" "$r1" "${r2:--}" "$k" "$threads" "$hash" 20 2> "logs/$s.k$k.time"
done
python3 scripts/mitekmer.py analyze --fasta RiceTE.fa --samples samples.tsv \
  --workdir work --k "$k" --reference "$ref" --out "results/full"
