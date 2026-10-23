#!/usr/bin/env bash
# Smoke-test subsample: first N reads of each file in samples.tsv -> reads_sub/
# usage: subsample.sh [N_READS=1000000]
set -euo pipefail
n=${1:-1000000}
here=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$here/reads_sub"
tail -n +2 "$here/samples.tsv" | while IFS=$'\t' read -r s r1 r2; do
  for r in "$r1" "$r2"; do
    out="$here/reads_sub/$(basename "$r")"
    # head closes the pipe early; gzip exits with SIGPIPE, which is expected
    { gzip -dc "$here/$r" || true; } | head -n $((n*4)) | gzip > "$out"
  done
done
