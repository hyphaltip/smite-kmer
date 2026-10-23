#!/usr/bin/env bash
# Build one jellyfish index per strain (R1+R2 together), histogram, and TE k-mer depths.
# usage: count.sh STRAIN R1 R2 K [THREADS] [HASH_SIZE] [MIN_QUAL]
set -euo pipefail
strain=$1 r1=$2 r2=$3 k=$4
threads=${5:-8} hash=${6:-4G} minq=${7:-20}
here=$(cd "$(dirname "$0")/.." && pwd)
work=${WORKDIR:-$here/work}
fasta=${TE_FASTA:-$here/RiceTE.fa}
mkdir -p "$work"
jf="$work/$strain.k$k.jf"

if [[ ! -s $jf ]]; then
  # --quality-start 33: jellyfish defaults to 64, these libraries are Phred+33
  jellyfish count -m "$k" -s "$hash" -t "$threads" -C \
    --quality-start 33 --min-quality "$minq" \
    -F 2 -o "$jf.tmp" <(gzip -dc "$r1") <(gzip -dc "$r2")
  mv "$jf.tmp" "$jf"
fi
jellyfish histo -t "$threads" -h 1000 "$jf" > "$work/$strain.k$k.histo"
python3 "$here/scripts/mitekmer.py" query --index "$jf" --fasta "$fasta" --k "$k" \
  --out "$work/$strain.k$k.depth.tsv"
