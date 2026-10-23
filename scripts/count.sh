#!/usr/bin/env bash
# Build one jellyfish index per strain (R1+R2 together), histogram, and TE k-mer depths.
# usage: count.sh STRAIN R1 R2 K [THREADS] [HASH_SIZE] [MIN_QUAL]
# R1 may be a single file, or a directory of read pairs (R2 is then ignored, pass "-").
# In a directory, R1 files match *pair1*, *_R1* or *_1.f*; the mate is found by
# swapping pair1->pair2, _R1->_R2 or _1.->_2. All pairs are counted into one index.
set -euo pipefail
strain=$1 r1=$2 r2=$3 k=$4
threads=${5:-8} hash=${6:-4G} minq=${7:-20}
here=$(cd "$(dirname "$0")/.." && pwd)
work=${WORKDIR:-$here/work}
fasta=${TE_FASTA:-$here/RiceTE.fa}
mkdir -p "$work"

# Fill R1S/R2S with matched file lists.
R1S=() R2S=()
if [[ -d $r1 ]]; then
  while IFS= read -r f; do
    base=$(basename "$f"); dir=$(dirname "$f")
    case $base in
      *pair1*) mate=${base/pair1/pair2} ;;
      *_R1*)   mate=${base/_R1/_R2} ;;
      *_1.f*)  mate=${base/_1./_2.} ;;
      *) continue ;;
    esac
    [[ -s $dir/$mate ]] || { echo "count.sh: no mate $dir/$mate for $f" >&2; exit 1; }
    R1S+=("$f"); R2S+=("$dir/$mate")
  done < <(find -L "$r1" -maxdepth 1 -type f \( -name '*.fastq*' -o -name '*.fq*' \) | sort)
  (( ${#R1S[@]} )) || { echo "count.sh: no read pairs found in $r1" >&2; exit 1; }
  echo "count.sh: $strain: ${#R1S[@]} read pairs from $r1" >&2
else
  R1S=("$r1") R2S=("$r2")
fi
# Read-file decompression: .gz or plain
cat_reads() { for f in "$@"; do if [[ $f == *.gz ]]; then gzip -dc "$f"; else cat "$f"; fi; done; }

jf="$work/$strain.k$k.jf"

if [[ ! -s $jf ]]; then
  # In a SLURM job write to node-local $SCRATCH, then move to shared storage
  tmp=${SCRATCH:+$SCRATCH/}$strain.k$k.jf.tmp
  [[ -n ${SCRATCH:-} ]] || tmp=$jf.tmp
  # --quality-start 33: jellyfish defaults to 64, these libraries are Phred+33
  jellyfish count -m "$k" -s "$hash" -t "$threads" -C \
    --quality-start 33 --min-quality "$minq" \
    -F 2 -o "$tmp" <(cat_reads "${R1S[@]}") <(cat_reads "${R2S[@]}")
  mv "$tmp" "$jf"
fi
jellyfish histo -t "$threads" -h 1000 "$jf" > "$work/$strain.k$k.histo"
python3 "$here/scripts/mitekmer.py" query --index "$jf" --fasta "$fasta" --k "$k" \
  --out "$work/$strain.k$k.depth.tsv"
