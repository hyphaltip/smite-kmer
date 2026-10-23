#!/usr/bin/env bash
# Profile one accession from streamed reads and keep only small tables.
# No FASTQ is written to disk. Each k gets its own jellyfish index in $SCRATCH,
# which is deleted after the tables are made.
# usage: profile_stream.sh STRAIN OUTDIR KS R1_LIST R2_LIST [THREADS] [HASH_SIZE] [MIN_QUAL]
#   KS       comma-separated k values, e.g. 21,31
#   R*_LIST  comma-separated URLs (http/https/ftp, or ENA paths like ftp.sra.ebi.ac.uk/...)
#            or local files (.gz or plain)
# Outputs in OUTDIR: STRAIN.k<k>.histo  .depth.tsv  .stats, and gcbias/STRAIN.k<k>.gcbias.tsv
# These names match count.sh, so `mitekmer.py analyze --workdir OUTDIR` reads them.
# KEEP_INDEX=1 moves each jellyfish index to OUTDIR/STRAIN.k<k>.jf instead of deleting it.
# The index is not much smaller than the FASTQ for 10-20x data; see the size notes in the README.
# Both read files go into one stream; k-mer counts do not depend on pairing.
# A failed download fails the job. Do not retry inside the stream: a retry would repeat reads.
set -euo pipefail
strain=$1 outdir=$2 ks=$3 r1s=$4 r2s=$5
threads=${6:-8} hash=${7:-2G} minq=${8:-20}
here=${PROJ_ROOT:?set PROJ_ROOT to the project directory}
fasta=${TE_FASTA:-$here/RiceTE.fa}
scratch=${SCRATCH:?SCRATCH is not set}
mkdir -p "$outdir/gcbias"
tmp=$(mktemp -d "$scratch/$strain.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

stream() {
  local u list
  IFS=, read -ra list <<< "$1"
  for u in "${list[@]}"; do
    case $u in
      http://*|https://*|ftp://*) curl -sS --fail -L "$u" | gzip -dc ;;
      ftp.*)                      curl -sS --fail -L "https://$u" | gzip -dc ;;
      *.gz)                       gzip -dc "$u" ;;
      *)                          cat "$u" ;;
    esac
  done
}

IFS=, read -ra KS <<< "$ks"
fifos=() pids=()
for k in "${KS[@]}"; do
  mkfifo "$tmp/in.k$k"
  fifos+=("$tmp/in.k$k")
  jellyfish count -m "$k" -s "$hash" -t "$threads" -C --quality-start 33 --min-quality "$minq" \
    -o "$tmp/$strain.k$k.jf" "$tmp/in.k$k" &
  pids+=($!)
done
{ stream "$r1s"; stream "$r2s"; } | tee "${fifos[@]}" > /dev/null
for p in "${pids[@]}"; do wait "$p"; done

for k in "${KS[@]}"; do
  jf="$tmp/$strain.k$k.jf"
  jellyfish histo -t "$threads" -h 1000 "$jf" > "$outdir/$strain.k$k.histo"
  jellyfish stats "$jf" > "$outdir/$strain.k$k.stats"
  python3 "$here/scripts/mitekmer.py" query --index "$jf" --fasta "$fasta" --k "$k" \
    --out "$outdir/$strain.k$k.depth.tsv"
  python3 "$here/scripts/mitekmer.py" gcbias --index "$jf" --k "$k" \
    --histo "$outdir/$strain.k$k.histo" --out "$outdir/gcbias/$strain.k$k.gcbias.tsv"
  if [[ -n ${KEEP_INDEX:-} ]]; then mv "$jf" "$outdir/$strain.k$k.jf.tmp" && mv "$outdir/$strain.k$k.jf.tmp" "$outdir/$strain.k$k.jf"; else rm -f "$jf"; fi
done
echo "profile_stream: $strain done for k=$ks (index kept: ${KEEP_INDEX:-no})"
