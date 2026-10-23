#!/usr/bin/env bash
# Download the ENA FASTQ files for one accession of a panel run list.
# usage: download_panel.sh RUNS_TSV ACCESSION OUTDIR   (needs aria2c; env CONNECTIONS, default 4)
# RUNS_TSV columns: sample_alias run_accession fastq_ftp fastq_bytes (from select_panel.py)
set -euo pipefail
runs=$1 acc=$2 outdir=$3
dest="$outdir/$acc"
mkdir -p "$dest"
awk -F'\t' -v a="$acc" 'NR>1 && $1==a {print $3"\t"$4}' "$runs" | while IFS=$'\t' read -r urls sizes; do
  IFS=';' read -ra U <<< "$urls"
  IFS=';' read -ra B <<< "$sizes"
  for i in "${!U[@]}"; do
    f=$(basename "${U[$i]}")
    if [[ -f $dest/$f && $(stat -c %s "$dest/$f") == "${B[$i]}" ]]; then continue; fi
    # -x/-s: connections per server / splits; -c resumes; --max-tries covers transient ENA errors
    aria2c --console-log-level=warn --summary-interval=0 -x "${CONNECTIONS:-4}" -s "${CONNECTIONS:-4}" -k 1M -c --max-tries=5 --retry-wait=10 \
      --file-allocation=none -d "$dest" -o "$f" "https://${U[$i]}"
    got=$(stat -c %s "$dest/$f")
    [[ $got == "${B[$i]}" ]] || { echo "size mismatch $f: got $got expected ${B[$i]}" >&2; exit 1; }
  done
done
echo "done $acc: $(ls "$dest" | wc -l) files"
