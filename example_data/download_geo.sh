#!/usr/bin/env bash
# Resumable, checksummed GEO supplementary downloads.
set -uo pipefail
OUT=${1:-data/raw}; shift || true
mkdir -p "$OUT"
fetch() { # url dest
  local url=$1 dest=$2
  for i in 1 2 3 4 5; do
    curl -sS -L --fail -C - --retry 5 --retry-all-errors --max-time 7200 -o "$dest" "$url" && break
    echo "retry $i $url" >&2; sleep 60
  done
  echo -e "$(shasum -a 256 "$dest" | cut -d' ' -f1)\t$(wc -c < "$dest")\t$dest\t$url" >> "$OUT/MANIFEST.tsv"
}
for spec in "$@"; do  # GSE:file
  g=${spec%%:*}; f=${spec#*:}; s=${g:0:${#g}-3}nnn
  if [ "$f" = "matrix" ]; then
    fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/$s/$g/matrix/${g}_series_matrix.txt.gz" "$OUT/${g}_series_matrix.txt.gz"
  else
    fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/$s/$g/suppl/$f" "$OUT/$f"
  fi
done
