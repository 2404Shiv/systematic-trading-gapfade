#!/usr/bin/env bash
set -euo pipefail
mkdir -p data

while read T; do
  [ -z "$T" ] && continue
  # Stooq symbol: dots -> dashes; lowercase; append .us
  url_sym=$(echo "${T//./-}.us" | tr '[:upper:]' '[:lower:]')
  # Save file: dots -> underscores
  out="data/${T//./_}.csv"

  echo "Downloading $T -> $out"
  curl -sL "https://stooq.com/q/d/l/?s=${url_sym}&i=d" -o "$out"

  # Require at least header + 500 rows; otherwise drop it
  rows=$(wc -l < "$out" || echo 0)
  if [ "$rows" -lt 501 ]; then
    echo "  (skip: too few rows)"
    rm -f "$out"
  fi
done < scripts/sp500.txt

# Benchmark SPY
echo "Downloading SPY benchmark"
curl -sL 'https://stooq.com/q/d/l/?s=spy.us&i=d' -o data/SPY.csv
