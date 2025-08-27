set -euo pipefail
while IFS= read -r sym; do
  low=$(printf '%s' "$sym" | tr '[:upper:]' '[:lower:]')
  url="https://stooq.com/q/d/l/?s=${low}.us&i=d"
  out="data/${sym}.csv"

  if curl -fsSL "$url" -o "$out"; then
    # keep only sane files with a proper CSV header
    head -n1 "$out" | grep -q '^Date,' || rm -f "$out"
  else
    rm -f "$out"
  fi
sleep 0.2
done < scripts/sp500.txt
