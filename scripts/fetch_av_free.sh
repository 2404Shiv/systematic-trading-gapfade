set -euo pipefail
list="${1:-scripts/test5.txt}"
: "${AVKEY:?Set AVKEY to your Alpha Vantage API key (export AVKEY=...)}"

mkdir -p data

fetch_one () {
  local sym="$1"
  local url="https://www.alphavantage.co/query?function=TIME_SERIES_DAILY&symbol=${sym}&outputsize=full&datatype=csv&apikey=${AVKEY}"
  local raw="data/${sym}_raw.csv"
  local out="data/${sym}.csv"

  echo "Fetching $sym (free DAILY)..."
  curl -sS "$url" -o "$raw" || { echo "curl failed for $sym"; return 1; }
if head -n1 "$raw" | grep -q '^timestamp,open,high,low,close,volume'; then
    {
      echo "Date,Open,High,Low,Close,Volume"
      # DAILY: timestamp,open,high,low,close,volume  -> fields: $1,$2,$3,$4,$5,$6
      awk -F, 'NR>1 && $2!="" {print $1","$2","$3","$4","$5","$6}' "$raw"
    } > "${out}.tmp" && mv "${out}.tmp" "$out"
  else
    echo "WARN: $sym returned non-CSV / rate-limited / error; skipping."
    rm -f "$out"
  fi
  rm -f "$raw"
}

while read -r sym; do
  [[ -z "$sym" ]] && continue
  fetch_one "$sym"
  sleep 12
done < "$list"
