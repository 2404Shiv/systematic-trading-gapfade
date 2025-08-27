set -euo pipefail

LIST="${1:-scripts/test5.txt}"     # file containing tickers, one per line
OUTDIR="data"
mkdir -p "$OUTDIR"

UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36'
COOKIE=/tmp/yf.cookie
HTML=/tmp/yf.html
CRUMBF=/tmp/yf.crumb

get_crumb() {
  rm -f "$COOKIE" "$HTML" "$CRUMBF"
curl -s --compressed -A "$UA" -c "$COOKIE" \
       'https://finance.yahoo.com/quote/SPY/history?p=SPY' > "$HTML"

  # extract crumb from JSON; unescape slashes
  sed -n 's/.*"CrumbStore":{"crumb":"\([^"]*\)".*/\1/p' "$HTML" | sed 's/\\u002F/\//g' > "$CRUMBF"
  CRUMB=$(cat "$CRUMBF" 2>/dev/null || true)
  [[ -n "${CRUMB:-}" ]] || return 1
  echo "$CRUMB"
}
CRUMB=""
for t in 1 2 3 4 5; do
  if CRUMB=$(get_crumb); then break; fi
  sleep 2
done
if [[ -z "$CRUMB" ]]; then
  echo "Failed to get Yahoo crumb after retries" >&2
  exit 1
fi

NOW=$(date +%s)

fetch_one() {
  local sym="$1"
  [[ -z "$sym" ]] && return 0
  echo "Fetching $sym from Yahoo (AdjClose)…"

  local url="https://query1.finance.yahoo.com/v7/finance/download/${sym}?period1=0&period2=${NOW}&interval=1d&events=history&includeAdjustedClose=true&crumb=${CRUMB}"
  local raw="${OUTDIR}/${sym}_yh.csv"
  local out="${OUTDIR}/${sym}.csv"
if ! curl -s --compressed -A "$UA" -b "$COOKIE" -o "$raw" "$url"; then
    echo "WARN: download failed for $sym"; rm -f "$raw"; return 0
  fi

  # expect: Date,Open,High,Low,Close,Adj Close,Volume
  if head -n1 "$raw" | grep -q '^Date,Open,High,Low,Close,Adj Close,Volume'; then
    {
      echo "Date,Open,High,Low,Close,Volume"
      # factor f = AdjClose / Close; rebuild adjusted OHLC; set Close = Adj Close
      awk -F, 'NR>1 && $2!="" && $7!="" {
        f = ($5==0 ? 1 : $6/$5);
        printf "%s,%.6f,%.6f,%.6f,%.6f,%s\n",$1,$2*f,$3*f,$4*f,$6,$7
      }' "$raw"
    } > "${out}.tmp" && mv "${out}.tmp" "$out"
  else
    echo "WARN: $sym bad header / rate-limited; skipping"
    rm -f "$out"
  fi
  rm -f "$raw"
  sleep 1.5
}

while IFS= read -r sym; do
  fetch_one "$sym"
done < "$LIST"
