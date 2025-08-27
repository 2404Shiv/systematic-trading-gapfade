#!/usr/bin/env bash
set -euo pipefail

mkdir -p data
rm -f data/SPY.csv data/SPY_yahoo_raw.csv data/SPY_av.csv /tmp/yf.cookie /tmp/yf.crumb || true
now=$(date +%s)

have_good_csv() { [[ -s data/SPY.csv ]] && head -n1 data/SPY.csv | grep -q '^Date,'; }

# ---------- Try Yahoo (crumb + cookie) ----------
crumb="$(curl -sL -c /tmp/yf.cookie 'https://finance.yahoo.com/quote/SPY/history?p=SPY' \
  | sed -n 's/.*"CrumbStore":{"crumb":"\([^"]*\)".*/\1/p' | sed 's/\\u002F/\//g' || true)"
if [[ -n "${crumb:-}" ]]; then
  curl -fSL --retry 4 --retry-delay 2 -b /tmp/yf.cookie \
    "https://query1.finance.yahoo.com/v7/finance/download/SPY?period1=0&period2=${now}&interval=1d&events=history&includeAdjustedClose=true&crumb=${crumb}" \
    -o data/SPY_yahoo_raw.csv || true

  if [[ -s data/SPY_yahoo_raw.csv ]]; then
    # to: Date,Open,High,Low,Close,Volume (skip header, drop nulls, sort by date)
    awk -F, 'NR==1{next} $2!="" && $5!="" && $7!="" && $0 !~ /null/ {print $1","$2","$3","$4","$5","$7}' \
      data/SPY_yahoo_raw.csv | sort -t, -k1,1 > data/SPY.csv || true
  fi
fi

# ---------- AlphaVantage fallback (needs AVKEY) ----------
if ! have_good_csv && [[ -n "${AVKEY:-}" ]]; then
  curl -fSL "https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&symbol=SPY&datatype=csv&apikey=${AVKEY}" \
    -o data/SPY_av.csv || true
  if [[ -s data/SPY_av.csv ]]; then
    # AV csv: timestamp,open,high,low,close,adjusted_close,volume...
    awk -F, 'NR==1{next} {print $1","$2","$3","$4","$5","$6}' data/SPY_av.csv \
      | sort -t, -k1,1 > data/SPY.csv || true
  fi
fi

# ---------- Stooq fallback with retries and header-check ----------
if ! have_good_csv; then
  for n in 1 2 3; do
    curl -fSL 'https://stooq.com/q/d/l/?s=spy.us&i=d' -o data/SPY.csv || true
    if head -n1 data/SPY.csv | grep -q '^Date,'; then
      break
    fi
    echo "Stooq returned bad shape; retrying in 60s..." >&2
    rm -f data/SPY.csv
    sleep 60
  done
fi

# ---------- Final check ----------
if ! have_good_csv; then
  echo "Failed to build a valid data/SPY.csv" >&2
  exit 1
fi

echo "OK: data/SPY.csv ready"
