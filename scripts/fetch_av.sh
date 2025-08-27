set -euo pipefail
list="${1:-scripts/test5.txt}"
: "${AVKEY:?Set AVKEY to your Alpha Vantage API key}"

mkdir -p data
while read -r sym; do
  [[ -z "$sym" ]] && continue
  out="data/$sym.csv"
  echo "Fetching $sym from Alpha Vantage..."
  curl -s "https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&symbol=${sym}&outputsize=full&datatype=csv&apikey=${AVKEY}" \
  | awk -F, 'NR==1{next} $2!="" && $5!="" {print $1","$2","$3","$4","$5","$6}' > "${out}.tmp" \
  && mv "${out}.tmp" "${out}"
  # free tier is 5 calls/min: sleep ≈12s to be safe
  sleep 12
done < "$list"
