#!/usr/bin/env zsh
set -euo pipefail
mkdir -p runs
out="runs/sweep_results.csv"
print "z_min,top_frac,thr_bps,total_ret,sharpe,max_dd,win_rate" > $out

for zmin in 0.6 0.8 1.0; do
  for top in 0.40 0.50 0.60; do
    for thr in 20 25 30 35; do
      echo "=== z_min=${zmin} top_frac=${top} gap_threshold=${thr} ==="
      ./gap_fade data/AAPL.csv data/MSFT.csv data/AMZN.csv data/GOOGL.csv \
        --gap_threshold=$thr --spread=1 --impact=0.5 --gross=1 --lookback=20 \
        --benchmark=data/SPY.csv --beta_neutral=1 --beta_lookback=80 \
        --z_min=$zmin --z_max=3.0 --top_frac=$top \
        --target_vol=0.06 --ewma_lambda=0.97 \
        --stop_bps=10 --both_sides=0 --max_weight=0.02 \
      | awk -v z=$zmin -v t=$top -v thr=$thr '
          /^Total return:/ {ret=$3}
          /^Sharpe:/      {sh=$2}
          /^Max drawdown:/ {dd=$3}
          /^Win rate:/    {wr=$3}
          END{
            gsub("%","",ret); gsub("%","",dd); gsub("%","",wr);
            print z","t","thr","ret","sh","dd","wr
          }' >> $out
    done
  done
done

echo "\nTop by Sharpe with DD<10% and Win≥30%:"
awk -F, 'NR>1 && $6+0<10 && $7+0>=30 {print}' $out | sort -t, -k5,5nr | head -20
