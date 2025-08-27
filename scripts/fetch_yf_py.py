#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd

try:
    import yfinance as yf
except ImportError:
    sys.exit("Please: pip install yfinance")

# ---------- inputs ----------
arg = sys.argv[1] if len(sys.argv) > 1 else "scripts/test5.txt"
outdir = Path("data")
outdir.mkdir(exist_ok=True)

# read tickers: file with one per line OR a single symbol passed directly
symbols: list[str]
p = Path(arg)
if p.exists():
    with p.open() as f:
        symbols = [ln.strip().upper() for ln in f if ln.strip() and not ln.startswith("#")]
else:
    symbols = [arg.strip().upper()]

# Known bad pre-inception rows (Yahoo quirk)
BAD_STARTS = {"SPY": "1993-01-29"}  # use 1993-01-29 for SPY to be safe

for t in symbols:
    sym = t.upper()
    print(f"Fetching {sym} via yfinance …")
    try:
       df = yf.download(
    sym,
    period="max",          # <= add this (or start="1980-01-01")
    interval="1d",
    auto_adjust=False,
    progress=False,
    actions=False,
    group_by="column",
    threads=False
)
    except Exception as e:
        print(f"  ERROR {sym}: {e}")
        continue

    if df is None or df.empty:
        print(f"  WARN {sym}: empty dataframe")
        continue

    # flatten MultiIndex if returned (paranoia)
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = df.columns.get_level_values(0)

    required = ["Open", "High", "Low", "Close", "Adj Close", "Volume"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"  WARN {sym}: missing columns {missing} – skipping")
        continue

    # drop very early bogus rows (e.g. SPY pre-inception)
    if sym in BAD_STARTS:
        df = df.loc[df.index >= BAD_STARTS[sym]]

    # require positive OHLC/Adj Close and Volume
    valid = (df[["Open","High","Low","Close","Adj Close"]].to_numpy() > 0).all(axis=1) & (df["Volume"] > 0).to_numpy()
    df = df.loc[valid]
    if df.empty:
        print(f"  WARN {sym}: no valid rows after cleaning")
        continue

    # adjust OHLC to Adj Close; set Close := Adj Close
    factor = (df["Adj Close"] / df["Close"]).replace([np.inf, -np.inf, np.nan], 1.0).fillna(1.0)
    out = pd.DataFrame({
        "Date": df.index.to_series().dt.strftime("%Y-%m-%d"),
        "Open":  (df["Open"]  * factor).to_numpy(),
        "High":  (df["High"]  * factor).to_numpy(),
        "Low":   (df["Low"]   * factor).to_numpy(),
        "Close": df["Adj Close"].to_numpy(),
        "Volume": pd.to_numeric(df["Volume"], errors="coerce").fillna(0).round().astype("int64").to_numpy(),
    })
    # nice rounding
    out[["Open","High","Low","Close"]] = out[["Open","High","Low","Close"]].round(6)

    out_path = outdir / f"{sym}.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path}  rows={len(out)}")

print("Done.")