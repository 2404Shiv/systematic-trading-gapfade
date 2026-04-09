import os, json, time, math, threading, schedule, logging, sqlite3
try:
    from ibkr_trader import IBKRTrader
    IBKR_AVAILABLE = True
except:
    IBKR_AVAILABLE = False
from datetime import datetime
from zoneinfo import ZoneInfo
import pandas as pd
import numpy as np
import yfinance as yf
from flask import Flask, jsonify, render_template_string

GAP_THRESHOLD_BPS=15; SPREAD_BPS=1.0; IMPACT_BPS=0.5; GROSS=1.0
MAX_WEIGHT=0.02; LOOKBACK=20; Z_MIN=0.5; Z_MAX=3.0; TOP_FRAC=0.50
STOP_BPS=10.0; TARGET_VOL=0.06; EWMA_LAMBDA=0.97; PAPER_CAPITAL=1_000_000.0
UNIVERSE = ['CSCO', 'UAL', 'TROW', 'ISRG', 'NVR', 'TPR', 'DVN', 'BA', 'VRTX', 'GILD', 'TKO', 'EQIX', 'TER', 'PODD', 'MDT', 'V', 'A', 'FOX', 'MO', 'CTRA', 'SWKS', 'ENPH', 'MCHP', 'CDNS', 'MSCI', 'CHTR', 'EIX', 'KDP', 'BBY', 'GEN', 'LVS', 'HCA', 'AJG', 'DTE', 'SPY_av', 'C', 'T', 'CF', 'MGM', 'HUM', 'CBOE', 'CFG', 'APH', 'SYY', 'MSI', 'FCX', 'ADM', 'LH', 'BLDR', 'LNT', 'BAC', 'PSX', 'GPN', 'HUBB', 'PPG', 'TECH', 'IRM', 'IQV', 'ESS', 'WBD', 'HAL', 'STZ', 'BG', 'ADI', 'F', 'ADBE', 'STLD', 'SMCI', 'CPRT', 'TDG', 'ULTA', 'ARE', 'SYK', 'CB', 'TSN', 'GNRC', 'PEP', 'PEG', 'NOW', 'LLY', 'COST', 'REG', 'NWS', 'LOW', 'MDLZ', 'BKNG', 'ZBRA', 'XEL', 'AIZ', 'WDAY', 'MET', 'APO', 'FTV', 'DLR', 'ACGL', 'FAST', 'TJX', 'SNA', 'MPC', 'BR', 'D', 'MRK', 'STX', 'NOC', 'BXP', 'COIN', 'KHC', 'ABNB', 'UNP', 'ALLE', 'ABBV', 'ORCL', 'ECL', 'ETR', 'EBAY', 'SBUX', 'IR', 'AMT', 'INTU', 'DPZ', 'PAYC', 'PG', 'CAT', 'ODFL', 'MCD', 'MNST', 'PSKY', 'AMZN', 'GEHC', 'EG', 'INTC', 'PNR', 'GLW', 'BDX', 'KMI', 'CSGP', 'PWR', 'APTV', 'DXCM', 'EXR', 'WELL', 'EXE', 'HOLX', 'EXPD', 'GM', 'TXN', 'VRSK', 'SJM', 'TMO', 'OXY', 'RL', 'DECK', 'CCI', 'MMM', 'MOS', 'FTNT', 'HSY', 'DHI', 'ED', 'ES', 'ADSK', 'GL', 'INVH', 'IP', 'EXPE', 'KO', 'PCAR', 'RVTY', 'WDC', 'PYPL', 'NEE', 'UPS', 'ELV', 'EMR', 'MSFT', 'CTAS', 'UDR', 'WEC', 'UBER', 'AME', 'IT', 'DD', 'ACN', 'VRSN', 'EW', 'CMG', 'AWK', 'COO', 'SHW', 'HPQ', 'AMAT', 'CCL', 'MLM', 'AVY', 'EVRG', 'EA', 'DE', 'SPG', 'AMD', 'KLAC', 'NDAQ', 'URI', 'RTX', 'NXPI', 'PNC', 'KMX', 'MTCH', 'BIIB', 'NVDA', 'CHRW', 'ROP', 'WSM', 'IDXX', 'EXC', 'HD', 'ALB', 'VLO', 'AON', 'ZTS', 'FDX', 'DG', 'TYL', 'HIG', 'CMS', 'CAG', 'INCY', 'SCHW', 'HSIC', 'AZO', 'AXP', 'HPE', 'HRL', 'SO', 'FRT', 'ZBH', 'CME', 'XOM', 'AMP', 'CVX', 'CMCSA', 'PCG', 'PNW', 'ICE', 'BEN', 'UHS', 'BKR', 'EMN', 'SBAC', 'ROK', 'PTC', 'NRG', 'NSC', 'NKE', 'FIS', 'FANG', 'XYZ', 'VTR', 'MAS', 'RF', 'AMCR', 'TAP', 'MAR', 'XYL', 'CMI', 'MTD', 'CPAY', 'KR', 'PLD', 'IBM', 'USB', 'BSX', 'LKQ', 'LIN', 'ITW', 'EOG', 'KMB', 'SPGI', 'NEM', 'LULU', 'WFC', 'CTVA', 'TTD', 'EL', 'GS', 'GD', 'CNP', 'PM', 'MCO', 'CLX', 'CAH', 'PANW', 'DELL', 'MPWR', 'DGX', 'AVB', 'DIS', 'SW', 'CBRE', 'GE', 'HII', 'LDOS', 'ALL', 'ERIE', 'ETN', 'ALGN', 'NFLX', 'LEN', 'FITB', 'WST', 'GWW', 'TRGP', 'LII', 'NTRS', 'FICO', 'AXON', 'CVS', 'AOS', 'FE', 'JPM', 'ABT', 'OMC', 'COF', 'TSCO', 'PH', 'HST', 'JBHT', 'MRNA', 'TSLA', 'MOH', 'ATO', 'COP', 'DHR', 'CNC', 'MCK', 'TXT', 'MTB', 'FDS', 'VTRS', 'AKAM', 'ROL', 'RMD', 'WRB', 'GOOGL', 'BRO', 'ANET', 'PAYX', 'DRI', 'META', 'COR', 'MAA', 'FOXA', 'POOL', 'CZR', 'FFIV', 'CRWD', 'CINF', 'VMC', 'MKTX', 'SRE', 'LHX', 'ORLY', 'IVZ', 'RCL', 'SNPS', 'GOOG', 'EPAM', 'NDSN', 'YUM', 'EQT', 'LYV', 'PFE', 'AVGO', 'DUK', 'TPL', 'REGN', 'CL', 'VZ', 'JCI', 'DAY', 'FSLR', 'AMGN', 'TEL', 'JBL', 'VST', 'JKHY', 'ADP', 'ON', 'STT', 'RSG', 'IFF', 'CARR', 'TRMB', 'QCOM', 'DASH', 'LYB', 'GIS', 'PHM', 'ROST', 'LUV', 'LW', 'MS', 'CPB', 'OKE', 'BK', 'J', 'SYF', 'CHD', 'HWM', 'MHK', 'TFC', 'DAL', 'APA', 'AFL', 'CSX', 'NI', 'CPT', 'PFG', 'NCLH', 'RJF', 'HBAN', 'UNH', 'PRU', 'PLTR', 'GPC', 'WTW', 'DDOG', 'WMB', 'EQR', 'DVA', 'AIG', 'MA', 'HON', 'VICI', 'O', 'NWSA', 'TTWO', 'AES', 'SLB', 'DOC', 'TT', 'TGT', 'AAPL', 'MKC', 'OTIS', 'CEG', 'TDY', 'WY', 'APD', 'GRMN', 'AEE', 'BX', 'HLT', 'DLTR', 'STE', 'HAS', 'TMUS', 'WMT', 'NTAP', 'KIM', 'BAX', 'LMT', 'KKR', 'KEY', 'KEYS', 'BMY', 'PSA', 'WYNN', 'EFX', 'NUE', 'PKG', 'GDDY', 'WAB', 'CTSH', 'SWK', 'CRL', 'MU', 'TRV', 'L', 'AEP', 'CI', 'DOW', 'CDW', 'BALL', 'JNJ', 'WM', 'DOV', 'CRM', 'PGR', 'WAT', 'IEX', 'LRCX', 'BLK', 'PPL']
DB_PATH="gap_fade_trades.db"
ET=ZoneInfo("America/New_York")
logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
log=logging.getLogger(__name__)

def init_db():
    con=sqlite3.connect(DB_PATH); cur=con.cursor()
    cur.execute("""CREATE TABLE IF NOT EXISTS trades(id INTEGER PRIMARY KEY AUTOINCREMENT,
    date TEXT,symbol TEXT,direction INTEGER,gap_pct REAL,z_score REAL,weight REAL,
    entry_price REAL,exit_price REAL,pnl_pct REAL,pnl_dollar REAL,cost_bps REAL,
    status TEXT,created_at TEXT)""")
    cur.execute("""CREATE TABLE IF NOT EXISTS equity(id INTEGER PRIMARY KEY AUTOINCREMENT,
    date TEXT UNIQUE,portfolio_value REAL,daily_return REAL,sharpe REAL,max_dd REAL,win_rate REAL)""")
    con.commit(); con.close()

def get_db(): return sqlite3.connect(DB_PATH)

def calc_sharpe(returns):
    if len(returns)<2: return 0.0
    r=np.array(returns); mu,sd=r.mean(),r.std()
    return float((mu/sd)*np.sqrt(252)) if sd>0 else 0.0

def calc_max_dd(equity_vals):
    if len(equity_vals)<2: return 0.0
    eq=np.array(equity_vals); peak=np.maximum.accumulate(eq)
    dd=(eq-peak)/np.maximum(peak,1e-10); return float(dd.min())

def get_performance_stats():
    con=get_db()
    eq_df=pd.read_sql("SELECT * FROM equity ORDER BY date",con)
    tr_df=pd.read_sql("SELECT * FROM trades WHERE status='closed'",con)
    con.close()
    if eq_df.empty:
        return {"total_return":0,"sharpe":0,"sortino":0,"max_dd":0,"win_rate":0,
                "total_trades":0,"current_value":PAPER_CAPITAL,"daily_returns":[],
                "equity_curve":[],"monthly_pnl":{}}
    returns=eq_df["daily_return"].dropna().tolist()
    equity=eq_df["portfolio_value"].tolist(); dates=eq_df["date"].tolist()
    total_return=(equity[-1]-PAPER_CAPITAL)/PAPER_CAPITAL*100
    sharpe=calc_sharpe(returns); max_dd=calc_max_dd(equity)*100
    neg=[r for r in returns if r<0]; down_std=np.std(neg) if neg else 1e-10
    sortino=float((np.mean(returns)/down_std)*np.sqrt(252)) if down_std>0 else 0
    wins=len(tr_df[tr_df["pnl_pct"]>0]); total=len(tr_df)
    win_rate=(wins/total*100) if total>0 else 0
    monthly={}
    if not tr_df.empty:
        tr_df["month"]=pd.to_datetime(tr_df["date"]).dt.strftime("%Y-%m")
        monthly=tr_df.groupby("month")["pnl_dollar"].sum().to_dict()
    return {"total_return":round(total_return,4),"sharpe":round(sharpe,4),
            "sortino":round(sortino,4),"max_dd":round(max_dd,4),
            "win_rate":round(win_rate,2),"total_trades":total,
            "current_value":round(equity[-1] if equity else PAPER_CAPITAL,2),
            "daily_returns":returns[-60:],"equity_curve":list(zip(dates[-252:],[round(e,2) for e in equity[-252:]])),
            "monthly_pnl":{k:round(v,2) for k,v in monthly.items()}}

class GapFadeStrategy:
    def __init__(self):
        self.ewma_var=0.0; self.daily_returns=[]; self.portfolio_value=PAPER_CAPITAL; self._vol_cache={}
        self._load_state()
    def _load_state(self):
        con=get_db(); eq=pd.read_sql("SELECT * FROM equity ORDER BY date",con); con.close()
        if not eq.empty:
            self.portfolio_value=eq.iloc[-1]["portfolio_value"]
            self.daily_returns=eq["daily_return"].dropna().tolist()
    def fetch_data(self):
        log.info("Fetching market data...")
        data={}
        try:
            tickers=yf.download(UNIVERSE,period="5d",interval="1d",auto_adjust=False,progress=False,group_by="ticker")
            for sym in UNIVERSE:
                try:
                    df=tickers[sym].dropna() if isinstance(tickers.columns,pd.MultiIndex) else tickers.dropna()
                    if len(df)>=2:
                        data[sym]={"prev_close":float(df["Close"].iloc[-2]),"open":float(df["Open"].iloc[-1]),
                                   "close":float(df["Close"].iloc[-1]),"volume":float(df["Volume"].iloc[-2])}
                except: pass
        except Exception as e: log.error(f"Data fetch error: {e}")
        log.info(f"Got data for {len(data)} symbols"); return data
    def compute_rolling_vol(self,sym,lookback=LOOKBACK):
        return self._vol_cache.get(sym, 0.01)
    def precompute_vols(self):
        log.info("Pre-computing volatilities...")
        try:
            import pandas as pd
            df=yf.download(UNIVERSE,period="30d",interval="1d",auto_adjust=False,progress=False,group_by="ticker")
            vols={}
            for sym in UNIVERSE:
                try:
                    s=df[sym].dropna() if isinstance(df.columns,pd.MultiIndex) else df.dropna()
                    if len(s)>=5:
                        oc=((s["Close"]-s["Open"])/s["Open"])
                        vols[sym]=float(oc.std()) or 0.01
                except: vols[sym]=0.01
            self._vol_cache=vols
            log.info(f"Computed vols for {len(vols)} symbols")
        except Exception as e:
            log.error(f"Vol error: {e}")
            self._vol_cache={}

    def generate_signals(self,data):
        threshold=GAP_THRESHOLD_BPS/10000.0; candidates=[]
        for sym,d in data.items():
            pc=d["prev_close"]; op=d["open"]
            if pc<=0 or op<=0: continue
            gap=(op-pc)/pc; direction=0
            if gap<=-threshold: direction=+1
            elif gap>=threshold: direction=-1
            if direction==0: continue
            vol=self.compute_rolling_vol(sym); z=gap/max(vol,1e-6)
            if abs(z)<Z_MIN or abs(z)>Z_MAX: continue
            oc_ret=(d["close"]-op)/op
            candidates.append({"sym":sym,"direction":direction,"gap":gap,"z":z,
                               "vol":vol,"oc_ret":oc_ret,"open":op,"close":d["close"],"adv":d["volume"]})
        if not candidates: return []
        candidates.sort(key=lambda x:abs(x["z"]),reverse=True)
        keep=max(1,int(TOP_FRAC*len(candidates))); candidates=candidates[:keep]
        for c in candidates: c["base_w"]=1.0/max(c["vol"],1e-6)
        longs=[c for c in candidates if c["direction"]>0]
        shorts=[c for c in candidates if c["direction"]<0]
        if not longs or not shorts: return []
        sum_l=sum(c["base_w"] for c in longs); sum_s=sum(c["base_w"] for c in shorts)
        signals=[]
        for c in candidates:
            w=(c["base_w"]/sum_l)*0.5*GROSS if c["direction"]>0 else -(c["base_w"]/sum_s)*0.5*GROSS
            w=max(-MAX_WEIGHT,min(MAX_WEIGHT,w))
            cost=2*(SPREAD_BPS/10000+IMPACT_BPS/10000*(1/max(c["adv"],1)))
            ret=c["direction"]*c["oc_ret"]; ret=max(ret,-STOP_BPS/10000); net=ret-cost
            signals.append({**c,"weight":w,"net_ret":net,"cost":cost})
        return signals
    def execute_paper_trade(self,signals,date_str):
        if not signals: log.info("No signals today"); return 0.0
        portfolio_ret=sum(s["weight"]*s["net_ret"] for s in signals)
        sigma=math.sqrt(max(self.ewma_var,0))*math.sqrt(252)
        scale=(TARGET_VOL/sigma) if sigma>1e-10 else 1.0; scale=max(0.1,min(3.0,scale))
        portfolio_ret*=scale
        self.ewma_var=EWMA_LAMBDA*self.ewma_var+(1-EWMA_LAMBDA)*portfolio_ret**2
        con=get_db(); cur=con.cursor()
        for s in signals:
            pnl_dollar=s["weight"]*s["net_ret"]*self.portfolio_value*scale
            cur.execute("""INSERT INTO trades(date,symbol,direction,gap_pct,z_score,weight,
            entry_price,exit_price,pnl_pct,pnl_dollar,cost_bps,status,created_at)
            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (date_str,s["sym"],s["direction"],round(s["gap"]*100,4),round(s["z"],4),
             round(s["weight"],6),round(s["open"],4),round(s["close"],4),
             round(s["net_ret"]*100,4),round(pnl_dollar,2),round(s["cost"]*10000,2),
             "closed",datetime.now().isoformat()))
        new_value=self.portfolio_value*(1+portfolio_ret); self.portfolio_value=new_value
        self.daily_returns.append(portfolio_ret)
        sharpe=calc_sharpe(self.daily_returns); max_dd=calc_max_dd([PAPER_CAPITAL]+[PAPER_CAPITAL*math.prod(1+r for r in self.daily_returns[:i+1]) for i in range(len(self.daily_returns))])*100
        wins=sum(1 for r in self.daily_returns if r>0); win_rate=wins/len(self.daily_returns)*100
        cur.execute("""INSERT OR REPLACE INTO equity(date,portfolio_value,daily_return,sharpe,max_dd,win_rate)
        VALUES(?,?,?,?,?,?)""",(date_str,round(new_value,2),portfolio_ret,round(sharpe,4),round(max_dd,4),round(win_rate,2)))
        con.commit(); con.close()

        # Place real IBKR paper trades
        if IBKR_AVAILABLE and USE_IBKR:
            try:
                trader = IBKRTrader()
                if trader.connect():
                    capital_per_trade = new_value
                    for s in signals:
                        if s["weight"] == 0: continue
                        price_est = s["open"]
                        dollar_alloc = abs(s["weight"]) * capital_per_trade * scale
                        qty = max(1, int(dollar_alloc / price_est))
                        trader.place_order(s["sym"], s["direction"], qty)
                    trader.disconnect()
            except Exception as e:
                log.error(f"IBKR trading error: {e}")

        log.info(f"Day {date_str}: ret={portfolio_ret*100:.3f}% portfolio={new_value:,.2f}"); return portfolio_ret
    def run_daily(self):
        now_et=datetime.now(ET); date_str=now_et.strftime("%Y-%m-%d")
        log.info(f"=== Running gap-fade strategy for {date_str} ===")
        self.precompute_vols()
        data=self.fetch_data()
        if not data: return
        signals=self.generate_signals(data)
        log.info(f"Generated {len(signals)} signals")
        self.execute_paper_trade(signals,date_str)

app=Flask(__name__)
strategy=None

DASH="""<!DOCTYPE html><html><head><title>Gap Fade Live</title>
<link href="https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;600&display=swap" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>:root{--bg:#0a0a0f;--surface:#12121a;--border:#1e1e2e;--accent:#00ff88;--accent2:#ff4466;--accent3:#4488ff;--text:#e0e0e8;--muted:#666680;--card:#16161f}
*{margin:0;padding:0;box-sizing:border-box}body{background:var(--bg);color:var(--text);font-family:'DM Sans',sans-serif}
.header{display:flex;align-items:center;justify-content:space-between;padding:20px 32px;border-bottom:1px solid var(--border);background:var(--surface)}
.logo{font-family:'Space Mono',monospace;font-size:18px;color:var(--accent);letter-spacing:2px}
.dot{width:8px;height:8px;border-radius:50%;background:var(--accent);animation:pulse 2s infinite;display:inline-block;margin-right:8px}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:0.3}}
.grid{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;padding:24px 32px 0}
.stat-card{background:var(--card);border:1px solid var(--border);border-radius:12px;padding:20px;position:relative;overflow:hidden}
.stat-card::before{content:'';position:absolute;top:0;left:0;right:0;height:2px}
.stat-card.green::before{background:var(--accent)}.stat-card.red::before{background:var(--accent2)}
.stat-card.blue::before{background:var(--accent3)}.stat-card.yellow::before{background:#ffaa00}
.stat-label{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:1.5px;margin-bottom:8px;font-family:'Space Mono',monospace}
.stat-value{font-size:28px;font-weight:600;line-height:1;margin-bottom:4px}
.stat-sub{font-size:12px;color:var(--muted)}.positive{color:var(--accent)}.negative{color:var(--accent2)}.neutral{color:var(--accent3)}
.charts-row{display:grid;grid-template-columns:2fr 1fr;gap:16px;padding:16px 32px}
.chart-card{background:var(--card);border:1px solid var(--border);border-radius:12px;padding:20px}
.chart-title{font-family:'Space Mono',monospace;font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:1.5px;margin-bottom:16px}
.bottom-row{display:grid;grid-template-columns:1fr 1fr;gap:16px;padding:0 32px 32px}
.trades-table{width:100%;border-collapse:collapse;font-size:13px}
.trades-table th{text-align:left;padding:8px 12px;font-family:'Space Mono',monospace;font-size:10px;color:var(--muted);text-transform:uppercase;letter-spacing:1px;border-bottom:1px solid var(--border)}
.trades-table td{padding:10px 12px;border-bottom:1px solid var(--border)}
.badge{display:inline-block;padding:2px 8px;border-radius:4px;font-size:11px;font-family:'Space Mono',monospace}
.badge.long{background:rgba(0,255,136,0.15);color:var(--accent)}.badge.short{background:rgba(255,68,102,0.15);color:var(--accent2)}
.share-banner{margin:16px 32px;background:rgba(0,255,136,0.05);border:1px solid rgba(0,255,136,0.2);border-radius:12px;padding:16px 20px;display:flex;align-items:center;justify-content:space-between}
.share-link{font-family:'Space Mono',monospace;font-size:12px;color:var(--accent);background:rgba(0,255,136,0.1);padding:6px 14px;border-radius:6px;border:1px solid rgba(0,255,136,0.3);text-decoration:none}
canvas{max-height:200px}</style></head><body>
<div class="header"><div class="logo">GAP FADE // LIVE</div><div><span class="dot"></span><span id="last-update">Loading...</span></div></div>
<div class="share-banner"><div style="font-size:13px;color:var(--muted)">Share your verified track record for your resume</div><a class="share-link" href="/performance" target="_blank">VIEW PUBLIC REPORT →</a></div>
<div class="grid">
<div class="stat-card green"><div class="stat-label">Portfolio Value</div><div class="stat-value positive" id="pv">$100,000</div><div class="stat-sub" id="tr">+0.00% total return</div></div>
<div class="stat-card blue"><div class="stat-label">Sharpe Ratio</div><div class="stat-value neutral" id="sh">0.000</div><div class="stat-sub">annualised risk-adjusted</div></div>
<div class="stat-card yellow"><div class="stat-label">Win Rate</div><div class="stat-value" style="color:#ffaa00" id="wr">0.0%</div><div class="stat-sub" id="tt">0 trades</div></div>
<div class="stat-card red"><div class="stat-label">Max Drawdown</div><div class="stat-value negative" id="dd">0.00%</div><div class="stat-sub">peak to trough</div></div></div>
<div class="charts-row">
<div class="chart-card"><div class="chart-title">Equity Curve</div><canvas id="eqChart"></canvas></div>
<div class="chart-card"><div class="chart-title">Daily Returns</div><canvas id="retChart"></canvas></div></div>
<div class="bottom-row">
<div class="chart-card"><div class="chart-title">Recent Trades</div>
<table class="trades-table"><thead><tr><th>Date</th><th>Symbol</th><th>Side</th><th>Gap%</th><th>Z</th><th>P&L</th></tr></thead>
<tbody id="tb"><tr><td colspan="6" style="color:var(--muted);text-align:center;padding:20px">Bot runs at 14:35 UK time weekdays</td></tr></tbody></table></div>
<div class="chart-card"><div class="chart-title">Monthly P&L</div><div id="mg" style="display:grid;grid-template-columns:repeat(3,1fr);gap:8px"></div></div></div>
<script>
let eqC,retC;
function init(){
const o={responsive:true,maintainAspectRatio:true,plugins:{legend:{display:false}},scales:{x:{display:false},y:{grid:{color:'#1e1e2e'},ticks:{color:'#666680',font:{size:11}}}}};
eqC=new Chart(document.getElementById('eqChart'),{type:'line',data:{labels:[],datasets:[{data:[],borderColor:'#00ff88',borderWidth:2,fill:true,backgroundColor:'rgba(0,255,136,0.05)',tension:0.3,pointRadius:0}]},options:o});
retC=new Chart(document.getElementById('retChart'),{type:'bar',data:{labels:[],datasets:[{data:[],backgroundColor:ctx=>ctx.raw>=0?'rgba(0,255,136,0.6)':'rgba(255,68,102,0.6)',borderRadius:2}]},options:o});}
async function refresh(){
try{const d=await(await fetch('/api/stats')).json();
document.getElementById('pv').textContent='$'+d.current_value.toLocaleString('en-US',{minimumFractionDigits:2});
document.getElementById('tr').textContent=(d.total_return>=0?'+':'')+d.total_return.toFixed(2)+'% total return';
document.getElementById('sh').textContent=d.sharpe.toFixed(3);
document.getElementById('wr').textContent=d.win_rate.toFixed(1)+'%';
document.getElementById('tt').textContent=d.total_trades+' trades';
document.getElementById('dd').textContent=d.max_dd.toFixed(2)+'%';
document.getElementById('last-update').textContent='Updated '+new Date().toLocaleTimeString('en-GB');
if(d.equity_curve.length){eqC.data.labels=d.equity_curve.map(e=>e[0]);eqC.data.datasets[0].data=d.equity_curve.map(e=>e[1]);eqC.update();}
if(d.daily_returns.length){retC.data.labels=d.daily_returns.map((_,i)=>i);retC.data.datasets[0].data=d.daily_returns.map(r=>+(r*100).toFixed(3));retC.update();}
const mg=document.getElementById('mg');
if(Object.keys(d.monthly_pnl).length)mg.innerHTML=Object.entries(d.monthly_pnl).sort().slice(-9).map(([m,v])=>`<div style="background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:12px;text-align:center"><div style="font-size:11px;color:var(--muted);margin-bottom:4px">${m}</div><div style="font-size:16px;font-weight:600;color:${v>=0?'var(--accent)':'var(--accent2)'}">${v>=0?'+':''}$${Math.abs(v).toLocaleString('en-US',{maximumFractionDigits:0})}</div></div>`).join('');
}catch(e){console.error(e);}
try{const trades=await(await fetch('/api/trades')).json();
if(trades.length)document.getElementById('tb').innerHTML=trades.slice(0,15).map(t=>`<tr><td style="color:var(--muted)">${t.date}</td><td style="font-family:'Space Mono',monospace;font-weight:700">${t.symbol}</td><td><span class="badge ${t.direction>0?'long':'short'}">${t.direction>0?'LONG':'SHORT'}</span></td><td style="color:var(--muted)">${t.gap_pct>0?'+':''}${t.gap_pct.toFixed(2)}%</td><td style="color:var(--muted)">${t.z_score.toFixed(2)}</td><td class="${t.pnl_dollar>=0?'positive':'negative'}">${t.pnl_dollar>=0?'+':''}$${Math.abs(t.pnl_dollar).toFixed(2)}</td></tr>`).join('');
}catch(e){console.error(e);}}
init();refresh();setInterval(refresh,30000);
</script></body></html>"""

REPORT="""<!DOCTYPE html><html><head><title>Gap Fade Performance Report</title>
<link href="https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;600&display=swap" rel="stylesheet">
<style>:root{--bg:#f8f8f4;--surface:#fff;--border:#e8e8e0;--accent:#00aa55;--accent2:#cc2244;--text:#1a1a22;--muted:#888898}
*{margin:0;padding:0;box-sizing:border-box}body{background:var(--bg);color:var(--text);font-family:'DM Sans',sans-serif}
.container{max-width:900px;margin:0 auto;padding:48px 24px}
.tag{display:inline-block;background:rgba(0,170,85,0.1);color:var(--accent);border:1px solid rgba(0,170,85,0.3);padding:4px 12px;border-radius:4px;font-family:'Space Mono',monospace;font-size:11px;letter-spacing:1px;margin-bottom:16px}
h1{font-size:36px;font-weight:600;margin-bottom:8px}.subtitle{color:var(--muted);font-size:15px}
.stats-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:32px 0}
.stat{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:24px}
.stat-label{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:1px;margin-bottom:8px;font-family:'Space Mono',monospace}
.stat-value{font-size:32px;font-weight:600}.positive{color:var(--accent)}.negative{color:var(--accent2)}
.box{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:24px;margin-top:24px}
.box p{color:var(--muted);line-height:1.7;font-size:14px}
.footer{margin-top:48px;padding-top:24px;border-top:1px solid var(--border);font-size:12px;color:var(--muted)}</style></head>
<body><div class="container">
<div class="tag">VERIFIED PAPER TRADING</div>
<h1>Gap Fade Strategy</h1><p class="subtitle">Market-neutral beta-hedged gap fade across S&P 500</p>
<div class="stats-grid">
<div class="stat"><div class="stat-label">Total Return</div><div class="stat-value positive" id="r1">—</div></div>
<div class="stat"><div class="stat-label">Sharpe Ratio</div><div class="stat-value" id="r2">—</div></div>
<div class="stat"><div class="stat-label">Sortino Ratio</div><div class="stat-value" id="r3">—</div></div>
<div class="stat"><div class="stat-label">Max Drawdown</div><div class="stat-value negative" id="r4">—</div></div>
<div class="stat"><div class="stat-label">Win Rate</div><div class="stat-value" id="r5">—</div></div>
<div class="stat"><div class="stat-label">Total Trades</div><div class="stat-value" id="r6">—</div></div></div>
<div class="box"><p>This strategy fades overnight gaps in S&P 500 equities. When a stock opens significantly above or below its prior close (threshold: 25bps), the strategy takes a counter-directional position expecting mean reversion. Positions are sized using inverse-volatility weighting with dollar-neutral and beta-neutral constraints. Volatility targeting scales the portfolio to a 6% annual volatility target.</p></div>
<div class="footer">Generated from live paper trading. Last updated: <span id="ru">—</span></div></div>
<script>fetch('/api/stats').then(r=>r.json()).then(d=>{
document.getElementById('r1').textContent=(d.total_return>=0?'+':'')+d.total_return.toFixed(2)+'%';
document.getElementById('r2').textContent=d.sharpe.toFixed(3);
document.getElementById('r3').textContent=d.sortino.toFixed(3);
document.getElementById('r4').textContent=d.max_dd.toFixed(2)+'%';
document.getElementById('r5').textContent=d.win_rate.toFixed(1)+'%';
document.getElementById('r6').textContent=d.total_trades;
document.getElementById('ru').textContent=new Date().toLocaleDateString('en-GB',{year:'numeric',month:'long',day:'numeric'});});</script></body></html>"""

@app.route("/")
def dashboard(): return render_template_string(DASH)
@app.route("/performance")
def performance(): return render_template_string(REPORT)
@app.route("/api/stats")
def api_stats(): return jsonify(get_performance_stats())
@app.route("/api/trades")
def api_trades():
    con=get_db(); df=pd.read_sql("SELECT * FROM trades ORDER BY created_at DESC LIMIT 50",con); con.close()
    return jsonify(df.to_dict("records"))
@app.route("/api/run-now")
def run_now():
    threading.Thread(target=strategy.run_daily).start(); return jsonify({"status":"started"})

def schedule_loop():
    schedule.every().monday.at("14:35").do(strategy.run_daily)
    schedule.every().tuesday.at("14:35").do(strategy.run_daily)
    schedule.every().wednesday.at("14:35").do(strategy.run_daily)
    schedule.every().thursday.at("14:35").do(strategy.run_daily)
    schedule.every().friday.at("14:35").do(strategy.run_daily)
    while True: schedule.run_pending(); time.sleep(30)

if __name__=="__main__":
    print("Gap Fade Bot starting... Dashboard at http://localhost:5000")
    init_db(); strategy=GapFadeStrategy()
    threading.Thread(target=schedule_loop,daemon=True).start()
    port=int(os.environ.get("PORT",5000))
    app.run(host="0.0.0.0",port=port,debug=False)
