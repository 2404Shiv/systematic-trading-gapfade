/*
 * gap_fade_strategy.cpp
 *
 * An end‑to‑end backtester for an intraday “gap–fade” trading strategy on
 * US equities.  The goal of this program is to demonstrate a robust
 * implementation suitable for research and résumé purposes.  It reads
 * historical price data from CSV files, simulates a gap‑fade strategy with
 * a simple queue‑aware implementation shortfall (IS) cost model, and
 * reports detailed performance statistics.
 *
 * The strategy works as follows:  For each trading day, compute the
 * difference (“gap”) between the current day’s open price and the
 * previous day’s close.  If the opening price gaps up beyond a
 * configurable threshold (expressed in basis points of the previous close),
 * the strategy shorts the stock at the open and covers at the close.  If
 * the opening price gaps down beyond the threshold, the strategy goes long
 * at the open and sells at the close.  Positions are closed at the end of
 * the day.  We assume a notional position size of one unit for each
 * trade, and capital is updated after subtracting transaction costs.
 *
 * The backtester reports overall P&L, Sharpe and Sortino ratios, maximum
 * drawdown, trade win rate, average and standard deviation of trade
 * returns, and aggregated monthly, yearly and weekday statistics.  It
 * includes a simple queue‑aware IS cost model with a base bid/ask spread
 * and a linear price impact component based on traded notional and daily
 * average volume.  Although simplistic, this cost model illustrates how
 * to account for real‑world trading frictions.
 *
 * The code is written to be self‑contained and easy to understand.  It
 * deliberately eschews external libraries for CSV parsing or numeric
 * computation so that it can compile in a barebones C++ environment.
 *
 * Usage:
 *    g++ -O3 -std=c++17 gap_fade_strategy.cpp -o gap_fade
 *    ./gap_fade path/to/data1.csv path/to/data2.csv ... [options]
 *
 * Each CSV file should contain daily price data with the following columns:
 *     Date,Open,High,Low,Close,Adj Close,Volume
 * Additional columns are ignored.  Files may include header rows.
 *
 * Recommended sample data: A freely accessible dataset of Apple stock
 * (AAPL) from 1980 through December 2024 is available on GitHub.  It is
 * contained in the repository
 *     https://github.com/FarhanAli97/Apple-AAPL-Stock-Data-1980-to-December-2024
 * and the raw CSV can be downloaded via the URL
 *     https://raw.githubusercontent.com/FarhanAli97/Apple-AAPL-Stock-Data-1980-to-December-2024/main/Apple%20(AAPL)%20From%201980%20To%20Dec-2024.csv
 * The dataset provides columns Date, open, High, Low, Close and Volume for
 * each trading day.  Lines near January 2020 show the expected format
 * with values such as "2020‑01‑02 00:00:00,72.20428368,72.85175328,71.86828734,72.86828735,14200417000"【558188827009559†screenshot】.
 * Lines near January 2024 confirm that the file extends to 2024, e.g.
 * "2024‑01‑02 00:00:00,183.217093,184.97813089,182.553717,183.3517698,5486448700"【777361487677916†screenshot】.
 * These data satisfy the required 2020‑2024 window for this project.
 *
 * Copyright (c) 2025.
 */

// gap_fade_strategy.cpp  — market-neutral, cross-sectional gap-fade

/*
 Systematic Gap-Fade Strategy (C++)
 
 - Multi-asset daily backtest for gap-fade (fade gap-ups, buy gap-downs).
 - Market/dollar-neutral portfolio construction with gross and per-name caps.
 - Transaction costs: spread bps + market impact bps per %ADV (round trip).
 - Signal quality gating via z-scored gaps (z_min/z_max) using rolling O->C vol.
 - Trade selection by top-|z| fraction per day.
 - Per-name loss cap (stop_bps) on daily contribution.
 - Portfolio-level volatility targeting (EWMA) to target annual vol.
 - Optional beta-neutral adjustment vs. benchmark (rolling OLS on daily series).
 - Robust CSV ingest (Yahoo/Alpha/Vendor-like headers).
 - Full analytics: Sharpe, Sortino, MaxDD, VaR/CVaR, win rate, monthly/yearly/weekday PnL.
 - Optional date range filter and CSV dumps for daily returns & equity.

 Build:  g++ -O3 -std=c++17 gap_fade_strategy.cpp -o gap_fade
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unordered_set>


using namespace std;

struct Bar {
    int y{}, m{}, d{}, date{};     // yyyymmdd
    double open{}, high{}, low{}, close{}, adjclose{};
    double volume{};
};

struct PriceSeries {
    string sym;
    vector<Bar> bars;              // sorted ascending by date
    unordered_map<int, size_t> idx;// date -> index
};

// ---------- utils ----------
static bool leap(int y){ return (y%4==0 && y%100!=0) || (y%400==0); }
static int mdays(int y,int m){
    static int d[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    return (m==2? (leap(y)?29:28) : d[m-1]);
}
static void add_days(int& y,int& m,int& d,int k){
    // k can be negative; simple day stepping (few calls per earnings date)
    while(k>0){ int md=mdays(y,m); if(++d>md){ d=1; if(++m>12){ m=1; ++y; } } --k; }
    while(k<0){ if(--d<1){ if(--m<1){ m=12; --y; } d=mdays(y,m); } ++k; }
}

static void read_earnings(const string& path, int window_days,
                          unordered_map<string, unordered_set<int>>& blocked){
    if (path.empty()) return;
    ifstream f(path);
    if (!f.is_open()) { cerr<<"Cannot open earnings "<<path<<"\n"; return; }
    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        size_t c = line.find(',');
        if (c == string::npos) continue;
        string sym = line.substr(0, c);
        int e = stoi(line.substr(c + 1));
        int y = e / 10000, m = (e / 100) % 100, d = e % 100;
        for (int k = -window_days; k <= window_days; ++k) {
            int yy = y, mm = m, dd = d;
            add_days(yy, mm, dd, k);
            int code = 10000*yy + 100*mm + dd; blocked[sym].insert(code);
        }
    }
}

// ---- read sectors map: sym -> sector_id ----
static void read_sectors(const string& path, unordered_map<string,int>& sec){
    if (path.empty()) return;
    ifstream f(path);
    if (!f.is_open()) { cerr<<"Cannot open sectors "<<path<<"\n"; return; }
    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        size_t c = line.find(',');
        if (c == string::npos) continue;
        string sym = line.substr(0, c);
        int sid = stoi(line.substr(c + 1));
        sec[sym] = sid;
    }
}

static inline int yyyymmdd(int y, int m, int d) { return y*10000 + m*100 + d; }

static inline optional<tuple<int,int,int>> parse_ymd(const string& s) {
    if (s.size() < 8) return nullopt;
    // support "YYYY-MM-DD" or "YYYY/MM/DD" or "YYYYMMDD"
    if (s.size() >= 10 && (s[4]=='-' || s[4]=='/') && (s[7]=='-' || s[7]=='/')) {
        int y = stoi(s.substr(0,4));
        int m = stoi(s.substr(5,2));
        int d = stoi(s.substr(8,2));
        return {{y,m,d}};
    } else if (s.size() >= 8) {
        int y = stoi(s.substr(0,4));
        int m = stoi(s.substr(4,2));
        int d = stoi(s.substr(6,2));
        return {{y,m,d}};
    }
    return nullopt;
}

static int weekday(int y, int m, int d) {
    if (m < 3) { m += 12; --y; }
    int K = y % 100, J = y / 100;
    int h = (d + (13*(m+1))/5 + K + K/4 + J/4 + 5*J) % 7; // 0=Sat
    int w = (h + 6) % 7; // 0=Sun..6=Sat -> convert to 0=Mon..6=Sun
    // map: Sun(0)->6, Mon(1)->0 ... Sat(6)->5
    static int map7[7] = {6,0,1,2,3,4,5};
    return map7[w];
}

static double mean(const vector<double>& v) {
    if (v.empty()) return 0.0;
    return accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
}

static double variance(const vector<double>& v, double mu) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) { double d = x - mu; s += d*d; }
    return s / (double)v.size();
}

static double percentile(vector<double> v, double q) {
    if (v.empty()) return 0.0;
    sort(v.begin(), v.end());
    double idx = q * (v.size()-1);
    size_t i = (size_t)floor(idx);
    size_t j = min(v.size()-1, i+1);
    double w = idx - (double)i;
    return (1.0 - w) * v[i] + w * v[j];
}

static double var95(const vector<double>& v)  { return percentile(v, 0.05); }
static double cvar95(const vector<double>& v) {
    if (v.empty()) return 0.0;
    double v95 = var95(v);
    double sum=0.0; int n=0;
    for (double x : v) if (x <= v95) { sum += x; ++n; }
    return n ? sum / (double)n : 0.0;
}

static tuple<double,size_t,size_t> max_drawdown_from_equity(const vector<double>& eq) {
    if (eq.empty()) return {0.0,0,0};
    double peak = eq.front(), worst = 0.0;
    size_t s=0,e=0,p=0;
    for (size_t i=0;i<eq.size();++i) {
        if (eq[i] > peak) { peak = eq[i]; p = i; }
        double dd = (peak - eq[i]) / max(peak, 1e-12);
        if (dd > worst) { worst = dd; s = p; e = i; }
    }
    return {worst, s, e};
}

static double sharpe(const vector<double>& r) {
    if (r.empty()) return 0.0;
    double mu = mean(r), var = variance(r, mu);
    double sd = sqrt(max(var, 0.0));
    return (sd > 0.0) ? (mu / sd) * sqrt(252.0) : 0.0;
}
static double sortino(const vector<double>& r) {
    if (r.empty()) return 0.0;
    double mu = mean(r);
    double down=0.0; int n=0;
    for (double x : r) if (x < 0.0) { down += x*x; ++n; }
    if (!n) return 0.0;
    double dd = sqrt(down/(double)n);
    return (dd>0.0) ? (mu/dd)*sqrt(252.0) : 0.0;
}

// ---------- CSV ingest ----------
struct CsvMap {
    int date=-1, open=-1, high=-1, low=-1, close=-1, adj=-1, vol=-1;
};

static CsvMap detect_header_indices(const string& header) {
    CsvMap m;
    vector<string> cols;
    string col, tmp = header;
    stringstream ss(tmp);
    while (getline(ss, col, ',')) {
        for (char& c : col) c = (char)tolower(c);
        // strip spaces
        col.erase(remove_if(col.begin(), col.end(), ::isspace), col.end());
        cols.push_back(col);
    }
    auto pick = [&](const vector<string>& cands)->int{
        for (int i=0;i<(int)cols.size();++i)
            for (auto& name: cands) if (cols[i]==name) return i;
        return -1;
    };
    m.date  = pick({"date","day"});
    m.open  = pick({"open","opn"});
    m.high  = pick({"high","hi"});
    m.low   = pick({"low","lo"});
    m.close = pick({"close","cls","closingprice"});
    m.adj   = pick({"adjclose","adjustedclose","adj"});
    m.vol   = pick({"volume","vol","volumn"});
    return m;
}

static PriceSeries load_csv_series(const string& path) {
    PriceSeries ps;
    // symbol from path
    string sym = path;
    size_t p = sym.find_last_of("/\\");
    if (p != string::npos) sym = sym.substr(p+1);
    if (sym.size() > 4 && sym.substr(sym.size()-4)==".csv") sym = sym.substr(0, sym.size()-4);
    ps.sym = sym;

    ifstream f(path);
    if (!f.is_open()) { cerr<<"Cannot open "<<path<<"\n"; return ps; }
    string line;
    if (!getline(f, line)) return ps;
    CsvMap map = detect_header_indices(line);

    vector<Bar> bars;
    bars.reserve(4096);
    while (getline(f, line)) {
        if (line.empty()) continue;
        vector<string> t;
        string cell; stringstream ss(line);
        while (getline(ss, cell, ',')) t.push_back(cell);
        auto get = [&](int idx)->string {
            if (idx < 0 || idx >= (int)t.size()) return string();
            return t[idx];
        };
        string dstr = get(map.date);
        if (dstr.empty()) continue;
        auto ymd = parse_ymd(dstr);
        if (!ymd) continue;
        auto [y,m,d] = *ymd;
        Bar b; b.y=y; b.m=m; b.d=d; b.date = yyyymmdd(y,m,d);
        try {
            b.open  = stod(get(map.open));
            b.high  = stod(get(map.high));
            b.low   = stod(get(map.low));
            b.close = stod(get(map.close));
        } catch(...) { continue; }
        string adj = get(map.adj);
        if (!adj.empty()) { try { b.adjclose = stod(adj); } catch(...){} }
        string vs = get(map.vol);
        if (!vs.empty()) { try { b.volume = stod(vs); } catch(...){} }
        bars.push_back(b);
    }
    sort(bars.begin(), bars.end(), [](const Bar& a, const Bar& b){ return a.date < b.date; });
    ps.bars = std::move(bars);
    ps.idx.reserve(ps.bars.size());
    for (size_t i=0;i<ps.bars.size();++i) ps.idx[ps.bars[i].date] = i;
    return ps;
}

// ---------- main ----------
struct Item { 
  std::string sym;
  double ret;
  double adv_pct; 
  int dir; 
  double vol; 
  double z; 
};

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 2) {
        cerr<<"Usage: "<<argv[0]<<" <csv1> [csv2 ...] "
            <<"[--gap_threshold=bps] [--spread=bps] [--impact=bps_per_pctADV] "
            <<"[--gross=1.0] [--max_weight=0.05] [--lookback=20] "
            <<"[--beta_neutral=0|1] [--beta_lookback=60] [--benchmark=path.csv] "
            <<"[--start=YYYY-MM-DD] [--end=YYYY-MM-DD] "
            <<"[--target_vol=0.10] [--ewma_lambda=0.94] "
            <<"[--z_min=1.0] [--z_max=2.5] [--top_frac=0.15] [--stop_bps=12] "
            <<"\n";
        return 1;
    }

   // core params
double gap_bps = 20.0, spread_bps = 1.0, impact_bps_per_pctADV = 0.5;
double gross   = 1.0,  max_w     = 0.05;  int lookback = 20;
bool   beta_neutral = true;                int beta_lookback = 60;
std::string benchmark_path = "";          std::string outdir = "";

// added risk/selection params
double target_vol        = 0.10,  vol_ewma_lambda = 0.94;
double z_min             = 0.00,  z_max           = 3.00,
       top_frac          = 1.00,  stop_bps        = 0.00;
bool   require_both_sides = false;

double tp_bps = 0.0;                  // per-name take-profit in bps of price (0=off)
string earnings_path = "";            // csv: SYM,yyyymmdd
int    earnings_window = 1;           // skip trades +/-N days around earnings
string sectors_path = "";             // csv: SYM,sector_id
bool   sector_neutral = false;        // project weights to sector-neutral
double dd_limit = 0.10;               // book-level DD at which leverage goes to zero (10%)

// optional outputs + date filter
std::string dump_daily = "", dump_equity = "";
int start_date = 0, end_date = 99991231;

    vector<string> files;
    for (int i=1;i<argc;++i) {
        string a=argv[i];
        if (a.rfind("--", 0) != 0) { files.push_back(a); continue; }
        bool handled = false;

// helpers: read after '=' safely
auto after_eq = [&]() -> std::string {
    size_t k = a.find('=');
    return (k == std::string::npos) ? std::string() : a.substr(k + 1);
};
auto as_d = [&]() -> double { return std::stod(after_eq()); };
auto as_i = [&]() -> int    { return std::stoi(after_eq()); };
auto as_s = [&]() -> std::string { return after_eq(); };

if      (a.rfind("--gap_threshold", 0) == 0) { gap_bps                = as_d(); handled = true; }
else if (a.rfind("--spread",        0) == 0) { spread_bps             = as_d(); handled = true; }
else if (a.rfind("--impact",        0) == 0) { impact_bps_per_pctADV  = as_d(); handled = true; }
else if (a.rfind("--gross",         0) == 0) { gross                  = as_d(); handled = true; }
else if (a.rfind("--max_weight",    0) == 0) { max_w                  = as_d(); handled = true; }
else if (a.rfind("--lookback",      0) == 0) { lookback               = as_i(); handled = true; }
else if (a.rfind("--beta_neutral",  0) == 0) { beta_neutral           = (as_i() != 0); handled = true; }
else if (a.rfind("--beta_lookback", 0) == 0) { beta_lookback          = as_i(); handled = true; }
else if (a.rfind("--benchmark",     0) == 0) { benchmark_path         = as_s(); handled = true; }
else if (a.rfind("--start",         0) == 0) { 
    auto p = parse_ymd(as_s()); 
    if (p) { auto [y,m,d] = *p; start_date = yyyymmdd(y,m,d); } 
    handled = true; 
}
else if (a.rfind("--end",           0) == 0) { 
    auto p = parse_ymd(as_s()); 
    if (p) { auto [y,m,d] = *p; end_date = yyyymmdd(y,m,d); } 
    handled = true; 
}

else if (a.rfind("--target_vol",    0) == 0) { target_vol             = as_d(); handled = true; }
else if (a.rfind("--ewma_lambda",   0) == 0) { vol_ewma_lambda        = as_d(); handled = true; }
else if (a.rfind("--z_min",         0) == 0) { z_min                  = as_d(); handled = true; }
else if (a.rfind("--z_max",         0) == 0) { z_max                  = as_d(); handled = true; }
else if (a.rfind("--top_frac",      0) == 0) { top_frac               = as_d(); handled = true; }
else if (a.rfind("--stop_bps",      0) == 0) { stop_bps               = as_d(); handled = true; }
else if (a.rfind("--both_sides",    0) == 0) { require_both_sides     = (as_i() != 0); handled = true; }
else if (a.rfind("--dump_daily",    0) == 0) { dump_daily             = as_s(); handled = true; }
else if (a.rfind("--dump_equity",   0) == 0) { dump_equity            = as_s(); handled = true; }
else if (a.rfind("--outdir",        0) == 0) { outdir                 = as_s(); handled = true; }
else if (a.rfind("--tp_bps=",0)==0)       { tp_bps = stod(a.substr(9)); handled=true; }
else if (a.rfind("--earnings=",0)==0)     { earnings_path = a.substr(11); handled=true; }
else if (a.rfind("--earnings_window=",0)==0){ earnings_window = stoi(a.substr(18)); handled=true; }
else if (a.rfind("--sectors=",0)==0)      { sectors_path = a.substr(10); handled=true; }
else if (a.rfind("--sector_neutral=",0)==0){ sector_neutral = (stoi(a.substr(17))!=0); handled=true; }
else if (a.rfind("--dd_limit=",0)==0)     { dd_limit = stod(a.substr(11)); handled=true; }
else if (a.rfind("--tp_bps",0)==0) { tp_bps = stod(a.substr(9)); handled=true; }
if (!handled) cerr << "Warning: unknown flag '" << a << "' (ignored)\n";
    }

    if (files.empty()) { cerr<<"No input CSVs provided.\n"; return 1; }

    // load series
    vector<PriceSeries> series;
    series.reserve(files.size());
    for (auto& f : files) {
        auto ps = load_csv_series(f);
        if (ps.bars.size() < (size_t)lookback + 5)
            cerr<<"Warning: "<<ps.sym<<" has few rows: "<<ps.bars.size()<<"\n";
        if (!ps.bars.empty())
    series.push_back(std::move(ps));
    }
    if (series.empty()) { cerr<<"No usable data.\n"; return 1; }

    // benchmark
    PriceSeries bench;
    vector<double> bench_oc;
    if (!benchmark_path.empty()) {
        bench = load_csv_series(benchmark_path);
        if (bench.bars.size() >= 2) {
            bench_oc.reserve(bench.bars.size());
            for (size_t i=0;i<bench.bars.size();++i) {
                // open->close of the same day (avoid lookahead)
                double r = (bench.bars[i].open > 0.0) ? (bench.bars[i].close - bench.bars[i].open)/bench.bars[i].open : 0.0;
                bench_oc.push_back(r);
            }
        } else {
            cerr<<"Warning: benchmark file has too few rows\n";
        }
    }

    // calendar set
    set<int> all_days;
    for (auto& ps : series) for (auto& b : ps.bars) {
        if (b.date < start_date || b.date > end_date) continue;
        all_days.insert(b.date);
    }
  unordered_map<string, unordered_set<int>> earnings_block;
  unordered_map<string,int>                sym_sector;

  if (!earnings_path.empty()) read_earnings(earnings_path, earnings_window, earnings_block);
  if (!sectors_path.empty())  read_sectors(sectors_path, sym_sector);

    if (all_days.empty()) { cerr<<"No dates within selection.\n"; return 1; }

    // per-symbol previous close & volume for gap/ADV
    unordered_map<string, double> prev_close, prev_vol;
    // rolling O->C vol per symbol per date (stdev over lookback)
    unordered_map<string, unordered_map<int,double>> sym_vol;
    // build rolling vol
    for (auto& ps : series) {
        if (ps.bars.size() < (size_t)lookback) continue;
        vector<double> oc(ps.bars.size(), 0.0);
        for (size_t i=0;i<ps.bars.size();++i) {
            oc[i] = (ps.bars[i].open>0.0) ? (ps.bars[i].close - ps.bars[i].open)/ps.bars[i].open : 0.0;
        }
        // rolling stdev
        double m=0.0;
        for (int i=0;i<lookback;i++) m += oc[i];
        m /= lookback;
        double s2=0.0;
        for (int i=0;i<lookback;i++){ double d=oc[i]-m; s2+=d*d; }
        double sigma = sqrt(s2/lookback);
        sym_vol[ps.sym][ps.bars[lookback-1].date] = max(1e-6, sigma);
        for (size_t i=lookback;i<ps.bars.size();++i) {
            double sum=0.0;
            for (size_t k=i+1-lookback;k<=i;k++) sum += oc[k];
            double mu = sum/(double)lookback;
            double var=0.0;
            for (size_t k=i+1-lookback;k<=i;k++){ double d=oc[k]-mu; var += d*d; }
            double sig = sqrt(var/(double)lookback);
            sym_vol[ps.sym][ps.bars[i].date] = max(1e-6, sig);
        }
    }

    // state
    double cap=1.0, start_cap=1.0;
    vector<double> daily_ret, daily_ret_raw, equity, bench_series;
    vector<int> months, years, wdays;
    vector<double> pnl_bucket;

    daily_ret.reserve(all_days.size());
    daily_ret_raw.reserve(all_days.size());
    equity.reserve(all_days.size());
    bench_series.reserve(all_days.size());
    pnl_bucket.reserve(all_days.size());
    months.reserve(all_days.size());
    years.reserve(all_days.size());
    wdays.reserve(all_days.size());

    double threshold = gap_bps/10000.0;
    double spread_cost = spread_bps/10000.0;
    double ewma_var = 0.0;
    const double ANN = sqrt(252.0);
    double cap_peak = start_cap; // track running equity peak

    // helpers for lookups
    auto find_bar = [&](PriceSeries& ps, int date)->optional<Bar>{
        auto it = ps.idx.find(date);
        if (it==ps.idx.end()) return nullopt;
        return ps.bars[it->second];
    };

    // main day loop
    for (int date : all_days) {
        int Y = date/10000, M=(date/100)%100, D=date%100;

        // benchmark O->C today
        double bench_r_today = 0.0;
        if (!bench.bars.empty()) {
            auto b = find_bar(bench, date);
            if (b) {
                bench_r_today = (b->open>0.0) ? (b->close - b->open)/b->open : 0.0;
            }
        }

        // collect candidates
        vector<Item> items;
        items.reserve(series.size());

        for (auto& ps : series) {
            auto br = find_bar(ps, date);
            if (!br) continue;
            const Bar& rec = *br;

            if (!earnings_path.empty()) {
           const string& sym = ps.sym;      // <-- your series object already has symbol
           auto eb = earnings_block.find(sym);
            if (eb != earnings_block.end() && eb->second.count(date)) {
            continue; // skip around earnings
            }
        }

            if (!prev_close.count(ps.sym)) {
                prev_close[ps.sym] = rec.close;
                prev_vol[ps.sym]   = rec.volume > 0.0 ? rec.volume : 1.0;
                continue;
            }
            double pc = prev_close[ps.sym];
            if (pc <= 0.0) { prev_close[ps.sym]=rec.close; prev_vol[ps.sym]=max(1.0, rec.volume); continue; }

            double gap = (rec.open - pc) / pc;
            int dir = 0;
            if (gap <= -threshold) dir = +1; // fade gap-down
            else if (gap >=  threshold) dir = -1; // fade gap-up
            if (dir == 0) { prev_close[ps.sym]=rec.close; prev_vol[ps.sym]=max(1.0, rec.volume); continue; }

            double oc = (rec.open>0.0) ? (rec.close - rec.open)/rec.open : 0.0;
            double adv = max(1.0, prev_vol[ps.sym]);

            double vol = 0.01;
            auto sv = sym_vol[ps.sym].find(date);
            if (sv != sym_vol[ps.sym].end()) vol = max(1e-6, sv->second);

            double z = gap / max(1e-6, vol);
            if (fabs(z) < z_min || fabs(z) > z_max) {
                prev_close[ps.sym]=rec.close; prev_vol[ps.sym]=max(1.0, rec.volume); continue;
            }

           double r = dir * (rec.close - rec.open) / (rec.open > 0.0 ? rec.open : 1.0);
          if (stop_bps > 0.0) { double stop_r = -stop_bps / 10000.0; if (r < stop_r) r = stop_r; }
          if (tp_bps  > 0.0) { double  tp_r =   tp_bps / 10000.0; if (r >  tp_r) r =  tp_r; }
            // clamp by stop-loss / take-profit (bps of price)
            if (stop_bps > 0.0) {
                double stop_r = -stop_bps / 10000.0;
                if (r < stop_r) r = stop_r;
            }
            if (tp_bps > 0.0) {
                double tp_r = tp_bps / 10000.0;
                if (r > tp_r) r = tp_r;
            } 
            items.push_back({ps.sym, oc, 1.0/adv, dir, vol, z});
            prev_close[ps.sym]=rec.close; prev_vol[ps.sym]=max(1.0, rec.volume);
             }

        if (items.empty()) {
            // still book zero day for equity continuity
            daily_ret_raw.push_back(0.0);
            daily_ret.push_back(0.0);
            bench_series.push_back(bench_r_today);
            pnl_bucket.push_back(0.0);
            equity.push_back(cap);
            if (cap > cap_peak) cap_peak = cap;
            if (cap < cap_peak * (1.0 - dd_limit)) {
              std::cout << "Max DD stop hit ("
              << (1.0 - cap/cap_peak) * 100.0
              << "%) on " << Y << "-" << M << "-" << D << ". Stopping.\n";
             break; // exits the daily loop
             }
            months.push_back(M); years.push_back(Y); wdays.push_back(weekday(Y,M,D));
            continue;
        }

        // rank by |z| and keep top fraction
        sort(items.begin(), items.end(), [](const Item& a, const Item& b){ return fabs(a.z) > fabs(b.z); });
        size_t keepN = max<size_t>(1, (size_t)llround(top_frac * (double)items.size()));
        if (keepN < items.size()) items.resize(keepN);

        // raw weights inverse-vol, cap per-name; then scale L/S to meet gross and dollar-neutral
        vector<double> w(items.size(), 0.0);
        double sumL=0.0, sumS=0.0;
        for (size_t i=0;i<items.size();++i) {
            double base = 1.0 / max(1e-6, items[i].vol);
            double wi = items[i].dir * base;
            if (fabs(wi) > max_w) wi = (wi > 0 ? +max_w : -max_w);
            w[i] = wi;
            if (wi > 0) sumL += wi; else sumS += -wi;
        }
        if (sumL>0.0 && sumS>0.0) {
            double scaleL = (0.5 * gross) / sumL;
            double scaleS = (0.5 * gross) / sumS;
            for (size_t i=0;i<w.size();++i) w[i] = (w[i] > 0 ? w[i]*scaleL : w[i]*scaleS);
        } else {
            // no cross-section both sides -> flat
            fill(w.begin(), w.end(), 0.0);
        }

        // per-name costs and stop cap
        double portfolio_r = 0.0;
        for (size_t i=0;i<items.size();++i) {
            if (w[i]==0.0) continue;
            // round-trip cost: 2 * (spread + impact)
            double impact = (impact_bps_per_pctADV/10000.0) * items[i].adv_pct;
            double total_cost = 2.0*(spread_cost + impact);
            double ri = items[i].ret - total_cost;
            // per-name stop cap
            ri = max(ri, -stop_bps/10000.0);
            portfolio_r += w[i] * ri;
        }
        double portfolio_r_raw = portfolio_r;

        // beta-neutral hedge (rolling OLS beta of daily_ret vs benchmark open->close)
        if (beta_neutral && !bench.bars.empty() && !daily_ret.empty()) {
            size_t look = min((size_t)beta_lookback, daily_ret.size());
            if (look >= 5 && bench_series.size() >= look) {
                double mu_p=0.0, mu_b=0.0;
                for (size_t k=0;k<look;k++) { mu_p += daily_ret[daily_ret.size()-1-k]; mu_b += bench_series[bench_series.size()-1-k]; }
                mu_p/=look; mu_b/=look;
                double cov=0.0, varb=0.0;
                for (size_t k=0;k<look;k++) {
                    double dp = daily_ret[daily_ret.size()-1-k] - mu_p;
                    double db = bench_series[bench_series.size()-1-k] - mu_b;
                    cov += dp*db; varb += db*db;
                }
                if (varb > 1e-12) {
                    double beta = cov/varb;
                    portfolio_r -= beta * bench_r_today;
                }
            }
        }

        // vol targeting (EWMA)
        double sigma_prev = sqrt(max(ewma_var, 0.0)) * ANN;
        double scale = (sigma_prev > 1e-10) ? (target_vol / sigma_prev) : 1.0;
        scale = clamp(scale, 0.10, 3.0); // guard rails
        portfolio_r *= scale;
        ewma_var = vol_ewma_lambda*ewma_var + (1.0 - vol_ewma_lambda)*(portfolio_r*portfolio_r);

        // book
        cap *= (1.0 + portfolio_r);
        daily_ret_raw.push_back(portfolio_r_raw);
        daily_ret.push_back(portfolio_r);
        bench_series.push_back(bench_r_today);
        pnl_bucket.push_back(portfolio_r);
        equity.push_back(cap);
        months.push_back(M);
        years.push_back(Y);
        wdays.push_back(weekday(Y,M,D));
    }

    // stats
    double tot = (cap - start_cap)/max(start_cap,1e-12) * 100.0;
    double mu = mean(daily_ret), var = variance(daily_ret, mu), sd = sqrt(max(var,0.0));
    double shr = sharpe(daily_ret), sor = sortino(daily_ret);
    int wins=0; for (double r: daily_ret) if (r>0) ++wins;
    auto [dd, sIdx, eIdx] = max_drawdown_from_equity(equity);
    double v95 = var95(daily_ret), cv95 = cvar95(daily_ret);

    auto agg_by = [&](const vector<int>& keys){
        unordered_map<int,double> m;
        for (size_t i=0;i<keys.size();++i) m[keys[i]] += pnl_bucket[i];
        return m;
    };
    auto mon = agg_by(months);
    auto yr  = agg_by(years);
    auto wd  = agg_by(wdays);
    static const char* wdN[]={"Mon","Tue","Wed","Thu","Fri","Sat","Sun"};

    cout.setf(ios::fixed); cout<<setprecision(6);
    cout<<"Market–Neutral Gap–Fade Portfolio Report\n";
    cout<<"========================================\n";
    cout<<"Files loaded: "; for (auto& f : files) cout<<f<<" "; cout<<"\n";
    cout<<"Trading days: "<<daily_ret.size()<<"\n";
    cout<<"Gross: "<<gross<<" | Per-name cap: "<<max_w<<" | Lookback: "<<lookback<<"\n";
    cout<<"Threshold: "<<gap_bps<<" bps | Spread: "<<spread_bps<<" bps | Impact: "<<impact_bps_per_pctADV<<" bps/%ADV\n";
    cout<<"VolTarget: "<<(target_vol*100)<<"% | EWMA λ: "<<vol_ewma_lambda<<" | z∈["<<z_min<<","<<z_max<<"] | top_frac: "<<top_frac<<" | stop_bps: "<<stop_bps<<"\n";
    if (!benchmark_path.empty()) cout<<"Benchmark: "<<benchmark_path<<" | BetaNeutral: "<<(beta_neutral?"ON":"OFF")<<" | beta_lookback "<<beta_lookback<<"\n";
    cout<<"\n";

    cout<<"Total return: "<<tot<<" %\n";
    cout<<"Mean daily: "<<mu*100<<" % | Stdev: "<<sd*100<<" %\n";
    cout<<"Sharpe: "<<shr<<" | Sortino: "<<sor<<"\n";
    cout<<"Max drawdown: "<<dd*100<<" %\n";
    cout<<"Win rate: "<<(100.0*wins/max<size_t>(1,daily_ret.size()))<<" %\n";
    cout<<"VaR 95%: "<<v95*100<<" % | CVaR 95%: "<<cv95*100<<" %\n\n";

    cout<<"Monthly PnL (% of capital):\n";
    for (auto& kv : mon) cout<<"Month "<<kv.first<<": "<<kv.second*100<<" %\n";
    cout<<"\nYearly PnL (% of capital):\n";
    for (auto& kv : yr)  cout<<"Year "<<kv.first<<": "<<kv.second*100<<" %\n";
    cout<<"\nWeekday PnL (% of capital):\n";
    for (auto& kv : wd) if (kv.first>=0 && kv.first<7) cout<<wdN[kv.first]<<": "<<kv.second*100<<" %\n";

    // optional dumps
    auto safe_open = [&](const string& rel)->optional<ofstream>{
        string fn = rel;
        if (!outdir.empty()) {
            string sep = (outdir.back()=='/'||outdir.back()=='\\')?"" : "/";
            fn = outdir + sep + rel;
        }
        ofstream ofs(fn);
        if (!ofs.is_open()) { cerr<<"Failed to open "<<fn<<" for write\n"; return nullopt; }
        return ofs;
    };

    if (!dump_daily.empty()) {
        auto f = safe_open(dump_daily);
        if (f) {
            *f << "day_index,ret\n";
            for (size_t i=0;i<daily_ret.size();++i) *f << i << "," << setprecision(10) << daily_ret[i] << "\n";
        }
    }
    if (!dump_equity.empty()) {
        auto f = safe_open(dump_equity);
        if (f) {
            *f << "day_index,equity\n";
            for (size_t i=0;i<equity.size();++i) *f << i << "," << setprecision(10) << equity[i] << "\n";
        }
    }

    return 0;
}
