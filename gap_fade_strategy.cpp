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

// gap_fade_strategy.cpp  — market-neutral, cross-sectional gap-fade
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <deque>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::set;
using std::string;
using std::stringstream;
using std::unordered_map;
using std::vector;

struct PriceRecord {
    int y{0}, m{0}, d{0};
    double open{0}, high{0}, low{0}, close{0};
    long long volume{0};
};

static bool parse_date(const string& s, int& y, int& m, int& d) {
    if (s.size() < 10) return false;
    y = std::stoi(s.substr(0, 4));
    m = std::stoi(s.substr(5, 2));
    d = std::stoi(s.substr(8, 2));
    return true;
}

static string basename_noext(const string& path) {
    size_t sl = path.find_last_of("/\\");
    string base = (sl == string::npos) ? path : path.substr(sl + 1);
    size_t dot = base.find_last_of('.');
    if (dot != string::npos) base = base.substr(0, dot);
    return base;
}

// Read CSV (Date,Open,High,Low,Close,(Adj),Volume), ignore header if present.
static void read_csv(const string& filename, vector<PriceRecord>& out) {
    ifstream f(filename);
    if (!f.is_open()) { cerr << "Cannot open " << filename << endl; return; }
    string line;
    bool first=true;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (first) {
            first=false;
            // skip header if first token not starting with digit
            if (!line.empty() && !isdigit(line[0])) continue;
        }
        stringstream ss(line);
        string ds, op, hi, lo, cl, adj, vol;
        if (!std::getline(ss, ds, ',')) continue;
        std::getline(ss, op, ','); std::getline(ss, hi, ',');
        std::getline(ss, lo, ','); std::getline(ss, cl, ',');
        std::getline(ss, adj, ','); std::getline(ss, vol, ',');
        PriceRecord r;
        if (!parse_date(ds, r.y, r.m, r.d)) continue;
        try {
            r.open = std::stod(op); r.high=std::stod(hi);
            r.low = std::stod(lo);  r.close=std::stod(cl);
            r.volume = vol.empty()?0:std::stoll(vol);
        } catch (...) { continue; }
        out.push_back(r);
    }
}

static double sharpe(const vector<double>& r) {
    if (r.empty()) return 0;
    double mu=0; for (double x: r) mu+=x; mu/=r.size();
    double v=0; for (double x: r){ double d=x-mu; v+=d*d; } v/=r.size();
    double sd = std::sqrt(v);
    return (sd==0)?0:(mu/sd)*std::sqrt(252.0);
}
static double sortino(const vector<double>& r) {
    if (r.empty()) return 0;
    double mu=0; for (double x: r) mu+=x; mu/=r.size();
    double down=0; int n=0;
    for (double x: r) if (x<0){ down+=x*x; ++n; }
    if (n==0) return 0;
    double dd = std::sqrt(down/n);
    return (dd==0)?0:(mu/dd)*std::sqrt(252.0);
}
static std::tuple<double,size_t,size_t> max_dd(const vector<double>& eq){
    if (eq.empty()) return {0,0,0};
    double peak=eq[0], mdd=0; size_t s=0,e=0, p=0;
    for (size_t i=1;i<eq.size();++i){
        if (eq[i]>peak){ peak=eq[i]; p=i; }
        double d=(eq[i]-peak)/peak;
        if (d<mdd){ mdd=d; s=p; e=i; }
    }
    return {mdd,s,e};
}
static double percentile(vector<double> v, double q){
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    double pos = q*(v.size()-1);
    size_t i = (size_t)std::floor(pos);
    size_t j = std::min(v.size()-1, i+1);
    double w = pos - i;
    return v[i]*(1-w) + v[j]*w;
}
static double var95(const vector<double>& r){
    // losses -> left tail; VaR95 ~ 5th percentile of returns
    return percentile(r, 0.05);
}
static double cvar95(const vector<double>& r){
    if (r.empty()) return 0;
    double v = var95(r);
    double sum=0; int n=0;
    for (double x: r) if (x<=v){ sum+=x; ++n; }
    return (n? sum/n : v);
}
static int yyyymmdd(int y,int m,int d){ return y*10000 + m*100 + d; }
static int weekday(int y,int m,int d){
    if(m<3){ m+=12; y-=1; }
    int K=y%100, J=y/100;
    int h=(d + 13*(m+1)/5 + K + K/4 + J/4 + 5*J)%7; // 0=Sat
    return (h+5)%7; // 0=Mon
}

struct Metrics {
    double total_return{0};
    double mean{0}, stdev{0}, sharpe{0}, sortino{0};
    double max_drawdown{0};
    double VaR95{0}, CVaR95{0};
    double win_rate{0};
};

struct SeriesRec{ int y,m,d; double r; };

int main(int argc, char** argv){
    if (argc<2){
        cerr<<"Usage: "<<argv[0]<<" <csv1> [csv2 ...] "
            <<"[--gap_threshold=bps] [--spread=bps] [--impact=bps_per_pctADV] "
            <<"[--gross=1.0] [--max_weight=0.05] [--lookback=20]\n";
        return 1;
    }
    double gap_bps=20.0, spread_bps=1.0, impact_bps_per_pctADV=0.5;
    double gross=1.0, max_w=0.05; int lookback=20;
    bool beta_neutral = true;          // 1=on, 0=off
int  beta_lookback = 60;           // days for rolling beta regression
std::string benchmark_path = "";   // CSV path for benchmark (e.g., data/SPY.csv)
std::string benchmark_sym  = "SPY"; // label for report only

   vector<string> files;

for (int i = 1; i < argc; ++i) {
    string a = argv[i];

    if (a.rfind("--", 0) == 0) {            // it's a flag
        bool handled = false;

        // strategy / cost / sizing
        if      (a.rfind("--gap_threshold=", 0) == 0) { gap_bps = std::stod(a.substr(16)); handled = true; }
        else if (a.rfind("--spread=",         0) == 0) { spread_bps = std::stod(a.substr(9)); handled = true; }
        else if (a.rfind("--impact=",         0) == 0) { impact_bps_per_pctADV = std::stod(a.substr(9)); handled = true; }
        else if (a.rfind("--gross=",          0) == 0) { gross = std::stod(a.substr(8)); handled = true; }
        else if (a.rfind("--max_weight=",     0) == 0) { max_w = std::stod(a.substr(13)); handled = true; }
        else if (a.rfind("--lookback=",       0) == 0) { lookback = std::stoi(a.substr(11)); handled = true; }

        // beta-neutral controls
        else if (a.rfind("--beta_neutral=",   0) == 0) { beta_neutral  = (std::stoi(a.substr(15)) != 0); handled = true; }
        else if (a.rfind("--beta_lookback=",  0) == 0) { beta_lookback =  std::stoi(a.substr(16));       handled = true; }
        else if (a.rfind("--benchmark=",      0) == 0) { benchmark_path =  a.substr(12);                 handled = true; }

        if (!handled) {
            cerr << "Warning: unknown flag '" << a << "' (ignored)\n";
        }
        // (do NOT push flags into files)
    } else {
        files.push_back(a);                  // treat as a data file
    }
}

    // Load each file -> symbol -> sorted records
    unordered_map<string, vector<PriceRecord>> sym_recs;
    set<int> all_days;
    for (auto& f: files){
        auto& v = sym_recs[basename_noext(f)];
        read_csv(f, v);
        std::sort(v.begin(), v.end(), [](auto&a, auto&b){
            if(a.y!=b.y) return a.y<b.y;
            if(a.m!=b.m) return a.m<b.m;
            return a.d<b.d;
        });
        for (auto& r: v) all_days.insert(yyyymmdd(r.y,r.m,r.d));
        if (v.size()<lookback+2){
            cerr<<"Warning: "<<f<<" has too few rows.\n";
        }
    }

   vector<PriceRecord> benchmark_recs;
vector<double> benchmark_daily_returns;
unordered_map<int, PriceRecord> bench_day;   // lookup SPY bar by yyyymmdd

// Load benchmark if path provided
if (!benchmark_path.empty()) {
    read_csv(benchmark_path, benchmark_recs);
    std::sort(benchmark_recs.begin(), benchmark_recs.end(), [](const PriceRecord& a, const PriceRecord& b) {
        if (a.y != b.y) return a.y < b.y;
        if (a.m != b.m) return a.m < b.m;
        return a.d < b.d;
    });

    for (size_t i = 0; i < benchmark_recs.size(); ++i) {
    int dc = yyyymmdd(benchmark_recs[i].y, benchmark_recs[i].m, benchmark_recs[i].d);
    bench_day[dc] = benchmark_recs[i];  // store today's SPY bar

    if (i > 0) {
        double r = (benchmark_recs[i].close / benchmark_recs[i-1].close) - 1.0;
        benchmark_daily_returns.push_back(r);    // C->C series if you need it later
    }
}

if (all_days.empty()) { 
    cerr << "No data.\n"; 
    return 1; 
}

    if (all_days.empty()){ cerr<<"No data.\n"; return 1; }

    // Build per-symbol maps date->record and rolling vol (open->close returns)
    unordered_map<string, unordered_map<int, PriceRecord>> sym_day;
    unordered_map<string, unordered_map<int, double>> sym_ocret; // daily O->C return
    unordered_map<string, unordered_map<int, double>> sym_vol20; // rolling stdev

    for (auto& kv: sym_recs){
        const string& sym = kv.first;
        const auto& v = kv.second;
        vector<double> oc;
        oc.reserve(v.size());
        for (size_t i=0;i<v.size();++i){
            sym_day[sym][yyyymmdd(v[i].y,v[i].m,v[i].d)] = v[i];
            double r = (v[i].open>0)? (v[i].close - v[i].open)/v[i].open : 0.0;
            sym_ocret[sym][yyyymmdd(v[i].y,v[i].m,v[i].d)] = r;
            oc.push_back(r);
            if (i+1>= (size_t)lookback){
                double mu=0; for (size_t k=i+1-lookback;k<=i;++k) mu+=oc[k];
                mu/=lookback;
                double v2=0; for (size_t k=i+1-lookback;k<=i;++k){ double d=oc[k]-mu; v2+=d*d; }
                sym_vol20[sym][yyyymmdd(v[i].y,v[i].m,v[i].d)] = std::sqrt(v2/lookback);
            }
        }
    }

    // Iterate over calendar days, build cross-section, compute portfolio return
    double cap=1.0, start_cap=1.0;
    vector<double> daily_ret; daily_ret.reserve(all_days.size());
    vector<double> equity; equity.reserve(all_days.size());
    vector<int> month_keys, year_keys, wday_keys;
    vector<double> pnl_keys;
    vector<double> daily_ret_raw;   daily_ret_raw.reserve(all_days.size());
vector<double> bench_series;    bench_series.reserve(all_days.size());  // SPY O->C aligned by day


    double threshold = gap_bps/10000.0;
    double spread_cost = spread_bps/10000.0;

    int prev_day = -1;
    // Maintain previous close & volume per symbol to compute gaps/ADV
    unordered_map<string, double> prev_close, prev_vol;

    for (int daycode : all_days){
        int Y = daycode/10000, M=(daycode/100)%100, D=daycode%100;

        double bench_r_today = 0.0;
if (!benchmark_path.empty()) {
    auto itb = bench_day.find(daycode);
    if (itb != bench_day.end()) {
        const auto &b = itb->second;
        bench_r_today = (b.open > 0.0) ? (b.close - b.open) / b.open : 0.0;  // SPY O->C
    }
}

        // Build raw weights from signals this day
        struct Item { string sym; double ret; double adv_pct; int dir; double vol; };
        vector<Item> items;
        for (auto& kv: sym_day){
            const string& sym = kv.first;
            auto it = kv.second.find(daycode);
            if (it==kv.second.end()) continue; // no data today

            const PriceRecord& rec = it->second;
            // Need prev close & volume for gap and adv; if missing, set after
           if (!prev_close.count(sym)) { 
    prev_close[sym] = rec.close; 
    prev_vol[sym]   = rec.volume; 
    continue; 
}

            double gap = (rec.open - prev_close[sym]) / (prev_close[sym]==0?1:prev_close[sym]);
            int dir = 0;
            if (gap <= -threshold) dir = +1;       // long gap-down
            else if (gap >=  threshold) dir = -1;  // short gap-up
            if (dir==0) { prev_close[sym]=rec.close; prev_vol[sym]=rec.volume; continue; }

            // one-name open->close return
            double r = dir * (rec.close - rec.open) / (rec.open==0?1:rec.open);

            // crude ADV fraction (we trade 1 "unit" per name; here normalized)
            double adv = std::max(1.0, (double)prev_vol[sym]);
            double adv_pct = 1.0/adv;

            // rolling vol for scaling; fallback to 1%
            double vol = 0.01;
            auto volit = sym_vol20[sym].find(daycode);
            if (volit!=sym_vol20[sym].end()) vol = std::max(1e-6, volit->second);

            items.push_back({sym, r, adv_pct, dir, vol});
            prev_close[sym]=rec.close; prev_vol[sym]=rec.volume;
        }

        // Build raw weights inversely proportional to vol, cap per-name
        vector<double> w(items.size(), 0.0);
        double sumL=0, sumS=0;
        for (size_t i=0;i<items.size();++i){
            double base = 1.0 / std::max(1e-6, items[i].vol);
            double wi = items[i].dir * base;
            // cap
            if (std::fabs(wi) > max_w) wi = (wi>0? +max_w : -max_w);
            w[i] = wi;
            if (wi>0) sumL += wi; else sumS += -wi;
        }

        // If both sides present, scale to be dollar-neutral and to gross
        double gross_now = 0.0;
        if (sumL>0 && sumS>0){
            double scaleL = (0.5*gross) / sumL;
            double scaleS = (0.5*gross) / sumS;
            for (size_t i=0;i<w.size();++i){
                if (w[i]>0) w[i] *= scaleL; else w[i] *= scaleS;
                gross_now += std::fabs(w[i]);
            }
        } else {
            // No cross-section today → flat
            w.assign(w.size(), 0.0);
        }

      // Compute portfolio return with costs per name
double portfolio_r = 0.0;
for (size_t i=0; i<items.size(); ++i) {
    if (w[i]==0.0) continue;
    double impact = (impact_bps_per_pctADV/10000.0) * items[i].adv_pct;
    double total_cost = 2.0*(spread_cost + impact); // enter+exit
    double ri = items[i].ret - total_cost;
    portfolio_r += w[i] * ri;
}
double portfolio_r_raw = portfolio_r;
// ===== Beta-neutral adjustment (step 4) =====
if (beta_neutral && !benchmark_daily_returns.empty() && bench_day.count(daycode)) {
    double bench_r = (bench_day[daycode].close / bench_day[daycode].open) - 1.0;

    // Compute covariance with benchmark (rolling beta)
    double cov = 0.0, var_b = 0.0;
    size_t lookback_n = std::min(daily_ret.size(), (size_t)beta_lookback);
    if (lookback_n > 1) {
        double mu_p = 0.0, mu_b = 0.0;
        for (size_t k=0; k<lookback_n; ++k) {
            mu_p += daily_ret[daily_ret.size()-1-k];
            mu_b += benchmark_daily_returns[benchmark_daily_returns.size()-1-k];
        }
        mu_p /= lookback_n;
        mu_b /= lookback_n;
        for (size_t k=0; k<lookback_n; ++k) {
            double dp = daily_ret[daily_ret.size()-1-k] - mu_p;
            double db = benchmark_daily_returns[benchmark_daily_returns.size()-1-k] - mu_b;
            cov   += dp * db;
            var_b += db * db;
        }
        if (var_b > 0.0) {
            double beta = cov / var_b;
            portfolio_r -= beta * bench_r;   // adjust return
        }
    }
}
// ===== End beta-neutral adjustment =====

cap *= (1.0 + portfolio_r);
daily_ret_raw.push_back(portfolio_r_raw);   // unhedged
daily_ret.push_back(portfolio_r);           // hedged
bench_series.push_back(bench_r_today);      // benchmark
pnl_keys.push_back(portfolio_r);
equity.push_back(cap);
month_keys.push_back(M);
year_keys.push_back(Y);
wday_keys.push_back(weekday(Y,M,D));

    // Summary stats
    double tot = (cap - start_cap)/start_cap * 100.0;
    double mu=0; for (double r: daily_ret) mu+=r; mu/=std::max<size_t>(1,daily_ret.size());
    double var=0; for (double r: daily_ret){ double d=r-mu; var+=d*d; }
    var/=std::max<size_t>(1,daily_ret.size());
    double sd = std::sqrt(var);
    int wins=0; for (double r: daily_ret) if (r>0) ++wins;
    auto [dd, sIdx, eIdx] = max_dd(equity);
    double v95 = var95(daily_ret);
    double cv95 = cvar95(daily_ret);

    // Aggregations
    auto agg = [&](const vector<int>& keys){
        unordered_map<int,double> m; for (size_t i=0;i<keys.size();++i) m[keys[i]]+=pnl_keys[i]; return m;
    };
    auto monthly = agg(month_keys);
    auto yearly  = agg(year_keys);
    auto wday    = agg(wday_keys);

    cout<<"Market-Neutral Gap-Fade Portfolio Report\n";
    cout<<"=======================================\n";
    cout<<"Files loaded: "; for (auto& f: files) cout<<f<<" "; cout<<"\n";
    cout<<"Trading days: "<<daily_ret.size()<<"\n";
    cout<<"Gross leverage: "<<gross<<" | Per-name cap: "<<max_w<<" | Lookback: "<<lookback<<"\n";
    cout<<"Threshold: "<<gap_bps<<" bps | Spread: "<<spread_bps<<" bps | Impact: "<<impact_bps_per_pctADV<<" bps/%ADV\n\n";

    cout<<"Total return: "<<tot<<" %\n";
    cout<<"Mean daily: "<<mu*100<<" % | Stdev: "<<sd*100<<" %\n";
    cout<<"Sharpe: "<<sharpe(daily_ret)<<" | Sortino: "<<sortino(daily_ret)<<"\n";
    cout<<"Max drawdown: "<<dd*100<<" %\n";
    cout<<"Win rate: "<<(100.0*wins/std::max<size_t>(1,daily_ret.size()))<<" %\n";
    cout<<"VaR 95%: "<<v95*100<<" % | CVaR 95%: "<<cv95*100<<" %\n\n";

    cout<<"Monthly PnL (% of capital):\n";
    for (auto& kv: monthly) cout<<"Month "<<kv.first<<": "<<kv.second*100<<" %\n";
    cout<<"\nYearly PnL (% of capital):\n";
    for (auto& kv: yearly) cout<<"Year "<<kv.first<<": "<<kv.second*100<<" %\n";
    static const char* wdN[]={"Mon","Tue","Wed","Thu","Fri","Sat","Sun"};
    cout<<"\nWeekday PnL (% of capital):\n";
    for (auto& kv: wday) if (kv.first>=0 && kv.first<7) cout<<wdN[kv.first]<<": "<<kv.second*100<<" %\n";
    return 0;
}
