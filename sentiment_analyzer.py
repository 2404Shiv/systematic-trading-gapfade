"""
Free Sentiment Analysis Module
Uses NewsAPI (free tier) + keyword scoring — no AI API costs
Works on Railway free tier
"""
import os
import logging
import requests
from datetime import datetime, timedelta

log = logging.getLogger(__name__)

NEWSAPI_KEY = os.environ.get("NEWSAPI_KEY", "")

# Financial sentiment keywords
BULLISH = [
    "beat", "beats", "surge", "surges", "record", "growth", "profit",
    "upgrade", "strong", "soars", "rally", "raises", "bullish", "outperform",
    "exceeds", "positive", "gain", "gains", "up", "rise", "rises", "high",
    "buy", "opportunity", "recover", "recovery", "boost", "boosts", "expand",
    "expansion", "revenue", "earnings beat", "better than expected", "raised guidance"
]

BEARISH = [
    "miss", "misses", "fall", "falls", "loss", "losses", "downgrade", "weak",
    "cut", "cuts", "crash", "crashes", "layoff", "layoffs", "below", "decline",
    "declines", "down", "drop", "drops", "low", "sell", "risk", "warning",
    "disappoint", "disappoints", "disappointing", "worse than expected",
    "lowered guidance", "investigation", "lawsuit", "recall", "bankruptcy",
    "debt", "concern", "concerns", "selloff", "bearish", "underperform"
]

COMPANY_NAMES = {
    "AAPL": "Apple", "MSFT": "Microsoft", "GOOGL": "Google Alphabet",
    "AMZN": "Amazon", "META": "Meta Facebook", "TSLA": "Tesla",
    "NVDA": "Nvidia", "JPM": "JPMorgan Chase", "V": "Visa",
    "MA": "Mastercard", "UNH": "UnitedHealth", "JNJ": "Johnson Johnson",
    "PG": "Procter Gamble", "HD": "Home Depot", "BAC": "Bank of America",
    "XOM": "Exxon Mobil", "CVX": "Chevron", "ABBV": "AbbVie",
    "MRK": "Merck", "PEP": "PepsiCo", "KO": "Coca Cola",
    "AVGO": "Broadcom", "COST": "Costco", "MCD": "McDonalds",
    "QCOM": "Qualcomm", "GS": "Goldman Sachs", "BLK": "BlackRock",
    "CAT": "Caterpillar", "AXP": "American Express", "LOW": "Lowes",
    "SBUX": "Starbucks", "INTU": "Intuit", "DE": "John Deere",
}

def get_headlines(symbol: str) -> list:
    """Fetch recent headlines for a stock"""
    if not NEWSAPI_KEY:
        return []
    
    company = COMPANY_NAMES.get(symbol, symbol)
    yesterday = (datetime.utcnow() - timedelta(days=1)).strftime("%Y-%m-%d")
    
    try:
        res = requests.get(
            "https://newsapi.org/v2/everything",
            params={
                "q": f"{company} stock",
                "from": yesterday,
                "sortBy": "publishedAt",
                "language": "en",
                "pageSize": 10,
                "apiKey": NEWSAPI_KEY
            },
            timeout=5
        )
        data = res.json()
        if data.get("status") == "ok":
            return [a["title"].lower() for a in data.get("articles", [])]
    except Exception as e:
        log.error(f"News error {symbol}: {e}")
    return []

def score_headlines(headlines: list) -> float:
    """
    Score headlines from -1.0 (very bearish) to +1.0 (very bullish)
    """
    if not headlines:
        return 0.0
    
    bullish_count = 0
    bearish_count = 0
    
    for headline in headlines:
        for word in BULLISH:
            if word in headline:
                bullish_count += 1
        for word in BEARISH:
            if word in headline:
                bearish_count += 1
    
    total = bullish_count + bearish_count
    if total == 0:
        return 0.0
    
    # Score between -1 and 1
    score = (bullish_count - bearish_count) / total
    return round(score, 3)

def should_trade(symbol: str, gap_pct: float, direction: int) -> dict:
    """
    Main function — decides whether to trade based on sentiment
    
    direction: +1 = long (fading gap down), -1 = short (fading gap up)
    
    Logic:
    - Gap DOWN + BULLISH news → LONG (fade + sentiment agree) ✅
    - Gap DOWN + BEARISH news → SHORT (bet on continuation) ✅  
    - Gap UP + BEARISH news → SHORT (fade + sentiment agree) ✅
    - Gap UP + BULLISH news → LONG (bet on continuation) ✅
    - Neutral news → follow gap fade signal ✅
    - Conflicting weak signals → follow gap fade signal ✅
    """
    headlines = get_headlines(symbol)
    score = score_headlines(headlines)
    
    STRONG_THRESHOLD = 0.3  # Score above this = strong signal
    
    # Determine action
    gap_is_down = gap_pct < 0
    gap_is_up = gap_pct > 0
    strongly_bullish = score > STRONG_THRESHOLD
    strongly_bearish = score < -STRONG_THRESHOLD
    
    if gap_is_down and strongly_bullish:
        # Gap down + good news = strong fade signal — LONG
        action = "LONG"
        trade = True
        reason = f"Gap down + bullish sentiment ({score:.2f}) — strong fade signal"
    
    elif gap_is_down and strongly_bearish:
        # Gap down + bad news = continuation — SHORT instead
        action = "SHORT"
        trade = True
        reason = f"Gap down + bearish sentiment ({score:.2f}) — betting on continuation"
    
    elif gap_is_up and strongly_bearish:
        # Gap up + bad news = strong fade signal — SHORT
        action = "SHORT"
        trade = True
        reason = f"Gap up + bearish sentiment ({score:.2f}) — strong fade signal"
    
    elif gap_is_up and strongly_bullish:
        # Gap up + good news = continuation — LONG instead
        action = "LONG"
        trade = True
        reason = f"Gap up + bullish sentiment ({score:.2f}) — betting on continuation"
    
    else:
        # Neutral or weak sentiment — follow gap fade
        action = "LONG" if direction > 0 else "SHORT"
        trade = True
        reason = f"Neutral sentiment ({score:.2f}) — following gap fade signal"
    
    result = {
        "score": score,
        "action": action,
        "trade": trade,
        "reasoning": reason,
        "headlines_count": len(headlines),
        "direction": 1 if action == "LONG" else -1
    }
    
    log.info(f"{symbol}: score={score:.2f} | {action} | {reason}")
    return result

def analyze_batch(signals: list) -> list:
    """Filter and enhance signals with sentiment"""
    if not NEWSAPI_KEY:
        log.warning("No NEWSAPI_KEY set — skipping sentiment, trading all signals")
        return signals
    
    enhanced = []
    for s in signals:
        sym = s["sym"]
        sentiment = should_trade(sym, s["gap"] * 100, s["direction"])
        
        if not sentiment["trade"]:
            log.info(f"SKIP {sym}: {sentiment['reasoning']}")
            continue
        
        # Override direction if sentiment says opposite
        if sentiment["direction"] != s["direction"]:
            log.info(f"DIRECTION CHANGE {sym}: {s['direction']} → {sentiment['direction']}")
            s["direction"] = sentiment["direction"]
        
        s["sentiment_score"] = sentiment["score"]
        s["sentiment_action"] = sentiment["action"]
        s["sentiment_reason"] = sentiment["reasoning"]
        
        enhanced.append(s)
    
    log.info(f"Sentiment: {len(signals)} signals → {len(enhanced)} trades")
    return enhanced
