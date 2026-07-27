#!/usr/bin/env python3
"""qdata — load price history from public sources or your own files. STDLIB ONLY (urllib).

The ONLY module here that touches the network. Everything else in this skill is offline, so if
your environment forbids outbound calls, use `from_csv` and the rest still works.

WHAT YOU CAN ACTUALLY GET (measured, not assumed — see `availability()`):

  daily      ~20 years for index futures, decades for indices          → statistically usable
  hourly     ~2 years                                                   → usable
  5-minute   ~60 calendar days (~40-50 sessions) — HARD LIMIT           → NOT enough for most claims
  1-minute   ~7 calendar days (~5 sessions)                             → diagnostic only

The intraday limit cannot be worked around by requesting older windows: the provider returns empty
for any intraday range beyond its retention. This is the binding constraint on evaluating any
intraday strategy (opening-range breakout, session fades, scalps) from free data.

Consequences, stated plainly:
  - A DAILY strategy can be evaluated on a sample long enough to matter.
  - An INTRADAY strategy cannot, from free sources. Either supply your own history (broker or
    vendor export via `from_csv`), or run a forward collector and wait (`append_csv`).
  - `qstats.required_n()` will usually tell you the free intraday sample is an order of magnitude
    too small. Print that number instead of running the backtest anyway.
"""
import csv
import datetime as _dt
import io
import json
import os
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"}
YAHOO = "https://query1.finance.yahoo.com/v8/finance/chart/"

# measured 2026-07-26; re-check with availability() rather than trusting these
RETENTION = {"1m": "~7 calendar days", "2m": "~60 days", "5m": "~60 days",
             "15m": "~60 days", "30m": "~60 days", "1h": "~2 years",
             "1d": "~20 years (use range=20y; range=max silently returns MONTHLY bars)",
             "1wk": "decades", "1mo": "decades"}


def _get(url, timeout=25):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()


def fetch(symbol, interval="1d", range_="20y", session_only=False, tz_offset_hours=0):
    """Fetch OHLCV bars. Returns a list of dicts: t (epoch s), dt (UTC datetime), o,h,l,c,v.

    symbol examples: 'NQ=F' (Nasdaq-100 futures), '^NDX', '^IXIC', 'QQQ', 'ES=F', 'BTC-USD'.
    NOTE: interval='1d' with range_='max' returns MONTHLY bars — use range_='20y' for real daily.
    """
    url = f"{YAHOO}{symbol}?range={range_}&interval={interval}"
    j = json.loads(_get(url))
    res = j["chart"]["result"][0]
    ts = res.get("timestamp") or []
    q = res["indicators"]["quote"][0]
    rows = []
    for i, t in enumerate(ts):
        o, h, l, c = q["open"][i], q["high"][i], q["low"][i], q["close"][i]
        if None in (o, h, l, c):
            continue
        d = _dt.datetime.fromtimestamp(t, _dt.timezone.utc) + _dt.timedelta(hours=tz_offset_hours)
        rows.append({"t": t, "dt": d, "o": float(o), "h": float(h), "l": float(l),
                     "c": float(c), "v": float(q.get("volume", [0] * len(ts))[i] or 0)})
    if session_only:
        rows = [r for r in rows if r["dt"].weekday() < 5]
    return rows


def availability(symbol="NQ=F", intervals=("1d", "1h", "5m", "1m")):
    """Measure what this provider actually returns for a symbol, right now.
    Run this BEFORE promising an analysis — the answer decides whether the question is answerable."""
    out = []
    for iv in intervals:
        rng = {"1d": "20y", "1h": "2y", "5m": "60d", "1m": "7d"}.get(iv, "60d")
        try:
            rows = fetch(symbol, iv, rng)
            days = len({r["dt"].date() for r in rows})
            wk = len({r["dt"].date() for r in rows if r["dt"].weekday() < 5})
            out.append({"interval": iv, "range": rng, "bars": len(rows),
                        "calendar_days": days, "weekday_sessions": wk,
                        "first": str(rows[0]["dt"])[:16] if rows else None,
                        "last": str(rows[-1]["dt"])[:16] if rows else None,
                        "documented_retention": RETENTION.get(iv, "?")})
        except Exception as e:
            out.append({"interval": iv, "range": rng, "error": f"{type(e).__name__}: {e}"[:70]})
    return out


def to_returns(rows, field="c", log=False):
    """Bar-to-bar returns from a fetched series."""
    out = []
    for a, b in zip(rows, rows[1:]):
        if a[field] and a[field] > 0:
            out.append(math.log(b[field] / a[field]) if log else (b[field] / a[field] - 1.0))
    return out


def daily_returns_by_session(rows, field="c"):
    """Collapse intraday bars into one return per session (close-to-close)."""
    by = {}
    for r in rows:
        by.setdefault(r["dt"].date(), []).append(r)
    days = sorted(by)
    closes = [by[d][-1][field] for d in days]
    return days[1:], [closes[i + 1] / closes[i] - 1.0 for i in range(len(closes) - 1)]


def group_sessions(rows, start=(9, 30), end=(16, 0), tz_offset_hours=-4):
    """Group bars into trading sessions filtered to a wall-clock window.
    tz_offset_hours converts UTC to the exchange clock (-4 = US Eastern in summer, -5 in winter).
    Returns {date: [bars]} with weekends removed."""
    out = {}
    for r in rows:
        d = r["dt"] + _dt.timedelta(hours=tz_offset_hours)
        if d.weekday() >= 5:
            continue
        if start <= (d.hour, d.minute) < end:
            out.setdefault(d.date(), []).append({**r, "local": d})
    return {k: v for k, v in sorted(out.items()) if v}

# ---------------------------------------------------------------- your own files


def from_databento_csv(path, price_scale=1e-9):
    """Load a Databento OHLCV CSV. Two conversions that WILL bite if forgotten:
      - prices are fixed-point integers where 1 unit = 1e-9, so divide by 1e9;
      - ts_event is nanoseconds since epoch, not seconds.
    Verified against the catalogue: NQ ends at 28,306.50 and GC at 4,055.70, matching the
    live quotes."""
    rows = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        for r in rd:
            try:
                ns = int(r["ts_event"])
                ts = ns // 1_000_000_000
                rows.append({"t": ts,
                             "dt": _dt.datetime.fromtimestamp(ts, _dt.timezone.utc),
                             "o": int(r["open"]) * price_scale,
                             "h": int(r["high"]) * price_scale,
                             "l": int(r["low"]) * price_scale,
                             "c": int(r["close"]) * price_scale,
                             "v": float(r.get("volume") or 0),
                             "instrument_id": r.get("instrument_id")})
            except (KeyError, ValueError, TypeError):
                continue
    rows.sort(key=lambda r: r["t"])
    return rows


def from_csv(path, t_col=0, o_col=1, h_col=2, l_col=3, c_col=4, v_col=None,
             has_header=True, dt_format=None):
    """Load YOUR history — a broker export, a vendor file, a TradingView download.
    This is the path that removes every retention limit above, and the one to prefer for any
    serious intraday claim."""
    rows = []
    with open(path, newline="") as f:
        rd = csv.reader(f)
        if has_header:
            next(rd, None)
        for r in rd:
            if not r or len(r) <= max(x for x in (t_col, o_col, h_col, l_col, c_col) if x is not None):
                continue
            raw = r[t_col].strip()
            try:
                d = (_dt.datetime.strptime(raw, dt_format) if dt_format
                     else (_dt.datetime.fromtimestamp(float(raw), _dt.timezone.utc)
                           if raw.replace(".", "").isdigit()
                           else _dt.datetime.fromisoformat(raw.replace("Z", "+00:00"))))
            except Exception:
                continue
            try:
                rows.append({"t": int(d.timestamp()), "dt": d, "o": float(r[o_col]),
                             "h": float(r[h_col]), "l": float(r[l_col]), "c": float(r[c_col]),
                             "v": float(r[v_col]) if v_col is not None and r[v_col] else 0.0})
            except (ValueError, IndexError):
                continue
    rows.sort(key=lambda r: r["t"])
    return rows


def append_csv(path, rows, key="t"):
    """Idempotent forward collector: merge new bars into a store, deduped by timestamp.
    The only way to build a long intraday sample when no vendor will sell you the history."""
    have = {}
    if os.path.exists(path):
        for r in from_csv(path):
            have[r[key]] = r
    for r in rows:
        have[r[key]] = r
    ks = sorted(have)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "open", "high", "low", "close", "volume"])
        for k in ks:
            r = have[k]
            w.writerow([r["t"], r["o"], r["h"], r["l"], r["c"], r["v"]])
    return {"path": path, "total_bars": len(ks), "added": len(rows),
            "first": str(_dt.datetime.fromtimestamp(ks[0], _dt.timezone.utc))[:16] if ks else None,
            "last": str(_dt.datetime.fromtimestamp(ks[-1], _dt.timezone.utc))[:16] if ks else None}


import math  # noqa: E402  (used by to_returns)


def report(symbol="NQ=F"):
    L = [f"=== data availability · {symbol} (measured now) ===",
         f"  {'interval':<9}{'bars':>8}{'cal.days':>10}{'sessions':>10}  window"]
    for a in availability(symbol):
        if "error" in a:
            L.append(f"  {a['interval']:<9}{'ERR':>8}  {a['error']}")
            continue
        L.append(f"  {a['interval']:<9}{a['bars']:>8}{a['calendar_days']:>10}"
                 f"{a['weekday_sessions']:>10}  {a['first']} → {a['last']}")
    L.append("  intraday retention is a hard provider limit; older intraday windows return empty.")
    L.append("  For an intraday strategy: supply your own history (from_csv) or collect forward "
             "(append_csv), and check qstats.required_n() before drawing conclusions.")
    return "\n".join(L)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--symbol", default="NQ=F")
    a = ap.parse_args()
    print(report(a.symbol))


# ═══════════════════════════════════════════════════════════════════════════════
# Databento (paid, historical CME/ICE/equities/options) and Binance (free crypto)
# ═══════════════════════════════════════════════════════════════════════════════
# The API key is read from the environment and NEVER passed as an argument, printed,
# or written to disk by this module. Set it once in your shell:
#     export DATABENTO_API_KEY='...'
# Costs are ESTIMATED BEFORE downloading; nothing is fetched until you pass confirm=True.

DB_HIST = "https://hist.databento.com/v0/"


DB_KEY_FILES = (os.path.expanduser("~/.databento.env"), os.path.expanduser("~/.env"))


def _db_key():
    """Read the key from the environment, else from a protected dotfile. Never from an argument,
    so it cannot end up in a shell history, a log, a traceback or a chat transcript."""
    key = os.environ.get("DATABENTO_API_KEY", "").strip()
    if key:
        return key
    for path in DB_KEY_FILES:
        if not os.path.exists(path):
            continue
        mode = os.stat(path).st_mode & 0o777
        if mode & 0o077:
            raise RuntimeError(f"{path} is readable by others (mode {mode:o}). "
                               f"Run: chmod 600 {path}")
        for line in open(path):
            line = line.strip()
            if line.startswith("DATABENTO_API_KEY"):
                return line.split("=", 1)[1].strip().strip("'\"")
    return ""


def _db_auth_header():
    key = _db_key()
    if not key:
        raise RuntimeError("No Databento key found. Put it in ~/.databento.env (chmod 600) as "
                           "DATABENTO_API_KEY=... or export it. Never pass it as an argument.")
    import base64
    token = base64.b64encode(f"{key}:".encode()).decode()
    return {"Authorization": f"Basic {token}", "User-Agent": "quant-rigor/1.3"}


def _db_get(endpoint, params, timeout=180):
    import urllib.parse
    url = DB_HIST + endpoint + "?" + urllib.parse.urlencode(params, doseq=True)
    req = urllib.request.Request(url, headers=_db_auth_header())
    return urllib.request.urlopen(req, timeout=timeout).read()


def db_cost(symbols, start, end, dataset="GLBX.MDP3", schema="ohlcv-1m",
            stype_in="continuous"):
    """Estimate the cost and size of a download BEFORE spending credit.
    symbols: e.g. ['NQ.c.0'] (continuous front month) or ['NQZ6'] with stype_in='raw_symbol'."""
    p = {"dataset": dataset, "symbols": ",".join(symbols), "schema": schema,
         "start": start, "end": end, "stype_in": stype_in, "mode": "historical"}
    cost = json.loads(_db_get("metadata.get_cost", p))
    size = json.loads(_db_get("metadata.get_record_count", dict(p)))
    try:
        bs = json.loads(_db_get("metadata.get_billable_size", dict(p)))
    except Exception:
        bs = None
    return {"symbols": symbols, "schema": schema, "start": start, "end": end,
            "usd_cost": cost, "record_count": size, "billable_bytes": bs,
            "billable_gb": (bs / 1e9) if isinstance(bs, (int, float)) else None}


def db_download(symbols, start, end, out_path, dataset="GLBX.MDP3", schema="ohlcv-1m",
                stype_in="continuous", encoding="csv", confirm=False):
    """Download to `out_path`. Refuses to run unless confirm=True — check db_cost() first."""
    est = db_cost(symbols, start, end, dataset, schema, stype_in)
    if not confirm:
        return {"downloaded": False, "estimate": est,
                "note": "pass confirm=True to actually download (this spends credit)"}
    p = {"dataset": dataset, "symbols": ",".join(symbols), "schema": schema,
         "start": start, "end": end, "stype_in": stype_in, "encoding": encoding}
    blob = _db_get("timeseries.get_range", p, timeout=1800)
    with open(out_path, "wb") as f:
        f.write(blob)
    return {"downloaded": True, "path": out_path, "bytes": len(blob), "estimate": est}


# ---------------------------------------------------------------- free crypto (Binance)

BINANCE_DUMP = "https://data.binance.vision/data/{market}/monthly/klines/{sym}/{iv}/{sym}-{iv}-{ym}.zip"


def binance_months(symbol, interval="1m", start_ym="2020-01", end_ym=None, market="spot"):
    """List the monthly dump URLs for a symbol. Free, no API key, no rate limit worth worrying
    about. Klines exist from 2017-08 for BTCUSDT. market: 'spot' | 'futures/um'."""
    y, m = (int(x) for x in start_ym.split("-"))
    if end_ym is None:
        t = _dt.datetime.now(_dt.timezone.utc)
        end_ym = f"{t.year}-{t.month:02d}"
    ey, em = (int(x) for x in end_ym.split("-"))
    out = []
    while (y, m) <= (ey, em):
        out.append({"ym": f"{y}-{m:02d}",
                    "url": BINANCE_DUMP.format(market=market, sym=symbol.upper(),
                                               iv=interval, ym=f"{y}-{m:02d}")})
        m += 1
        if m > 12:
            y, m = y + 1, 1
    return out


def binance_fetch(symbol, interval="1m", start_ym="2020-01", end_ym=None, market="spot",
                  out_csv=None, verbose=True):
    """Download and concatenate Binance monthly klines into OHLCV rows (and optionally a CSV).
    Kline columns: open_time, open, high, low, close, volume, close_time, ..."""
    import zipfile
    rows, missing = [], []
    for item in binance_months(symbol, interval, start_ym, end_ym, market):
        try:
            blob = _get(item["url"], timeout=120)
        except Exception:
            missing.append(item["ym"])
            continue
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            name = z.namelist()[0]
            for line in z.read(name).decode().splitlines():
                p = line.split(",")
                if len(p) < 6 or not p[0].strip().replace(".", "").isdigit():
                    continue
                raw = int(float(p[0]))
                # Binance switched kline timestamps from milliseconds to MICROSECONDS in 2025.
                # Normalise all three possible scales to seconds.
                ts = (raw // 1_000_000 if raw > 1e14 else
                      raw // 1_000 if raw > 1e11 else raw)
                rows.append({"t": ts,
                             "dt": _dt.datetime.fromtimestamp(ts, _dt.timezone.utc),
                             "o": float(p[1]), "h": float(p[2]), "l": float(p[3]),
                             "c": float(p[4]), "v": float(p[5])})
        if verbose:
            print(f"  {symbol} {item['ym']}: {len(rows)} righe cumulate")
    rows.sort(key=lambda r: r["t"])
    if out_csv:
        append_csv(out_csv, rows)
    return {"symbol": symbol, "interval": interval, "rows": len(rows),
            "first": str(rows[0]["dt"])[:16] if rows else None,
            "last": str(rows[-1]["dt"])[:16] if rows else None,
            "missing_months": missing, "csv": out_csv}


# ---------------------------------------------------------------- CFD -> underlying map

CFD_TO_FUTURES = {
    # what a CFD broker calls it        ->  the exchange-traded instrument that actually has history
    "XAUUSD (gold)":        ("GC / MGC (micro)", "GLBX.MDP3", "COMEX gold futures"),
    "XAGUSD (silver)":      ("SI / SIL",         "GLBX.MDP3", "COMEX silver futures"),
    "NAS100 / US100":       ("NQ / MNQ",         "GLBX.MDP3", "E-mini Nasdaq-100"),
    "US500 / SPX500":       ("ES / MES",         "GLBX.MDP3", "E-mini S&P 500"),
    "US30 / DJ30":          ("YM / MYM",         "GLBX.MDP3", "E-mini Dow"),
    "US2000":               ("RTY / M2K",        "GLBX.MDP3", "E-mini Russell 2000"),
    "GER40 / DAX":          ("FDAX / FDXM",      "XEUR",      "Eurex DAX futures (history from 2025)"),
    "UK100 / FTSE":         ("Z",                "IFLL",      "ICE Europe Financials FTSE"),
    "JP225 / Nikkei":       ("NKD / NIY",        "GLBX.MDP3", "CME Nikkei futures"),
    "USOIL / WTI":          ("CL / MCL",         "GLBX.MDP3", "NYMEX crude"),
    "UKOIL / Brent":        ("B",                "IFEU",      "ICE Brent"),
    "NATGAS":               ("NG / QG",          "GLBX.MDP3", "NYMEX natural gas"),
    "EURUSD":               ("6E / M6E",         "GLBX.MDP3", "CME euro FX futures"),
    "GBPUSD":               ("6B / M6B",         "GLBX.MDP3", "CME sterling futures"),
    "USDJPY":               ("6J",               "GLBX.MDP3", "CME yen futures (inverted quote)"),
    "AUDUSD":               ("6A / M6A",         "GLBX.MDP3", "CME Aussie futures"),
    "DXY (dollar index)":   ("DX",               "IFUS",      "ICE US Dollar Index"),
    "VIX":                  ("VX",               "XCBF",      "Cboe VIX futures"),
    "BTCUSD":               ("BTC / MBT",        "GLBX.MDP3", "CME bitcoin futures — or free Binance spot"),
    "ETHUSD":               ("ETH / MET",        "GLBX.MDP3", "CME ether futures — or free Binance spot"),
}


def cfd_note():
    """CFDs have no centralised history: they are OTC contracts, each broker quoting its own book.
    There is no vendor that sells 'CFD data' for this reason. Test on the underlying listed
    instrument and add YOUR broker's spread and financing as explicit costs — that is both the
    honest and the only reproducible approach."""
    L = ["CFDs are OTC: no vendor has a canonical price history, because each broker has its own.",
         "Use the listed underlying, then add your broker's spread/financing as measured costs.",
         "", f"{'CFD name':<24}{'listed instrument':<20}{'dataset':<12}note"]
    for k, (sym, ds, note) in CFD_TO_FUTURES.items():
        L.append(f"{k:<24}{sym:<20}{ds:<12}{note}")
    return "\n".join(L)
