#!/usr/bin/env python3
"""qmicro — microstructure & transaction-cost estimation. STDLIB ONLY.

Answers the question this package could not answer before: *what does it actually cost*,
estimated from the price data you already have, rather than assumed.

Added after studying (see references/08-microstructure-costs.md for provenance):
  - stefan-jansen/machine-learning-for-trading, ch.18 (cost taxonomy, spread estimation,
    cost cliff, commission-vs-slippage decomposition)
  - the Roll (1984) and Corwin-Schultz (2012) estimators as implemented there

WHY THIS MATTERS MORE THAN IT LOOKS
An assumed cost model silently decides verdicts. A strategy with a 40pt target and a 7pt
round-trip is a different economic proposition from the same strategy at 0.7pt, and no amount
of statistical machinery downstream can recover from getting this wrong. Estimate it, then
stress it — never assume it once and forget.
"""
import math
import random

# ---------------------------------------------------------------- 50. Roll (1984)


def roll_spread(closes, min_relative=1e-9):
    """Roll's serial-covariance spread estimator.

    With trades alternating between bid and ask, consecutive price changes have
    covariance -(s/2)^2, so s = 2*sqrt(-cov). Returns the spread in PRICE UNITS.

    Returns None when the covariance is positive: that means prices are trending
    (positive autocorrelation) and the estimator is undefined. Reporting a spread of 0
    there would be a silent lie — the caller must handle None.

    FAILURE MODE: on trending data the covariance flips positive and Roll breaks. It also
    conflates genuine mean reversion with bounce, which is exactly why a measured
    short-horizon reversion that matches the tick size is NOT a tradable edge.
    """
    n = len(closes)
    if n < 3:
        return None
    d = [closes[i] - closes[i - 1] for i in range(1, n)]
    m = sum(d) / len(d)
    var_d = sum((x - m) ** 2 for x in d) / (len(d) - 1)
    if var_d <= 0:
        return None                      # serie deterministica: nessuna informazione
    cov = sum((d[i] - m) * (d[i - 1] - m) for i in range(1, len(d))) / (len(d) - 1)
    if cov >= 0:
        return None
    s = 2.0 * math.sqrt(-cov)
    # Guardia FISICA, non statistica. L'arrotondamento in virgola mobile su una serie
    # deterministica produce differenze che si alternano, quindi covarianza fortemente
    # negativa: il test sul segno da solo non basta. Uno spread di 1e-15 su un prezzo di 100
    # non e' uno spread. Sotto min_relative del prezzo lo dichiaro non misurato.
    mid = sorted(closes)[len(closes) // 2]
    if mid and s / abs(mid) < min_relative:
        return None
    return s


def roll_spread_relative(closes):
    """Roll spread as a fraction of price (comparable across instruments and eras)."""
    s = roll_spread(closes)
    if s is None:
        return None
    mid = sorted(closes)[len(closes) // 2]
    return s / mid if mid else None


# ---------------------------------------------------------------- 51. Corwin-Schultz (2012)

_K = 3.0 - 2.0 * math.sqrt(2.0)


def corwin_schultz_spread(highs, lows, window=1):
    """Corwin-Schultz (2012) high-low spread estimator. Returns spread as a FRACTION.

    Insight: over a short interval the high is predominantly transacted at the ask and the
    low at the bid, so the high-low range contains the spread plus volatility. Comparing a
    single-period range with a two-period range separates the two, because volatility scales
    with time while the spread does not.

        beta  = ln(H_t/L_t)^2 + ln(H_{t-1}/L_{t-1})^2
        gamma = ln(max(H_t,H_{t-1}) / min(L_t,L_{t-1}))^2
        alpha = (sqrt(2*beta) - sqrt(beta)) / (3 - 2*sqrt(2)) - sqrt(gamma / (3 - 2*sqrt(2)))
        S     = 2*(e^alpha - 1) / (1 + e^alpha)

    Negative estimates are clamped to 0 (they occur when volatility swamps the spread).
    `window` averages beta and gamma over that many observations before solving, which is
    what the original paper recommends to reduce noise.

    WHY PREFER IT TO ROLL: it uses the high-low range rather than the serial covariance of
    closes, so it does not break on trending data and does not confuse mean reversion with
    bounce. On OHLC bars it is usually the better of the two. Report BOTH — when they
    disagree materially, that disagreement is information about the data, not a nuisance.
    """
    n = min(len(highs), len(lows))
    if n < 2:
        return []
    betas, gammas = [], []
    for t in range(1, n):
        if min(lows[t], lows[t - 1]) <= 0:
            betas.append(None); gammas.append(None); continue
        b = math.log(highs[t] / lows[t]) ** 2 + math.log(highs[t - 1] / lows[t - 1]) ** 2
        h2 = max(highs[t], highs[t - 1])
        l2 = min(lows[t], lows[t - 1])
        g = math.log(h2 / l2) ** 2
        betas.append(b); gammas.append(g)

    def _avg(seq, i, w):
        lo = max(0, i - w + 1)
        vals = [v for v in seq[lo:i + 1] if v is not None]
        return (sum(vals) / len(vals)) if vals else None

    out = []
    for i in range(len(betas)):
        b = _avg(betas, i, window) if window > 1 else betas[i]
        g = _avg(gammas, i, window) if window > 1 else gammas[i]
        if b is None or g is None or b < 0 or g < 0:
            out.append(None); continue
        alpha = (math.sqrt(2.0 * b) - math.sqrt(b)) / _K - math.sqrt(g / _K)
        s = 2.0 * (math.exp(alpha) - 1.0) / (1.0 + math.exp(alpha))
        out.append(max(0.0, s))
    return out


def corwin_schultz_median(highs, lows, window=1):
    """Median Corwin-Schultz estimate — the single number to quote."""
    vals = [v for v in corwin_schultz_spread(highs, lows, window) if v is not None]
    if not vals:
        return None
    vals.sort()
    return vals[len(vals) // 2]


# ---------------------------------------------------------------- 52. is the reversion real?


def plausibility_flag(spread_rel, mid_price, tick_size, max_ticks=20.0,
                     median_bar_range_rel=None, max_range_share=0.5):
    """E' credibile che quella stima sia davvero lo spread? DUE guardie, non una.

    SCOPERTA D'USO 1 (2026-07): Roll e Corwin-Schultz sono stimatori di MICROSTRUTTURA. Su
    barre GIORNALIERE di futures liquidi misurano volatilita', non spread: su ES daily, Roll
    ha dato 62bp contro uno spread reale di ~0.5bp, cioe' ~53 tick. Guardia A: sopra ~20
    tick impliciti la stima non e' uno spread.

    SCOPERTA D'USO 2 (2026-07-30, e la prima guardia da sola NON l'ha presa): la guardia a
    tick impliciti non discrimina quando il TICK E' MOLTO PIU' FINE DELLO SPREAD. Su BTC il
    tick e' $0.01 su un prezzo di $58 000: Roll ha dato 17.87bp = 4 986 tick impliciti, un
    numero enorme ma privo di significato, perche' su quell'asset anche lo spread VERO vale
    ~580 tick. Il conteggio in tick non e' una scala universale.

    Guardia B, FISICA e indipendente dal tick: uno spread non puo' superare di routine il
    RANGE della barra che lo contiene, perche' ogni scambio che lo attraversa crea un range
    almeno pari allo spread. `stima / range_mediano > ~0.5` e' impossibile.
      · ES daily : 62bp / 120bp = 0.52  → A la prende (53 tick), B al limite
      · BTC 1-min: 17.87bp / 6.15bp = 2.9 → B lo prende, A no (tick inutile)
    Le due guardie coprono casi DISGIUNTI. Servono entrambe, e il verdetto e' l'OR.

    Passa `median_bar_range_rel` (mediana di (high-low)/mid sulle stesse barre) per attivare
    la guardia B. Senza quello si ottiene solo la A, e il risultato e' marcato `partial`.
    """
    if spread_rel is None or not mid_price or not tick_size:
        return {"plausible": None, "note": "not computable"}
    ticks = (spread_rel * mid_price) / tick_size
    # APPLICABILITA' della guardia A, scoperta dai test (2026-07-30): se il tick e' cosi'
    # fine che anche uno spread di 1bp — plausibile su qualunque asset liquido — implica piu'
    # di max_ticks, allora il conteggio in tick non discrimina NULLA e usarlo per respingere
    # produce FALSI POSITIVI. Su BTC: 1bp di $58 000 = $5.80 = 580 tick da $0.01.
    ticks_at_1bp = (1e-4 * mid_price) / tick_size
    a_applicable = ticks_at_1bp <= max_ticks
    a_ok = (ticks <= max_ticks) if a_applicable else None
    out = {"implied_ticks": ticks, "max_ticks": max_ticks, "guard_a_ok": a_ok,
           "guard_a_applicable": bool(a_applicable), "ticks_at_1bp": ticks_at_1bp}
    b_ok = None
    if median_bar_range_rel and median_bar_range_rel > 0:
        share = spread_rel / median_bar_range_rel
        b_ok = share <= max_range_share
        out.update({"range_share": share, "max_range_share": max_range_share,
                    "guard_b_ok": b_ok})
    checks = [v for v in (a_ok, b_ok) if v is not None]
    out["plausible"] = bool(all(checks)) if checks else None
    out["coverage"] = ("full" if (a_ok is not None and b_ok is not None) else
                       ("B only" if b_ok is not None else
                        ("A only" if a_ok is not None else "none")))
    if not checks:
        out["note"] = ("NON VERIFICABILE: la guardia A non e' applicabile (tick troppo fine: "
                       f"1bp = {ticks_at_1bp:.0f} tick) e non hai passato "
                       "`median_bar_range_rel` per la B. Nessun giudizio.")
    elif b_ok is False:
        out["note"] = (f"NON CREDIBILE (guardia B): la stima e' {out['range_share']:.1f}x il "
                       "range mediano della barra. Fisicamente impossibile — stai misurando "
                       "volatilita', non spread.")
    elif a_ok is False:
        out["note"] = (f"NON CREDIBILE (guardia A): {ticks:.0f} tick impliciti. Barre troppo "
                       "grossolane per uno stimatore di microstruttura.")
    else:
        out["note"] = ("credibile su " + out["coverage"] +
                       ("" if out["coverage"] == "full" else
                        " — ⚠ una sola guardia attiva, l'altra non applicabile o non fornita"))
    return out


def sign_flip_parameter(f, lo, hi, tol=1e-9, max_iter=200):
    """Il valore del parametro a cui la CONCLUSIONE cambia segno. Bisezione su f(p).

    A cosa serve, ed e' la lezione piu' costosa del 2026-07-30: quando un verdetto dipende da
    un numero che non hai misurato tu, la cosa da pubblicare non e' il verdetto — e' la
    SOGLIA DI RIBALTAMENTO. "Positivo se il penalty sta sotto 40.9c" e' una consegna;
    "positivo" con un penalty preso in prestito non lo e'.

    Ritorna None se f non cambia segno in [lo, hi] (nessun ribaltamento nell'intervallo:
    la conclusione e' robusta li' dentro, e QUELLO e' il risultato da riportare).
    """
    fl, fh = f(lo), f(hi)
    if fl == 0:
        return lo
    if fh == 0:
        return hi
    if (fl > 0) == (fh > 0):
        return None
    a, b = lo, hi
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(b - a) < tol or fm == 0:
            return m
        if (fm > 0) == (fl > 0):
            a, fl = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def borrowed_parameter(name, value, source, f, lo, hi, measured=None):
    """Registra un parametro PRESO IN PRESTITO e dice se puo' decidere il verdetto.

    Il caso reale (2026-07-30): ho accettato un penalty di −9c dal README di un terzo. Era
    l'unico numero che decideva il verdetto. Misurato nel MIO contesto: −44.3c, 4.9x piu'
    grande, e per un ERRORE DI CATEGORIA (le sue posizioni erano chiudibili, le mie erano
    binari tenuti a scadenza senza uscita). Il verdetto si e' ribaltato.

    Regola: se il ribaltamento cade DENTRO l'intervallo plausibile del parametro, il verdetto
    NON si consegna. Si consegna la soglia, e si va a misurare.
    """
    flip = sign_flip_parameter(f, lo, hi)
    v_now = measured if measured is not None else value
    out = {"name": name, "borrowed_value": value, "source": source,
           "measured_value": measured, "used": v_now,
           "plausible_range": (lo, hi), "flip_at": flip,
           "verdict_at_used": ("positive" if f(v_now) > 0 else "negative")}
    if measured is None:
        out["status"] = "⚠ IN PRESTITO, NON MISURATO"
        if flip is not None:
            out["deliverable"] = (f"NON consegnare un verdetto: si ribalta a {name}="
                                 f"{flip:.4g}, dentro l'intervallo plausibile "
                                 f"[{lo:.4g}, {hi:.4g}]. Consegna la soglia e misura.")
        else:
            out["deliverable"] = (f"verdetto consegnabile: nessun ribaltamento per {name} in "
                                 f"[{lo:.4g}, {hi:.4g}], quindi il valore in prestito non "
                                 "puo' cambiare la conclusione.")
    else:
        out["status"] = "misurato nel contesto proprio"
        out["deliverable"] = f"verdetto basato su {name}={measured:.4g} misurato."
        if abs(measured - value) > 0.5 * abs(value):
            out["note"] = (f"⚠ il valore misurato differisce dal prestito di "
                           f"{abs(measured-value)/max(abs(value),1e-12):.1f}x — controlla se "
                           "c'e' una DIFFERENZA DI CATEGORIA fra i due contesti, non solo "
                           "di grandezza.")
    return out


def bounce_check(closes, highs, lows, tick_size, mid_price=None):
    """Is measured short-horizon mean reversion a signal, or is it the bid-ask bounce?

    THE TEST THAT SETTLES IT. Estimate the implied spread from the same data. If the
    implied spread is of the order of one tick, the negative autocorrelation you measured
    IS the cost you would pay trying to harvest it — you would be buying at the ask and
    selling at the bid to capture a reversion that exists only between those two prices.

    Returns a dict with both estimators, the ratio to the real tick, and a verdict.
    """
    roll_rel = roll_spread_relative(closes)
    cs = corwin_schultz_median(highs, lows, window=1)
    mid = mid_price or (sorted(closes)[len(closes) // 2] if closes else None)
    out = {"roll_relative": roll_rel, "corwin_schultz_relative": cs,
           "tick_size": tick_size, "mid_price": mid}
    for name, rel in (("roll", roll_rel), ("corwin_schultz", cs)):
        if rel is None or not mid or not tick_size:
            out[f"{name}_ticks"] = None
            continue
        out[f"{name}_ticks"] = (rel * mid) / tick_size
    for name, rel in (("roll", roll_rel), ("corwin_schultz", cs)):
        out[f"{name}_plausibility"] = plausibility_flag(rel, mid, tick_size)
    if any(out[f"{n}_plausibility"].get("plausible") is False
           for n in ("roll", "corwin_schultz")):
        out["verdict"] = ("STIMA NON CREDIBILE a questa frequenza — vedi *_plausibility. "
                          "Gli stimatori di spread richiedono barre fini; su daily di "
                          "strumenti liquidi misurano volatilita', non spread.")
        return out
    ticks = [v for k, v in out.items() if k.endswith("_ticks") and v is not None]
    if not ticks:
        out["verdict"] = "not computable"
    else:
        r = sum(ticks) / len(ticks)
        out["implied_ticks_mean"] = r
        out["verdict"] = (
            "reversion is CONSISTENT WITH BID-ASK BOUNCE — not harvestable" if 0.5 <= r <= 2.0
            else "implied spread exceeds the bounce — investigate liquidity, not signal" if r > 2.0
            else "implied spread below one tick — no meaningful reversion measured")
    return out


# ---------------------------------------------------------------- 53. the cost cliff


def cost_cliff(gross_sharpe, trades_per_year, cost_per_trade, sd_per_trade,
               periods_per_year=252):
    """How much of a gross Sharpe survives a given per-trade cost.

    The arithmetic that kills most intraday claims: per-trade friction is fixed while the
    predictable component shrinks with horizon, so the SAME alpha is valuable at low
    frequency and worthless at high frequency.

    cost_per_trade and sd_per_trade must be in the SAME units (both bp, or both currency).
    Returns gross/net Sharpe and the fraction of the edge consumed.
    """
    if sd_per_trade <= 0 or trades_per_year <= 0:
        return {"error": "sd_per_trade and trades_per_year must be positive"}
    mean_per_trade = gross_sharpe * sd_per_trade / math.sqrt(trades_per_year)
    net_mean = mean_per_trade - cost_per_trade
    net_sharpe = net_mean / sd_per_trade * math.sqrt(trades_per_year)
    return {"gross_sharpe": gross_sharpe, "net_sharpe": net_sharpe,
            "mean_per_trade": mean_per_trade, "cost_per_trade": cost_per_trade,
            "cost_share_of_edge_pct": (100.0 * cost_per_trade / mean_per_trade
                                       if mean_per_trade > 0 else float("inf")),
            "trades_per_year": trades_per_year,
            "annual_cost": cost_per_trade * trades_per_year,
            "survives": net_sharpe > 0}


def max_trade_frequency(gross_sharpe, target_net_sharpe, cost_per_trade, sd_per_trade):
    """The trade-frequency ceiling implied by a target net Sharpe.

    Solves for the largest trades_per_year at which net Sharpe still reaches the target.
    Use it BEFORE building an intraday system: if the ceiling is below the frequency your
    idea requires, the idea is dead on arrival and no backtest is needed.
    """
    if gross_sharpe <= target_net_sharpe or cost_per_trade <= 0:
        return float("inf")
    n = ((gross_sharpe - target_net_sharpe) * sd_per_trade / cost_per_trade) ** 2
    return n


def required_edge_per_trade(cost_per_trade, cost_ratio_target=0.20):
    """Minimum gross edge per trade so that costs are at most `cost_ratio_target` of it.

    A rule of thumb worth stating out loud: if round-trip cost exceeds ~20% of the target,
    the strategy is an execution problem wearing a strategy costume.
    """
    return cost_per_trade / cost_ratio_target


# ---------------------------------------------------------------- 54. where did the money go


def decompose_slippage(theoretical_entry, actual_entry, theoretical_exit, actual_exit,
                       side, commission_per_round_trip=0.0, multiplier=1.0):
    """Split the gap between theoretical and realised P&L into ENTRY slip, EXIT slip, commission.

    This is the diagnostic to run the first time realised P&L disagrees with the simulator.
    The three components have three different fixes:
      - entry slippage   -> order type, limit vs stop, queue position
      - exit slippage    -> same, plus whether the exit is a scheduled or a momentum exit
      - commission       -> broker/venue, or size per trade
    Confusing them wastes months. side is +1 for long, -1 for short.
    """
    s = 1.0 if side >= 0 else -1.0
    theo = s * (theoretical_exit - theoretical_entry) * multiplier
    real = s * (actual_exit - actual_entry) * multiplier - commission_per_round_trip
    entry_slip = -s * (actual_entry - theoretical_entry) * multiplier
    exit_slip = s * (actual_exit - theoretical_exit) * multiplier
    residual = real - theo - entry_slip - exit_slip + commission_per_round_trip
    return {"theoretical_pnl": theo, "realised_pnl": real, "total_gap": real - theo,
            "entry_slippage": entry_slip, "exit_slippage": exit_slip,
            "commission": -commission_per_round_trip, "unexplained_residual": residual,
            "gap_pct_of_theoretical": (100.0 * (real - theo) / abs(theo)) if theo else float("nan")}


def cost_stack(venue="retail_futures"):
    """Reference cost stacks. THESE ARE PRIORS, NOT MEASUREMENTS — always mark them so.

    Orders of magnitude only, to sanity-check a measured number rather than replace it.
    Round-trip, per contract, excluding market impact.
    """
    stacks = {
        "retail_futures":      {"commission": 4.50, "spread_ticks": 1.0,
                                "note": "discount FCM, e-mini/micro"},
        "prop_firm_futures":   {"commission": 3.00, "spread_ticks": 1.0,
                                "note": "prop account; check for monthly/activation fees separately"},
        "institutional_futures": {"commission": 0.50, "spread_ticks": 0.5,
                                  "note": "clearing member, passive fills achievable"},
        "hft_futures":         {"commission": 0.20, "spread_ticks": 0.0,
                                "note": "maker rebates possible; spread earned, not paid"},
    }
    return stacks.get(venue, stacks["retail_futures"])


def report(closes, highs, lows, tick_size, label="instrument"):
    """One-call summary: both spread estimators, bounce verdict, implied cost."""
    out = {"label": label, "n": len(closes)}
    out.update(bounce_check(closes, highs, lows, tick_size))
    lines = [f"{label}: n={len(closes)}"]
    for k in ("roll_relative", "corwin_schultz_relative"):
        v = out.get(k)
        lines.append(f"  {k:26s} {'n/a' if v is None else f'{1e4*v:8.2f} bp'}")
    for k in ("roll_ticks", "corwin_schultz_ticks"):
        v = out.get(k)
        lines.append(f"  {k:26s} {'n/a' if v is None else f'{v:8.2f} ticks'}")
    lines.append(f"  verdict: {out.get('verdict')}")
    out["text"] = "\n".join(lines)
    return out


if __name__ == "__main__":
    random.seed(7)
    px, p = [], 100.0
    for _ in range(4000):
        p *= math.exp(random.gauss(0, 0.001))
        px.append(p)
    obs = [p + (0.01 if random.random() < 0.5 else -0.01) for p in px]
    hi = [max(a, b) + 0.01 for a, b in zip(obs, obs[1:])]
    lo = [min(a, b) - 0.01 for a, b in zip(obs, obs[1:])]
    print(report(obs[:-1], hi, lo, tick_size=0.02, label="synthetic (true spread 0.02)")["text"])

# ---------------------------------------------------------------- 60. operabilità dell'ingresso


def unfillable_share(level, open_px, side, pnl=None):
    """L'ingresso AL LIVELLO era possibile, o il prezzo ci era gia' passato?

    IL CASO REALE, 2026-08-03. Breakout del range di Londra all'apertura di New York:
    t = +4.72 su NQ e +13.11 su ORO, 16 anni, ~4 000 operazioni, 16/17 anni positivi,
    jackknife 20/20, nessuna concentrazione, replicata su strumento indipendente.
    Avrebbe superato OGNI controllo statistico del protocollo.

    L'ha uccisa un controllo di ESECUZIONE: nel 36% dei casi su NQ e nel 44% su ORO la
    sessione apriva GIA' oltre il livello d'ingresso. Su quelle giornate entrare al livello
    e' impossibile — il prezzo ci e' passato durante la notte — e il backtest contabilizzava
    il gap notturno come profitto di un breakout mai eseguibile. Quelle operazioni
    inesistenti valevano il **106%** e il **96%** del guadagno totale. Sulle rotture vere
    intraday: t -0.35 e +0.69. Niente.

    > Nessuna quantita' di rigore statistico prende un'ipotesi di RIEMPIMENTO sbagliata.
    > Prima di credere a un backtest con ingressi a livello (breakout, limit, stop, ritocco
    > di un massimo/minimo), chiedersi sempre: *"a quel prezzo ci potevo davvero entrare, o
    > il mercato ci era gia' oltre quando la sessione e' iniziata?"*

    level    prezzo d'ingresso teorico (lista o scalare)
    open_px  prezzo alla PRIMA barra della finestra in cui si opera
    side     +1 long (si compra sopra il livello) · -1 short
    pnl      opzionale: P&L per operazione, per attribuire quanto del guadagno viene dalle
             operazioni non eseguibili — che e' il numero che decide
    """
    lv = list(level) if hasattr(level, "__iter__") else [level]
    op = list(open_px) if hasattr(open_px, "__iter__") else [open_px]
    sd = list(side) if hasattr(side, "__iter__") else [side] * len(lv)
    if not (len(lv) == len(op) == len(sd)):
        raise ValueError(f"lunghezze diverse: level {len(lv)}, open {len(op)}, side {len(sd)}")
    bad = [(o > l) if s > 0 else (o < l) for l, o, s in zip(lv, op, sd)]
    n = len(lv)
    out = {"n": n, "n_unfillable": sum(bad),
           "share_unfillable": (sum(bad) / n) if n else float("nan")}
    if pnl is not None:
        p = list(pnl)
        if len(p) != n:
            raise ValueError(f"pnl ha {len(p)} elementi, le operazioni sono {n}")
        tot = sum(p)
        pb = sum(v for v, b in zip(p, bad) if b)
        out["pnl_unfillable"] = pb
        out["pnl_fillable"] = tot - pb
        out["share_of_profit_unfillable"] = (pb / tot) if tot else float("nan")
        out["flag"] = bool(tot and (pb / tot) > 0.30)
        out["note"] = ("⚠ oltre il 30% del guadagno viene da ingressi IMPOSSIBILI: il risultato "
                       "misura il gap, non la strategia. Rifare col filtro di operabilita'."
                       if out["flag"] else
                       "gli ingressi non eseguibili non dominano il risultato")
    return out
