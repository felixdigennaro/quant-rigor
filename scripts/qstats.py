#!/usr/bin/env python3
"""qstats — uncertainty toolkit for strategy analysis. STDLIB ONLY (no numpy/scipy needed).

Every function answers one question that decides whether a backtest result means anything.

Usage as a library:
    import sys; sys.path.insert(0, '<skill>/scripts')
    from qstats import describe, paired, bootstrap, jackknife, headroom, required_n, measurability

CLI:
    python3 qstats.py --series 410 -322 234 ...            # one variant
    python3 qstats.py --a 410 -322 ... --b 710 -322 ...    # paired comparison (a=base, b=variant)
    python3 qstats.py --ranges 67 40 108 ... --dist 40     # measurability of a 40pt rule
"""
import argparse
import math
import random
import statistics as st

# ---------------------------------------------------------------- distributions


def _t_log_pdf(x, v):
    """log of the Student-t pdf. Uses lgamma: math.gamma overflows for v > ~340, which means
    any daily series longer than ~1.5 years would crash a naive implementation."""
    return (math.lgamma((v + 1) / 2) - math.lgamma(v / 2) - 0.5 * math.log(v * math.pi)
            - ((v + 1) / 2) * math.log1p(x * x / v))


def _t_pdf(x, v):
    lp = _t_log_pdf(x, v)
    return math.exp(lp) if lp > -700 else 0.0


def t_sf(t, v, steps=4000):
    """One-sided survival P(T > t) for Student-t. No scipy.

    Validated against standard t-tables (max error < 3e-5 for all df incl. 1 and 2).
    Small df needs care: the t has power-law tails (df=1 is Cauchy, tail ~1/x^2), so plain
    truncated integration loses a whole percent. Closed forms for df=1,2; adaptive range plus
    an analytic tail remainder for df>=3.
    """
    t = abs(t)
    if v <= 0:
        return float('nan')
    if abs(v - 1) < 1e-12:                       # Cauchy — exact
        return 0.5 - math.atan(t) / math.pi
    if abs(v - 2) < 1e-12:                       # exact
        return 0.5 * (1 - t / math.sqrt(t * t + 2))
    if v > 2000:                                 # indistinguishable from the normal at this point
        return 0.5 * math.erfc(t / math.sqrt(2.0))
    # df >= 3: integrate to `hi`, then add the analytic power-law tail beyond it.
    hi = max(t * 4.0, t + 60.0, 200.0 / max(1.0, v - 2) + t)
    n = steps if steps % 2 == 0 else steps + 1
    h = (hi - t) / n
    s = _t_pdf(t, v) + _t_pdf(hi, v)
    for i in range(1, n):
        s += _t_pdf(t + i * h, v) * (4 if i % 2 else 2)
    integral = s * h / 3
    # tail: pdf(x) ~ C x^-(v+1) with C = pdf(hi)*hi^(v+1)  =>  int_hi^inf = pdf(hi)*hi/v
    tail = _t_pdf(hi, v) * hi / v
    return max(0.0, min(1.0, integral + tail))


def p_two_sided(t, n):
    """Two-sided p-value for a t-statistic with n observations."""
    if n < 2:
        return float('nan')
    return min(1.0, 2 * t_sf(t, n - 1))


READING = [(1.0, "noise — do not act"),
           (1.6, "suggestive — a hypothesis, not a finding"),
           (2.0, "weak evidence — needs a structural argument"),
           (3.0, "evidence — still one-regime conditional"),
           (float('inf'), "strong — if the sample spans regimes")]


def verdict(t):
    if t is None or (isinstance(t, float) and (math.isnan(t) or math.isinf(t))):
        return "undefined — zero variance or n too small"
    a = abs(t)
    for thr, txt in READING:
        if a < thr:
            return txt
    return READING[-1][1]

# ---------------------------------------------------------------- core stats


def describe(x, label="series"):
    """mean, SD, SE, t, p, and the honest reading."""
    n = len(x)
    m = st.mean(x)
    sd = st.stdev(x) if n > 1 else 0.0
    se = sd / math.sqrt(n) if n > 1 else float('inf')
    t = m / se if (n > 1 and se > 0 and math.isfinite(se)) else float('nan')
    return {"label": label, "n": n, "total": sum(x), "mean": m, "sd": sd, "se": se,
            "t": t, "p": p_two_sided(t, n) if n > 1 else float('nan'),
            "verdict": verdict(t) if n > 1 else "n too small",
            "win_rate": 100 * sum(1 for v in x if v > 0) / n,
            "worst": min(x), "best": max(x),
            "skew": _skew(x)}


def _skew(x):
    n = len(x)
    if n < 3:
        return float('nan')
    m, sd = st.mean(x), st.stdev(x)
    if sd == 0:
        return 0.0
    return (n / ((n - 1) * (n - 2))) * sum(((v - m) / sd) ** 3 for v in x)


def paired(a, b, labels=("base", "variant")):
    """PAIRED difference test — the correct comparison for two variants on the same periods.
    Removes the common market factor; far more powerful than comparing two independent means."""
    if len(a) != len(b):
        raise ValueError("paired() needs equal-length, period-aligned series")
    d = [bi - ai for ai, bi in zip(a, b)]
    out = describe(d, f"{labels[1]} - {labels[0]}")
    out["improved"] = sum(1 for v in d if v > 0)
    out["worsened"] = sum(1 for v in d if v < 0)
    out["unchanged"] = sum(1 for v in d if v == 0)
    out["diffs"] = d
    return out


def bootstrap(x, stat=sum, reps=10000, seed=7, block=1):
    """Bootstrap band. block>1 = moving-block bootstrap (preserves autocorrelation).

    Returns the 5/50/95 percentiles of ``stat`` and ``share_opposite_sign_pct`` — the percentage
    of resamples whose statistic has the opposite sign to the observed one. Units: PERCENT
    (0-100), like every ``_pct`` field in this package. Do not multiply by 100 again.
    """
    rnd = random.Random(seed)
    # COERCED TO A LIST, and this is not cosmetic. The wrap-around branch below relies on
    # `+` meaning CONCATENATION. On a numpy array `+` is elementwise ADDITION: when the two
    # slices have different lengths it raises, and when they happen to have the SAME length
    # it silently adds the numbers and returns a plausible, wrong band. The crash is the
    # lucky case. Found 2026-08-10 by calling bootstrap(block=5) on an ndarray.
    x = list(x)
    n = len(x)
    vals = []
    for _ in range(reps):
        if block <= 1:
            s = [x[rnd.randrange(n)] for _ in range(n)]
        else:
            s = []
            while len(s) < n:
                i = rnd.randrange(n)
                s.extend(x[i:i + block] if i + block <= n else x[i:] + x[:(i + block) % n])
            s = s[:n]
        vals.append(stat(s))
    vals.sort()
    q = lambda p: vals[min(len(vals) - 1, int(p * len(vals)))]
    opp = 100.0 * sum(1 for v in vals if (v < 0) != (stat(x) < 0)) / len(vals)
    # every "share"/"ratio" field in this package is in PERCENT and named _pct.
    return {"p05": q(0.05), "p50": q(0.50), "p95": q(0.95),
            "share_opposite_sign_pct": opp,
            "share_opposite_sign": opp}   # deprecated alias, same units (percent)


def jackknife(x, stat=sum):
    """Drop-one stability: does the sign survive every single deletion?"""
    x = list(x)          # `x[:i] + x[i+1:]` must CONCATENATE — see the note in bootstrap()
    n = len(x)
    outs = []
    for i in range(n):
        outs.append((i, stat(x[:i] + x[i + 1:])))
    full = stat(x)
    same = sum(1 for _, v in outs if (v > 0) == (full > 0))
    worst = min(outs, key=lambda kv: kv[1])
    best = max(outs, key=lambda kv: kv[1])
    return {"full": full, "sign_survives": f"{same}/{n}",
            "worst_drop_index": worst[0], "worst_value": worst[1],
            "best_drop_index": best[0], "best_value": best[1]}


def concentration(x, top=3):
    """How concentrated is the result in a few periods?

    On a signed P&L series a single "share of total" is NOT interpretable: the k largest
    periods by magnitude can sum to the opposite sign of the total, giving shares below 0
    or above 100 that read like a bug. So three distinct quantities are returned:

    - ``gain_share_pct``  : top-k gains / sum of all gains   (0-100, the honest concentration)
    - ``loss_share_pct``  : top-k losses / sum of all losses  (0-100)
    - ``abs_share_pct``   : top-k by |value| / total          (signed, may exceed 100 or go
                            negative — meaningful only as a *fragility* flag)

    ``interpretable`` is False when the total is small relative to the biggest periods, which
    is exactly when a share-of-total figure must not be quoted.
    """
    gains = sorted((v for v in x if v > 0), reverse=True)
    losses = sorted((v for v in x if v < 0))
    by_abs = sorted(x, key=abs, reverse=True)[:top]
    tot, gsum, lsum = sum(x), sum(gains), sum(losses)
    pct = lambda a, b: (100.0 * a / b) if b else float('nan')
    return {"top_k": top,
            "top_k_sum": sum(by_abs),
            "top_gains_sum": sum(gains[:top]), "top_losses_sum": sum(losses[:top]),
            "gain_share_pct": pct(sum(gains[:top]), gsum),
            "loss_share_pct": pct(sum(losses[:top]), lsum),
            "abs_share_pct": pct(sum(by_abs), tot),
            "interpretable": bool(tot) and abs(tot) > abs(sum(by_abs)),
            # kept for backward compatibility with existing callers
            "share_pct": pct(sum(by_abs), tot)}

# ---------------------------------------------------------------- search / power


def headroom(sd_period, n, n_variants):
    """Best-of-N noise headroom: sigma_total * sqrt(2 ln N).

    WHAT THIS IS, STATED EXACTLY (measured, not asserted): sqrt(2 ln N) is the leading
    term of the asymptotic for the max of N iid normals. It is NOT the mean of that max,
    and it is NOT a fixed significance level -- the quantile it corresponds to DRIFTS
    with N:

        N=2 -> 22.5%   N=8 -> 15.4%   N=100 -> 11.3%   N=462 -> 10.1%   N=1e4 -> 8.5%

    So "beats the headroom" is a WEAK screen (roughly the 10-20% level), not a 5% test.
    A winner just above it is still beaten by noise in one search out of six to ten.
    For a stated significance level use `search_threshold()` below.

    Inverse use, and this one is strong: if 0/N beat the baseline while headroom was
    large, the deficit is STRUCTURAL, not sampling."""
    if n_variants < 2:
        return 0.0
    return sd_period * math.sqrt(n) * math.sqrt(2 * math.log(n_variants))


def _norm_sf(t):
    """P(Z > t) for standard normal, via erfc — accurate deep in the tail where
    1 - Phi(t) would lose every significant digit to cancellation."""
    return 0.5 * math.erfc(t / math.sqrt(2.0))


def _max_exceed_prob(t, n):
    """P(max of n iid standard normals > t) = 1 - Phi(t)^n.

    Computed as -expm1(n * log1p(-sf)) rather than 1 - (1-sf)**n: at t=4, sf is 3e-5 and
    the naive form throws away the answer's leading digits. The whole point of this
    function is the tail, so the tail is where it has to be right."""
    sf = _norm_sf(t)
    if sf >= 1.0:
        return 1.0
    return -math.expm1(n * math.log1p(-sf))


def search_threshold(n_variants, alpha=0.05):
    """EXACT critical t for the best of `n_variants` independent tries, at level `alpha`.

    Solves P(max of N iid standard normals > t) = alpha, i.e. 1 - Phi(t)^N = alpha.
    Use this when you want to say "significant at 5% after searching N" instead of the
    looser `headroom()` screen.

        search_threshold(8)    -> 2.49   (vs sqrt(2 ln 8)   = 2.04, which is the 15.4% level)
        search_threshold(462)  -> 3.69   (vs sqrt(2 ln 462) = 3.50, which is the 10.1% level)

    CAVEAT that cuts both ways and must not be assumed away: this assumes the N tries are
    INDEPENDENT. Real variant grids are correlated, so the effective N is smaller and the
    true threshold is LOWER -- with mean rho 0.3 over 462 alphas the effective breadth
    collapses to ~3.3 and the 5% threshold falls to ~2.16. The correction is large in both
    directions, so measure it (`qportfolio.effective_breadth`) rather than guess."""
    if n_variants < 2:
        return 1.959963984540054      # N=1: il classico 1,96 a due code
    lo, hi = 0.0, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _max_exceed_prob(mid, n_variants) > alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def search_quantile(n_variants):
    """Which significance level does sqrt(2 ln N) actually correspond to, for this N?

    Exists so the drift is visible instead of assumed: callers who reach for headroom()
    can print the level they are really screening at."""
    if n_variants < 2:
        return 1.0
    return _max_exceed_prob(math.sqrt(2 * math.log(n_variants)), n_variants)


def required_n(sd_period, effect_per_period, t_target=2.0):
    """Periods needed to detect `effect` at |t|=t_target."""
    if effect_per_period == 0:
        return float('inf')
    return (t_target * sd_period / abs(effect_per_period)) ** 2


def effective_n(x, max_lag=5):
    """n adjusted for autocorrelation: n / (1 + 2*sum(rho_k)). Overlapping/serial data counts less."""
    n = len(x)
    if n < 4:
        return n
    m = st.mean(x)
    var = sum((v - m) ** 2 for v in x)
    if var == 0:
        return n
    s = 0.0
    for k in range(1, min(max_lag, n - 2) + 1):
        cov = sum((x[i] - m) * (x[i + k] - m) for i in range(n - k))
        s += cov / var
    denom = 1 + 2 * s
    # cap at n: negative autocorrelation can push n_eff above n, which is technically true but
    # misleading in a report — never claim MORE independent observations than you have.
    return min(float(n), n / denom) if denom > 0 else float(n)

# ---------------------------------------------------------------- measurability


def measurability(bar_ranges, rule_distance):
    """THE gate (SKILL.md §2). Is the rule's distance resolvable by this data?"""
    med = st.median(bar_ranges)
    R = rule_distance / med if med else float('inf')
    share_ge = 100 * sum(1 for r in bar_ranges if r >= rule_distance) / len(bar_ranges)
    if R >= 3:
        v = "MEASURABLE"
    elif R >= 1:
        v = "NOT MEASURABLE — results are artefacts"
    else:
        v = "FICTION — the rule lives inside one bar"
    return {"median_bar_range": med, "rule_distance": rule_distance, "ratio_R": R,
            "pct_bars_swallowing_rule": share_ge, "verdict": v,
            "needed_distance_for_R3": 3 * med}


def refinement_factor(coarse_ranges, fine_ranges_per_coarse):
    """Measured path-refinement: sum(fine ranges)/coarse range, vs sqrt(n) random-walk null.
    ~sqrt(n) => path has no intrinsic scale => path-dependent rules never converge."""
    ratios = []
    for c, fines in zip(coarse_ranges, fine_ranges_per_coarse):
        if c > 0:
            ratios.append(sum(fines) / c)
    n_sub = st.mean([len(f) for f in fine_ranges_per_coarse]) if fine_ranges_per_coarse else 1
    return {"median_ratio": st.median(ratios) if ratios else float('nan'),
            "mean_ratio": st.mean(ratios) if ratios else float('nan'),
            "random_walk_null_sqrt_n": math.sqrt(n_sub),
            "interpretation": "≈sqrt(n) => no convergence; only rules >> finest bar range are testable"}

# ---------------------------------------------------------------- drawdown / risk


def equity_from_returns(x, start=0.0):
    eq, acc = [], start
    for v in x:
        acc += v
        eq.append(acc)
    return eq


def drawdown_from_pnl(pnl, start=0.0):
    """Drawdown of a PER-PERIOD P&L series (the usual input for a trade/session list).

    Accumulates first, then measures. Use this whenever your list is "profit of period i",
    not "account value at period i". Passing per-period P&L to :func:`drawdown` measures
    ``max(x) - min(x)`` instead of a drawdown — a silent and severe understatement.
    """
    eq, c = [], start
    for v in pnl:
        c += v
        eq.append(c)
    out = drawdown(eq, _checked=True)
    out["equity"] = eq
    out["total"] = c - start
    return out


def _looks_like_pnl(x):
    """Heuristic: is this a per-period P&L series rather than a cumulative equity curve?

    A P&L series changes sign often and its cumulative sum wanders far from its own values.
    An equity curve does the opposite. Used only to warn, never to silently reinterpret.
    """
    if len(x) < 4:
        return False
    # An equity curve does not change sign at a third of its own steps: that would mean the
    # account repeatedly goes to zero and back. Frequent sign changes = per-period P&L.
    flips = sum(1 for a, b in zip(x, x[1:]) if (a > 0) != (b > 0))
    return flips / (len(x) - 1) > 0.30


def drawdown(equity, _checked=False):
    """Max drawdown of an **equity curve** (cumulative account value per period).

    For a per-period P&L / returns list use :func:`drawdown_from_pnl` instead. If the input
    looks like per-period P&L, the result carries ``warning`` and ``input_looks_like_pnl=True``
    rather than failing silently.
    """
    if not equity:
        return {"series": [], "max_dd": 0.0, "periods_under_water": 0,
                "longest_under_water": 0, "input_looks_like_pnl": False}
    flag = (not _checked) and _looks_like_pnl(equity)
    peak = equity[0]
    dd, worst, under, longest, run = [], 0.0, 0, 0, 0
    for v in equity:
        peak = max(peak, v)
        d = v - peak
        dd.append(d)
        worst = min(worst, d)
        if d < 0:
            run += 1
            under += 1
            longest = max(longest, run)
        else:
            run = 0
    res = {"series": dd, "max_dd": worst, "periods_under_water": under,
           "longest_under_water": longest, "input_looks_like_pnl": flag}
    if flag:
        res["warning"] = ("input looks like per-period P&L, not an equity curve; "
                          "max_dd here is max(x)-min(x), not a drawdown — use drawdown_from_pnl()")
    return res


def var_cvar(x, alpha=0.95):
    """EMPIRICAL VaR and CVaR (never the normal formula — returns aren't normal)."""
    s = sorted(x)
    k = max(0, min(len(s) - 1, int((1 - alpha) * len(s))))
    v = s[k]
    tail = s[:k + 1]
    return {"alpha": alpha, "VaR_empirical": v,
            "CVaR_expected_shortfall": st.mean(tail) if tail else v,
            "worst": s[0]}


def kelly_fraction(mean_per_period, sd_period):
    """f* = mu/sigma^2 as a CEILING. Above it, more size lowers MEDIAN wealth."""
    if sd_period == 0:
        return float('inf')
    return mean_per_period / (sd_period ** 2)

# ---------------------------------------------------------------- reporting


def report(x, label="strategy", n_variants=1, costs_note=None):
    d = describe(x, label)
    bs = bootstrap(x)
    jk = jackknife(x)
    cc = concentration(x)
    dd = drawdown(equity_from_returns(x))
    vc = var_cvar(x)
    eff = effective_n(x)
    hr = headroom(d["sd"], d["n"], n_variants)
    L = [f"=== {label} ===",
         f"  n {d['n']} (effective {eff:.1f}) · total {d['total']:+.1f} · mean {d['mean']:+.2f} ± {d['se']:.2f} (SE)",
         f"  t {d['t']:+.2f} · p {d['p']:.3f} → {d['verdict']}",
         f"  win {d['win_rate']:.0f}% · worst {d['worst']:+.1f} · best {d['best']:+.1f} · skew {d['skew']:+.2f}",
         f"  max DD {dd['max_dd']:+.1f} · periods under water {dd['periods_under_water']}/{d['n']} (longest {dd['longest_under_water']})",
         f"  bootstrap total 5/50/95: {bs['p05']:+.0f} / {bs['p50']:+.0f} / {bs['p95']:+.0f} · opposite sign {bs['share_opposite_sign_pct']:.0f}%",
         f"  jackknife sign survives {jk['sign_survives']} · worst deletion → {jk['worst_value']:+.0f}",
         f"  top-{cc['top_k']} periods = {cc['share_pct']:.0f}% of result",
         f"  empirical VaR95 {vc['VaR_empirical']:+.1f} · CVaR {vc['CVaR_expected_shortfall']:+.1f}",
         f"  Kelly ceiling f* {kelly_fraction(d['mean'], d['sd']):.3f} (use a FRACTION of it)"]
    if n_variants > 1:
        L.append(f"  best-of-{n_variants} search headroom ±{hr:.0f} → anything inside is luck")
    if costs_note:
        L.append(f"  costs: {costs_note}")
    return "\n".join(L)


def _main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", nargs="+", type=float)
    ap.add_argument("--a", nargs="+", type=float, help="baseline per-period series")
    ap.add_argument("--b", nargs="+", type=float, help="variant per-period series")
    ap.add_argument("--variants", type=int, default=1, help="how many variants were searched")
    ap.add_argument("--ranges", nargs="+", type=float, help="bar ranges for measurability")
    ap.add_argument("--dist", type=float, help="rule characteristic distance")
    a = ap.parse_args()
    if a.series:
        print(report(a.series, "series", a.variants))
    if a.a and a.b:
        print(report(a.a, "baseline", a.variants))
        print()
        print(report(a.b, "variant", a.variants))
        print()
        pr = paired(a.a, a.b)
        print("=== PAIRED DIFFERENCE (the decisive test) ===")
        print(f"  mean {pr['mean']:+.2f} ± {pr['se']:.2f} · t {pr['t']:+.2f} · p {pr['p']:.3f} → {pr['verdict']}")
        print(f"  improved {pr['improved']} · worsened {pr['worsened']} · unchanged {pr['unchanged']}")
        jk = jackknife(pr["diffs"])
        print(f"  jackknife on the DIFFERENCE: sign survives {jk['sign_survives']} · worst {jk['worst_value']:+.0f}")
        if a.variants > 1:
            print(f"  headroom ±{headroom(describe(a.a)['sd'], len(a.a), a.variants):.0f}")
    if a.ranges and a.dist:
        m = measurability(a.ranges, a.dist)
        print("=== MEASURABILITY GATE ===")
        for k, v in m.items():
            print(f"  {k}: {v if not isinstance(v, float) else round(v, 2)}")


if __name__ == "__main__":
    _main()

# ---------------------------------------------------------------- 60. soglia raggiungibile


def reachable_threshold(window, threshold, ddof=0):
    """La soglia e' RAGGIUNGIBILE dallo stimatore, o il test non puo' sparare mai?

    IL BUG CHE ESISTE PER PRENDERE (2026-07-30): uno z-score su finestra W, calcolato per un
    punto DENTRO la finestra, ha un massimo ALGEBRICO:
        sd di popolazione (ddof=0):  |z|max = sqrt(W-1)
        sd campionaria    (ddof=1):  |z|max = (W-1)/sqrt(W)
    Con W=4 e ddof=0 il massimo e' 1.732: una soglia |z|>2 e' MATEMATICAMENTE IRRAGGIUNGIBILE.
    La strategia stava in posizione lo 0.0% del tempo, senza sollevare NESSUN errore. Un test
    che non puo' sparare non e' un test superato: e' un test rotto, e assomiglia a un
    risultato negativo.

    Da chiamare PRIMA di ogni backtest che usi una soglia su una statistica normalizzata.
    """
    if window < 2:
        return {"reachable": False, "max_abs": 0.0,
                "note": "finestra < 2: nessuna statistica definita"}
    m = math.sqrt(window - 1) if ddof == 0 else (window - 1) / math.sqrt(window)
    ok = abs(threshold) < m
    need = None
    if not ok:
        # finestra minima che rende raggiungibile la soglia
        w = max(2, int(math.ceil(threshold ** 2 + 1)) if ddof == 0 else 2)
        while True:
            mm = math.sqrt(w - 1) if ddof == 0 else (w - 1) / math.sqrt(w)
            if mm > abs(threshold):
                break
            w += 1
            if w > 10 ** 7:
                break
        need = w
    return {"reachable": bool(ok), "max_abs": m, "threshold": threshold,
            "window": window, "ddof": ddof, "min_window_needed": need,
            "note": ("ok" if ok else
                     f"IRRAGGIUNGIBILE: |z|max={m:.3f} < soglia {abs(threshold):.3f}. "
                     f"La regola non puo' sparare MAI. Serve finestra >= {need}.")}


# ---------------------------------------------------------------- 61. tasso di eventi


def event_rate_check(n_events, years, intended_per_year, tol=2.0, label="regola"):
    """Il tasso di eventi osservato corrisponde a quello che la regola INTENDE?

    DUE BUG PRESI COSI' (2026-07-30): un ORB e una reversione a zona rientravano 5-7 volte al
    giorno invece di una, perche' mancava il flag "gia' operato oggi" — 5 492 e 2 743 cambi di
    posizione l'anno invece di ~730. Nessun errore, solo numeri sbagliati a valle.
    E in senso opposto: uno z-score che non spara mai da 0 eventi.

    Da chiamare subito dopo aver generato i segnali e PRIMA di guardare il P&L.
    """
    rate = n_events / years if years > 0 else float("inf")
    ratio = rate / intended_per_year if intended_per_year > 0 else float("inf")
    ok = (1.0 / tol) <= ratio <= tol
    if n_events == 0:
        note = ("ZERO EVENTI: non e' un risultato negativo, e' una regola che non spara. "
                "Verifica che la soglia sia raggiungibile (reachable_threshold).")
    elif ratio > tol:
        note = (f"{ratio:.1f}x TROPPI eventi: cerca un vincolo mancante (una sola operazione "
                "per giorno/sessione? un flag 'gia' fatto'? un cooldown?).")
    elif ratio < 1.0 / tol:
        note = (f"{ratio:.2f}x troppo POCHI: la condizione e' piu' stretta di quanto credi.")
    else:
        note = "ok"
    return {"label": label, "n_events": n_events, "rate_per_year": rate,
            "intended": intended_per_year, "ratio": ratio, "ok": bool(ok), "note": note}


# ---------------------------------------------------------------- 62. scala logaritmica


def log_to_simple(x):
    """Converte una metrica calcolata in spazio LOG in composta. Sembra banale: non lo e'.

    ERRORE COMMESSO (2026-07-30): ho pubblicato drawdown calcolati su cumsum di log-rendimenti
    stampandoli come percentuali. −0.652 in log e' −48% reale, non −65%. Sovrastimava il
    drawdown di 17 punti su ogni riga di ogni tabella.
    Regola: qualunque numero nato in spazio log si converte PRIMA di pubblicarlo, e la scala
    si dichiara nell'etichetta.
    """
    return math.exp(x) - 1.0

# ---------------------------------------------------------------- 69. anni per decidere


# Decadenza post-pubblicazione misurata da McLean & Pontiff (2016) su 97 anomalie:
# ~26% dovuto a rumore statistico nel campione originale, ~32% al fatto che il mercato la
# arbitraggia dopo la pubblicazione. Il rendimento fuori campione e' circa il 42% di quello
# dichiarato. ⚠ E' un numero PRESO IN PRESTITO da una popolazione di anomalie azionarie
# americane: se la strategia in esame e' di altra natura, va rimisurato o dichiarato tale (G5).
PUBLICATION_DECAY = 0.42


def years_to_decide(sharpe, t_target=2.0, decay=None, years_available=None):
    """Quanti anni servono perche' questo Sharpe si distingua da zero? DA CHIAMARE PRIMA
    di scrivere il backtest, non dopo.

    Da `t = Sharpe x sqrt(anni)` segue `anni = (t_target / Sharpe)^2`. E' la domanda che
    rende molte analisi inutili prima ancora di cominciare: uno Sharpe di 0.359 — la mediana
    delle strategie pubblicate — richiede **31 anni** per t=2, e 176 anni dopo la decadenza
    post-pubblicazione. Scrivere quel backtest e' tempo speso per produrre un numero che non
    puo' decidere niente.

    Misurato il 2026-08-04 su 60 strategie pubblicate con Sharpe dichiarato
    (paperswithbacktest/awesome-systematic-trading):
        Sharpe mediano 0.359 · 9 su 60 gia' NEGATIVO nel paper
        decidibili con 20 anni al valore dichiarato ......... 22 su 60
        decidibili con 20 anni dopo la decadenza ............  0 su 60

    Per il calcolo rigoroso (che corregge per asimmetria, curtosi e numero di tentativi) usare
    `qvalidate.min_track_record_length` e `qvalidate.deflated_sharpe`. Questa e' la versione
    da un rigo, buona per il cancello d'ingresso.
    """
    s = float(sharpe)
    out = {"sharpe": s, "t_target": float(t_target)}
    if s <= 0:
        out["years"] = float("inf")
        out["note"] = "Sharpe non positivo: nessun orizzonte lo rende positivo"
        return out
    out["years"] = (float(t_target) / s) ** 2
    if decay is not None or decay is None:
        d = PUBLICATION_DECAY if decay is None else float(decay)
        sd = s * d
        out["decay_applied"] = d
        out["sharpe_after_decay"] = sd
        out["years_after_decay"] = (float(t_target) / sd) ** 2 if sd > 0 else float("inf")
    if years_available is not None:
        ya = float(years_available)
        out["years_available"] = ya
        out["decidable_at_claimed"] = out["years"] <= ya
        out["decidable_after_decay"] = out.get("years_after_decay", float("inf")) <= ya
        out["max_t_reachable"] = s * math.sqrt(ya)
        out["verdict"] = ("DECIDIBILE" if out["decidable_after_decay"] else
                          ("decidibile solo al valore DICHIARATO, non dopo la decadenza"
                           if out["decidable_at_claimed"] else
                           f"NON DECIDIBILE: servono {out['years']:.0f} anni "
                           f"({out.get('years_after_decay', float('inf')):.0f} dopo decadenza), "
                           f"ne hai {ya:.0f}"))
    return out


def literature_queue(entries, years_available, t_target=2.0, decay=None):
    """Ordina una lista di strategie pubblicate per DECIDIBILITA', non per Sharpe.

    `entries`: lista di dict con almeno `name` e `sharpe`.
    Il Sharpe pubblicato e' il massimo di una ricerca che non hai controllato, in un contesto
    che non hai verificato: serve a ORDINARE la coda, mai a decidere. Quello che decide e'
    se l'orizzonte disponibile basta a distinguerlo da zero.
    """
    out = []
    for e in entries:
        r = years_to_decide(e.get("sharpe", 0.0), t_target=t_target, decay=decay,
                            years_available=years_available)
        out.append(dict(e, **{"years_needed": r["years"],
                              "years_after_decay": r.get("years_after_decay", float("inf")),
                              "decidable": r.get("decidable_at_claimed", False),
                              "decidable_after_decay": r.get("decidable_after_decay", False)}))
    out.sort(key=lambda x: x["years_needed"])
    n_ok = sum(1 for x in out if x["decidable"])
    n_dec = sum(1 for x in out if x["decidable_after_decay"])
    return {"queue": out, "n": len(out), "n_decidable_claimed": n_ok,
            "n_decidable_after_decay": n_dec,
            "note": (f"{n_ok}/{len(out)} decidibili in {years_available:.0f} anni al valore "
                     f"dichiarato, {n_dec}/{len(out)} dopo la decadenza post-pubblicazione. "
                     "Il Sharpe pubblicato ordina la coda; non e' evidenza.")}
