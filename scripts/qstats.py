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
    Any winner inside this band is indistinguishable from luck.
    Inverse use: if 0/N beat the baseline while headroom was large, the deficit is STRUCTURAL."""
    if n_variants < 2:
        return 0.0
    return sd_period * math.sqrt(n) * math.sqrt(2 * math.log(n_variants))


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
