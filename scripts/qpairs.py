#!/usr/bin/env python3
"""qpairs — mean reversion, cointegration and pair diagnostics. STDLIB ONLY.

For any strategy whose thesis is "this spread comes back". Answers, in the order that
matters, the three questions that decide whether such a strategy is even worth backtesting:

  1. Does the series mean-revert at all?          -> hurst, adf_stat
  2. HOW FAST does it revert?                     -> ou_half_life        <- decides everything
  3. Is the relationship stable, or drifting?     -> engle_granger, kalman_hedge_ratio

WHY HALF-LIFE COMES BEFORE THE BACKTEST
Half-life sets the holding period, and the holding period sets the trade count, and the trade
count times the round-trip cost is the hurdle. A spread with a 3-day half-life traded on daily
bars gives ~80 round trips a year; the same spread with a 200-day half-life gives ~1 and can
never pay for a research programme. Compute it FIRST: it tells you the strategy's economics
before you write a single line of backtest.

Added after studying (see references/08-microstructure-costs.md):
  - cantaro86/Financial-Models-Numerical-Methods (OU process, Kalman regression)
  - the classical Engle-Granger two-step procedure
"""
import math

# ---------------------------------------------------------------- 55. OU / half-life


def ou_half_life(series):
    """Half-life of mean reversion from an AR(1) fit. THE most useful single number here.

    Fits  d_t = a + lam * s_{t-1} + eps  and returns  -ln(2)/lam.

    Reading it:
      lam >= 0        -> no reversion (returns inf). The thesis is wrong, stop here.
      half-life ~ 1-5 bars   -> reverts fast; likely inside the spread, check qmicro.bounce_check
      half-life ~ 5-50 bars  -> the tradable window
      half-life > 100 bars   -> too slow: few trades, huge exposure to the relationship breaking

    FAILURE MODE: an AR(1) fitted on a NON-stationary spread produces a finite, comfortable
    half-life that means nothing. Always pair it with adf_stat/hurst — half-life is only
    interpretable once stationarity is not rejected.
    """
    n = len(series)
    if n < 20:
        return {"half_life": float("nan"), "lambda": float("nan"), "n": n,
                "note": "need at least 20 observations"}
    y = [series[i] - series[i - 1] for i in range(1, n)]
    x = series[:-1]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxx = sum((v - mx) ** 2 for v in x)
    if sxx <= 0:
        return {"half_life": float("inf"), "lambda": 0.0, "n": n, "note": "degenerate"}
    lam = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx
    hl = (-math.log(2.0) / lam) if lam < 0 else float("inf")
    return {"half_life": hl, "lambda": lam, "n": n,
            "reverts": lam < 0,
            "reading": ("no mean reversion (lambda >= 0)" if lam >= 0
                        else "very fast — check it is not just the bid-ask bounce" if hl < 5
                        else "tradable window" if hl <= 50
                        else "too slow to trade cheaply")}


def ou_fit(series, dt=1.0):
    """Full Ornstein-Uhlenbeck fit: ds = theta*(mu - s)*dt + sigma*dW.

    Returns theta (speed), mu (level), sigma (diffusion) and the half-life. Useful when you
    want the equilibrium level, not just the speed — e.g. to size the entry threshold in
    units of the process's own stationary sd, which is sigma/sqrt(2*theta).
    """
    h = ou_half_life(series)
    lam = h["lambda"]
    if not (lam < 0):
        return {**h, "theta": float("nan"), "mu": float("nan"), "sigma": float("nan"),
                "stationary_sd": float("nan")}
    theta = -lam / dt
    n = len(series)
    y = [series[i] - series[i - 1] for i in range(1, n)]
    x = series[:-1]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    a = my - lam * mx
    mu = -a / lam
    resid = [(b - (a + lam * v)) for v, b in zip(x, y)]
    mr = sum(resid) / len(resid)
    var = sum((r - mr) ** 2 for r in resid) / max(1, len(resid) - 2)
    sigma = math.sqrt(var / dt)
    return {**h, "theta": theta, "mu": mu, "sigma": sigma,
            "stationary_sd": sigma / math.sqrt(2 * theta) if theta > 0 else float("nan")}


# ---------------------------------------------------------------- 56. Hurst exponent


def hurst(series, min_lag=2, max_lag=None):
    """Hurst exponent by the variance-of-differences method.

        H < 0.5  mean reverting      H = 0.5  random walk      H > 0.5  trending

    Cheap, model-free, and the right FIRST look at any series. But it is a weak instrument:
    on a few hundred points its standard error is large enough that 0.45 and 0.55 are not
    distinguishable. Use it to reject the obviously hopeless, never to confirm.

    DRIFT WARNING: this uses the second moment about zero, per the definition
    E[|X(t+tau)-X(t)|^2] ~ tau^(2H). A series with a strong deterministic drift will return
    H near 1 — which is CORRECT (it is trending) but is probably not what you meant if you
    wanted the Hurst of the fluctuations around the drift. De-trend first in that case.
    """
    n = len(series)
    max_lag = max_lag or max(min_lag + 2, min(100, n // 4))
    if n < 40 or max_lag <= min_lag:
        return {"hurst": float("nan"), "n": n, "note": "need at least 40 observations"}
    xs, ys = [], []
    for lag in range(min_lag, max_lag):
        d = [series[i] - series[i - lag] for i in range(lag, n)]
        if len(d) < 8:
            continue
        # MOMENTO SECONDO ATTORNO A ZERO, non varianza attorno alla media.
        # La definizione e' E[|X(t+tau) - X(t)|^2] ~ tau^(2H). Usare la varianza
        # sottrarrebbe la deriva, e un trend puro darebbe H~0 invece di H~1.
        m2 = sum(v * v for v in d) / len(d)
        if m2 <= 0:
            continue
        xs.append(math.log(lag)); ys.append(0.5 * math.log(m2))
    if len(xs) < 4:
        return {"hurst": float("nan"), "n": n, "note": "not enough usable lags"}
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((v - mx) ** 2 for v in xs)
    h = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sxx if sxx else float("nan")
    return {"hurst": h, "n": n, "lags_used": len(xs),
            "reading": ("mean reverting" if h < 0.45 else
                        "trending" if h > 0.55 else "indistinguishable from a random walk")}


# ---------------------------------------------------------------- 57. ADF / cointegration

# Dickey-Fuller critical values, constant no trend, large sample (Fuller 1976 tables).
_ADF_CV = {"1%": -3.43, "5%": -2.86, "10%": -2.57}
# Engle-Granger residual critical values, 2 variables (Phillips-Ouliaris / MacKinnon).
# MORE NEGATIVE than plain ADF because the residual comes from an estimated regression:
# using the plain ADF values here is a classic way to "find" cointegration that is not there.
_EG_CV = {"1%": -3.90, "5%": -3.34, "10%": -3.04}


def adf_stat(series, lags=1):
    """Augmented Dickey-Fuller t-statistic (constant, no trend). More negative = more stationary.

    Returns the statistic plus the critical values, and does NOT convert to a p-value:
    the distribution is non-standard and a fabricated p-value would invite exactly the kind
    of false precision this package exists to prevent. Compare to the printed thresholds.
    """
    n = len(series)
    if n < 20 + lags:
        return {"adf": float("nan"), "n": n, "note": "sample too small"}
    dy = [series[i] - series[i - 1] for i in range(1, n)]
    rows, target = [], []
    for t in range(lags, len(dy)):
        r = [1.0, series[t]]                      # costante + livello ritardato
        r += [dy[t - k] for k in range(1, lags + 1)]
        rows.append(r); target.append(dy[t])
    k = len(rows[0])
    # normal equations, Gaussian elimination (k is tiny)
    A = [[sum(rows[i][a] * rows[i][b] for i in range(len(rows))) for b in range(k)] + \
         [sum(rows[i][a] * target[i] for i in range(len(rows)))] for a in range(k)]
    for c in range(k):
        p = max(range(c, k), key=lambda r_: abs(A[r_][c]))
        if abs(A[p][c]) < 1e-14:
            return {"adf": float("nan"), "n": n, "note": "singular design"}
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r_ in range(k):
            if r_ != c and A[r_][c]:
                f = A[r_][c]
                A[r_] = [v - f * w for v, w in zip(A[r_], A[c])]
    beta = [A[i][k] for i in range(k)]
    resid = [target[i] - sum(beta[j] * rows[i][j] for j in range(k)) for i in range(len(rows))]
    dof = len(rows) - k
    if dof <= 0:
        return {"adf": float("nan"), "n": n, "note": "not enough degrees of freedom"}
    s2 = sum(r * r for r in resid) / dof
    # varianza del coefficiente sul livello ritardato = s2 * (X'X)^-1[1][1]
    B = [[sum(rows[i][a] * rows[i][b] for i in range(len(rows))) for b in range(k)] for a in range(k)]
    inv = [[1.0 if i == j else 0.0 for j in range(k)] for i in range(k)]
    for c in range(k):
        p = max(range(c, k), key=lambda r_: abs(B[r_][c]))
        B[c], B[p] = B[p], B[c]; inv[c], inv[p] = inv[p], inv[c]
        pv = B[c][c]
        if abs(pv) < 1e-14:
            return {"adf": float("nan"), "n": n, "note": "singular"}
        B[c] = [v / pv for v in B[c]]; inv[c] = [v / pv for v in inv[c]]
        for r_ in range(k):
            if r_ != c and B[r_][c]:
                f = B[r_][c]
                B[r_] = [v - f * w for v, w in zip(B[r_], B[c])]
                inv[r_] = [v - f * w for v, w in zip(inv[r_], inv[c])]
    se = math.sqrt(max(1e-300, s2 * inv[1][1]))
    stat = beta[1] / se
    return {"adf": stat, "n": n, "lags": lags, "critical_values": dict(_ADF_CV),
            "stationary_at_5pct": stat < _ADF_CV["5%"]}


def engle_granger(y, x, lags=1):
    """Engle-Granger two-step cointegration test. Returns hedge ratio, spread and residual ADF.

    Step 1: OLS  y = a + b*x        (b is the hedge ratio)
    Step 2: ADF on the residual, against the ENGLE-GRANGER critical values, not the plain
            ADF ones — the residual comes from an estimated regression, so the null
            distribution shifts. Using plain ADF values inflates apparent cointegration.

    ORDER MATTERS: regressing y on x and x on y gives different hedge ratios and can give
    different verdicts. Report both directions when they disagree; a relationship that is
    cointegrated in only one direction is fragile.
    """
    n = min(len(y), len(x))
    if n < 30:
        return {"note": "need at least 30 observations", "n": n}
    y, x = list(y[:n]), list(x[:n])
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((v - mx) ** 2 for v in x)
    if sxx <= 0:
        return {"note": "x is constant", "n": n}
    b = sum((a - mx) * (c - my) for a, c in zip(x, y)) / sxx
    a0 = my - b * mx
    spread = [y[i] - (a0 + b * x[i]) for i in range(n)]
    ad = adf_stat(spread, lags=lags)
    hl = ou_half_life(spread)
    stat = ad.get("adf", float("nan"))
    return {"hedge_ratio": b, "intercept": a0, "spread": spread,
            "adf": stat, "critical_values": dict(_EG_CV),
            "cointegrated_at_5pct": (stat < _EG_CV["5%"]) if stat == stat else False,
            "half_life": hl["half_life"], "hurst": hurst(spread).get("hurst"),
            "n": n,
            "verdict": ("cointegrated at 5% (Engle-Granger CVs)" if stat == stat and stat < _EG_CV["5%"]
                        else "NOT cointegrated at 5% — do not trade the spread")}


# ---------------------------------------------------------------- 58. Kalman hedge ratio


def kalman_hedge_ratio(y, x, delta=1e-4, var_eps=1.0):
    """Time-varying hedge ratio via a Kalman filter (random-walk state on [alpha, beta]).

    A static OLS hedge ratio assumes the relationship never changes. It always does. The
    Kalman version lets beta drift, which is both more honest and more dangerous: a filter
    that adapts fast enough will track ANY relationship, including one that has genuinely
    broken, and will keep telling you the spread is about to revert while it never does.

    `delta` controls how fast beta moves (smaller = stiffer). Report the path of beta, not
    just its last value: a beta that wanders monotonically is a relationship dying, not a
    spread oscillating.
    """
    n = min(len(y), len(x))
    if n < 10:
        return {"note": "need at least 10 observations", "n": n}
    vw = delta / (1.0 - delta)
    state = [0.0, 1.0]                                   # [alpha, beta]
    P = [[1.0, 0.0], [0.0, 1.0]]
    betas, alphas, errs = [], [], []
    for t in range(n):
        # predizione: random walk sullo stato
        for i in range(2):
            for j in range(2):
                P[i][j] += vw if i == j else 0.0
        H = [1.0, x[t]]
        yhat = state[0] * H[0] + state[1] * H[1]
        e = y[t] - yhat
        PH = [P[0][0] * H[0] + P[0][1] * H[1], P[1][0] * H[0] + P[1][1] * H[1]]
        S = H[0] * PH[0] + H[1] * PH[1] + var_eps
        K = [PH[0] / S, PH[1] / S]
        state = [state[0] + K[0] * e, state[1] + K[1] * e]
        P = [[P[i][j] - K[i] * PH[j] for j in range(2)] for i in range(2)]
        alphas.append(state[0]); betas.append(state[1]); errs.append(e)
    drift = betas[-1] - betas[len(betas) // 2]
    rng = max(betas) - min(betas)
    return {"beta": betas, "alpha": alphas, "residual": errs,
            "beta_final": betas[-1], "beta_range": rng, "beta_drift_second_half": drift,
            "n": n,
            "warning": ("beta drifts monotonically — the relationship may be breaking, "
                        "not oscillating") if rng > 0 and abs(drift) > 0.5 * rng else ""}


# ---------------------------------------------------------------- 59. one-call diagnostic


def pair_report(y, x, label_y="y", label_x="x", cost_per_round_trip_bp=None,
                periods_per_year=252):
    """Everything you need BEFORE backtesting a pair, in the order that can kill it early."""
    eg = engle_granger(y, x)
    if "hedge_ratio" not in eg:
        return {"error": eg.get("note"), "label": f"{label_y}~{label_x}"}
    sp = eg["spread"]
    hl = eg["half_life"]
    m = sum(sp) / len(sp)
    sd = math.sqrt(sum((v - m) ** 2 for v in sp) / max(1, len(sp) - 1))
    out = {"pair": f"{label_y} ~ {label_x}", "hedge_ratio": eg["hedge_ratio"],
           "adf": eg["adf"], "cointegrated_at_5pct": eg["cointegrated_at_5pct"],
           "half_life": hl, "hurst": eg["hurst"], "spread_sd": sd, "n": eg["n"]}
    if hl == hl and hl not in (float("inf"),) and hl > 0:
        out["trades_per_year_est"] = periods_per_year / (2.0 * hl)
        if cost_per_round_trip_bp is not None:
            out["annual_cost_bp"] = out["trades_per_year_est"] * cost_per_round_trip_bp
    lines = [f"{out['pair']}: n={out['n']}",
             f"  hedge ratio      {out['hedge_ratio']:+.4f}",
             f"  ADF (residuo)    {out['adf']:+.3f}   soglia 5% Engle-Granger {_EG_CV['5%']}",
             f"  cointegrata 5%   {out['cointegrated_at_5pct']}",
             f"  half-life        {hl:.1f} periodi",
             f"  Hurst            {out['hurst']:.3f}" if out["hurst"] == out["hurst"] else "  Hurst n/a"]
    if "trades_per_year_est" in out:
        lines.append(f"  scambi/anno ~    {out['trades_per_year_est']:.0f}")
        if "annual_cost_bp" in out:
            lines.append(f"  costo annuo      {out['annual_cost_bp']:.0f} bp")
    out["text"] = "\n".join(lines)
    return out


if __name__ == "__main__":
    import random
    rnd = random.Random(5)
    # coppia cointegrata sintetica: x random walk, y = 2x + spread OU
    x, s, xs, ys = 100.0, 0.0, [], []
    for _ in range(1500):
        x += rnd.gauss(0, 0.5)
        s = 0.92 * s + rnd.gauss(0, 0.4)      # half-life teorica = -ln2/ln(0.92) ~ 8.3
        xs.append(x); ys.append(2.0 * x + s)
    print(pair_report(ys, xs, "y", "x", cost_per_round_trip_bp=2.0)["text"])
    print(f"\n  half-life teorica attesa: {-math.log(2)/math.log(0.92):.1f} periodi")
