#!/usr/bin/env python3
"""qregime — "no edge works forever": the objection, made testable.

THE OBJECTION, which is serious and usually right in spirit:
> "Test any strategy over 16 years — even a big fund's — and you get a negative result,
>  because no edge survives every cycle and regime. So the full-sample test is the wrong tool."

WHERE IT IS RIGHT. `t = Sharpe·√years` assumes stationarity. An edge alive for 4 of 16 years
gets diluted by averaging over 16. Edge decay is documented (McLean & Pontiff: ~58% of anomaly
returns disappear post-publication). None of this is in dispute.

WHERE THE ARGUMENT BREAKS, and it is measurable:
**Dilution lowers t; it does not flip the sign.** 4 good years + 12 flat years = a POSITIVE
mean with a low t. To get a *negative* full-sample result you need years where the strategy
actually LOSES, not years where it merely stops working. So the objection, taken seriously,
makes three falsifiable predictions:

    P1  contiguous windows with strong performance must exist
    P2  there must be MORE of them than chance produces — otherwise "the good regime" is just
        the maximum of a search over windows, which exists in pure noise too
    P3  the regime must be identifiable IN ADVANCE, not in hindsight

P3 is the one that decides everything. A strategy that works in identifiable regimes is
tradeable; one that works in periods identifiable only afterwards is not.

TWO TRAPS THIS MODULE EXISTS TO AVOID:

1. **Picking "the best window" is a SEARCH.** 16 years scanned with 2-year windows gives ~14
   overlapping positions; the max of 14 correlated windows beats the mean by ~1.5-2σ even with
   no edge. So P1 alone proves nothing, and the null must be **block**-shuffled to preserve
   autocorrelation and volatility clustering — an iid shuffle makes every series look regimed.

2. **A switch that is merely out of the market part of the time looks like skill.** So P3 is
   always run against TWO baselines: always-on, and a random switch with the same duty cycle.
   Beating always-on but not the random switch means the merit is in being flat, not in the rule.

MEASURED, 2026-08-03 (mean-reversion on 16y NQ, 2 805 trading days):
    best 2-year window t = +0.06, p = 0.983 against the block null — 98% of *shuffled* series
    produce a better best-window. Not "an edge diluted across regimes": never an edge.
    And the trailing switch made it worse (-128 → -333), worse than random switching too.
Validated in both directions first: on a synthetic edge confined to a quarter of the sample,
the same test returns p = 0.000 and the switch improves +443 → +690, beating random. A test
that cannot detect a regime-limited edge would be worthless here — see `validate_power()`.
"""
import math

import numpy as np

__all__ = ["block_shuffle", "window_scan", "regime_switch", "regime_report", "validate_power"]


def _t(x):
    x = np.asarray(x, float)
    if len(x) < 10 or x.std(ddof=1) == 0:
        return float("nan")
    return float(x.mean() / (x.std(ddof=1) / math.sqrt(len(x))))


def block_shuffle(x, block, rng):
    """Shuffle in BLOCKS. Destroys the ordering of regimes while preserving autocorrelation
    and volatility clustering. With an iid shuffle the null is too easy to beat and every
    series appears to have regimes — the exact illusion this module is built to prevent."""
    x = np.asarray(x, float)
    n = len(x)
    nb = int(np.ceil(n / block))
    idx = rng.permutation(nb)
    return np.concatenate([x[i * block:(i + 1) * block] for i in idx])[:n]


def window_scan(x, w, block=None, n_null=400, seed=0, t_good=2.0):
    """P1 + P2: are there strong windows, and MORE of them than the block null gives?

    Returns the best window t, the share of windows above `t_good`, and the empirical p of
    each against the block-shuffled null. `block` defaults to w//4 (never smaller than the
    autocorrelation scale you care about).
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 3 * w:
        return {"ok": False, "note": f"servono almeno {3*w} osservazioni, ce ne sono {len(x)}"}
    block = block or max(20, w // 4)
    step = max(1, w // 4)

    def scan(y):
        v = np.array([_t(y[i:i + w]) for i in range(0, len(y) - w + 1, step)])
        return v[np.isfinite(v)]

    real = scan(x)
    if len(real) < 3:
        return {"ok": False, "note": "troppe poche finestre"}
    rng = np.random.default_rng(seed)
    bn, fn = [], []
    for _ in range(n_null):
        v = scan(block_shuffle(x, block, rng))
        if len(v) < 3:
            continue
        bn.append(float(v.max()))
        fn.append(float(np.mean(v > t_good)))
    bn, fn = np.array(bn), np.array(fn)
    best, frac = float(real.max()), float(np.mean(real > t_good))
    return {"ok": True, "n_win": len(real), "window": w, "block": block,
            "best": best, "frac_good": frac,
            "null_best_median": float(np.median(bn)),
            "null_best_q95": float(np.percentile(bn, 95)),
            "null_frac_median": float(np.median(fn)),
            "p_best": float(np.mean(bn >= best)), "p_frac": float(np.mean(fn >= frac))}


def regime_switch(x, look=250, seed=0, detector=None):
    """P3: does a switch using ONLY the past beat always-on AND a random switch?

    `detector`: optional array of the same length, causal, > 0 meaning "on". Defaults to
    "trailing mean of the strategy's own P&L is positive" — the simplest possible detector,
    and the one to beat before proposing a more elaborate one.

    The random control is not decoration. A switch that is simply flat part of the time
    changes the mean without any skill; if the rule does not beat random at the same duty
    cycle, the merit is in being out of the market, not in the timing.
    """
    x = np.asarray(x, float)
    n = len(x)
    if n < look + 60:
        return {"ok": False, "note": f"servono almeno {look+60} osservazioni"}
    on = np.zeros(n, bool)
    if detector is None:
        for i in range(look, n):
            on[i] = x[i - look:i].mean() > 0
    else:
        d = np.asarray(detector, float)
        on[look:] = np.nan_to_num(d[look:]) > 0
    rng = np.random.default_rng(seed)
    k = int(on[look:].sum())
    rnd = np.zeros(n, bool)
    if 0 < k <= n - look:
        rnd[rng.choice(np.arange(look, n), size=k, replace=False)] = True
    sel, srnd, always = x[look:][on[look:]], x[look:][rnd[look:]], x[look:]
    return {"ok": True, "duty": k / (n - look),
            "mean_rule": float(sel.mean()) if len(sel) else float("nan"),
            "t_rule": _t(sel),
            "mean_random": float(srnd.mean()) if len(srnd) else float("nan"),
            "t_random": _t(srnd),
            "mean_always": float(always.mean()), "t_always": _t(always),
            "beats_always": bool(len(sel) and sel.mean() > always.mean()),
            "beats_random": bool(len(sel) and len(srnd) and sel.mean() > srnd.mean())}


def regime_report(x, windows=(250, 500), look=250, n_null=400, seed=0, detector=None):
    """The three predictions, evaluated together, with a verdict DERIVED from them."""
    out = {"n": int(len(x)), "t_full": _t(x), "scans": {}}
    for w in windows:
        out["scans"][w] = window_scan(x, w, n_null=n_null, seed=seed)
    out["switch"] = regime_switch(x, look=look, seed=seed, detector=detector)
    ok = [s for s in out["scans"].values() if s.get("ok")]
    p1 = any(s["best"] > 2.0 for s in ok)
    p2 = any(s["p_best"] < 0.05 for s in ok)
    sw = out["switch"]
    p3 = bool(sw.get("ok") and sw.get("beats_always") and sw.get("beats_random"))
    out["P1_strong_windows_exist"] = p1
    out["P2_more_than_chance"] = p2
    out["P3_identifiable_ex_ante"] = p3
    if p2 and p3:
        v = "ROTAZIONE SOSTENUTA — regimi reali e riconoscibili in anticipo"
    elif p2 and not p3:
        v = "REGIMI REALI MA NON RICONOSCIBILI IN ANTICIPO — non operabile così"
    elif p1 and not p2:
        v = ("FINESTRE BUONE COMPATIBILI COL CASO — il 'buon regime' è il massimo di una "
             "ricerca su finestre, non un fatto")
    else:
        v = "NESSUN REGIME — non è un edge diluito, non c'è mai stato"
    out["verdict"] = v
    return out


def validate_power(n=2800, sd=1000.0, frac=0.25, effect=0.45, seed=1, n_null=200):
    """Does this test have the power to fire? Run it before trusting a negative.

    Builds two synthetic series — one with no edge, one with a REAL edge confined to `frac`
    of the sample — and returns both reports. The regime-limited one must come back with
    p_best ≈ 0 and a switch that beats both baselines; the null one must not. A test that
    cannot separate these two says nothing about your real series.
    """
    rng = np.random.default_rng(seed)
    flat = rng.normal(0, sd, n)
    reg = rng.normal(0, sd, n)
    a = int(n * 0.30)
    reg[a:a + int(n * frac)] += effect * sd
    return {"no_edge": regime_report(flat, n_null=n_null, seed=seed),
            "regime_edge": regime_report(reg, n_null=n_null, seed=seed)}
