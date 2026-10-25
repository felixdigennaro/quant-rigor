#!/usr/bin/env python3
"""qconfluence — conditional edge validation, for "size up when the setup is better".

WHY THIS EXISTS.
The idea is standard and sound: some setups are better than others, so trade them bigger.
The *procedure* people use to find which ones is not sound, and it fails in a way that is
invisible from the inside — it returns a confident positive when there is nothing there.

Measured, 2026-08-03, real case: a mean-reversion strategy on 16 years of NQ, ~22 context
variables, all singles and pairs = 968 candidate rules. Search for the best, report its t.

    best rule on real P&L .......... t = +3.38
    best rule on SHUFFLED P&L ...... t = +3.02  (median of 200 replications)
    95th percentile under the null .. t = +3.97

The search itself manufactures t around +3. The "discovery" was inside its own noise band.
Applied to gold with the rule frozen: t = +0.15 and t = -0.86 (sign flipped). Nothing there.

So the rule is: **the t of a searched rule is meaningless until compared to the t that the
SAME SEARCH produces on data with no edge.** Not compared to 2.0. Not compared to a textbook
threshold. Compared to its own null.

THE FIVE GATES, in order. Each can kill the finding, and killing early saves the rest.

  1. DOSE-RESPONSE, not a binary filter.
     Score each period by how many confluences align, then require the mean P&L to rise
     MONOTONICALLY across buckets. A top-tercile-vs-rest split is the weakest possible
     evidence and the easiest for noise to fake; a monotone gradient across 5 buckets is
     much harder. This is also the only form that answers the actual question — if you are
     going to size by confluence, you need the slope, not one bucket's mean.

  2. NULL CALIBRATION OF THE SEARCH.
     Re-run the entire search on shuffled P&L, many times. Your t is evidence only if it
     beats that distribution's 95th percentile.

  3. OUT-OF-SAMPLE.
     Select on the first part, apply frozen to the last part. No re-tuning, no second pick.

  4. INDEPENDENT REPLICATION.
     Same frozen rule, different instrument. This is what killed the real case above.

  5. CONCENTRATION.
     Does the effect live in 1-3 periods? Jackknife it. A confluence that only works on
     three days is a story about three days.

ONLY IF ALL FIVE PASS does sizing follow — and then the weight comes from the fitted slope,
capped, never from the ratio of bucket means (which is unbounded and noisy at the extremes).

THE COST OF BEING WRONG, which decides how careful to be:
`sizing_cost()` computes what confluence-sizing does to the Sharpe when the confluence is
pure noise. Sizing 1x/3x on a coin flip multiplies Sharpe by mean(w)/sqrt(mean(w^2)) = 0.89:
an 11% loss for nothing. The upside if real is larger, which is exactly why the gates matter.
"""
import math

import numpy as np

__all__ = ["dose_response", "search_null", "confluence_gate", "sizing_from_slope",
           "sizing_cost", "jackknife_sign"]


def _t(x):
    x = np.asarray(x, float)
    if len(x) < 3:
        return float("nan")
    s = x.std(ddof=1)
    return float(x.mean() / (s / math.sqrt(len(x)))) if s > 0 else float("nan")


# ---------------------------------------------------------------- gate 1: dose-response


def dose_response(score, pnl, n_buckets=5, min_per_bucket=20):
    """Does P&L rise MONOTONICALLY with the confluence score?

    `score`: per-period confluence strength (count of aligned conditions, or continuous).
    Returns bucket means, the OLS slope with its SE and t, and a monotonicity verdict.

    Why monotonicity and not just "top bucket beats bottom": a spurious rule typically
    produces one good bucket, not an ordered gradient. Requiring order is a much stronger
    test at zero extra cost, and it is the shape you need if you intend to size by it.
    """
    score = np.asarray(score, float)
    pnl = np.asarray(pnl, float)
    m = np.isfinite(score) & np.isfinite(pnl)
    score, pnl = score[m], pnl[m]
    if len(score) < 3 * min_per_bucket:
        return {"ok": False, "note": f"servono almeno {3*min_per_bucket} osservazioni, "
                                     f"ce ne sono {len(score)}"}
    # UN PUNTEGGIO DI CONFLUENZA È QUASI SEMPRE DISCRETO ("quante condizioni sono allineate"),
    # e su valori ripetuti np.quantile produce estremi duplicati → bucket VUOTI → medie nan →
    # il test di monotonia falliva SEMPRE, anche su un effetto vero iniettato apposta.
    # Il controllo positivo l'ha scoperto: un cancello che boccia anche il vero non è un cancello.
    uniq = np.unique(score)
    if len(uniq) <= max(n_buckets, 8):
        groups = [(u, score == u) for u in uniq]          # discreto: un gruppo per valore
    else:
        qs = np.quantile(score, np.linspace(0, 1, n_buckets + 1))
        qs[-1] += 1e-12
        groups = [(qs[i], (score >= qs[i]) & (score < qs[i + 1])) for i in range(n_buckets)]
    means, ns, ses, edges = [], [], [], []
    for edge, sel in groups:
        v = pnl[sel]
        if len(v) < min_per_bucket:                        # bucket troppo piccolo: si SCARTA,
            continue                                       # non si tiene con una media nan
        edges.append(float(edge))
        ns.append(int(sel.sum()))
        means.append(float(v.mean()))
        ses.append(float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 2 else float("nan"))
    if len(means) < 3:
        return {"ok": False, "note": f"solo {len(means)} bucket con almeno {min_per_bucket} "
                                     f"osservazioni: la monotonia non è verificabile"}
    mm = np.array(means, float)
    mono_up = bool(np.all(np.diff(mm) >= 0))
    mono_dn = bool(np.all(np.diff(mm) <= 0))
    # Spearman: monotonicity of the underlying relation, not just of the bucket means
    rs = _spearman(score, pnl)
    # OLS slope of pnl on standardised score
    z = (score - score.mean()) / (score.std(ddof=1) if score.std(ddof=1) > 0 else 1.0)
    slope = float(np.cov(z, pnl)[0, 1] / z.var(ddof=1)) if z.var(ddof=1) > 0 else float("nan")
    resid = pnl - (pnl.mean() + slope * z)
    se = float(resid.std(ddof=1) / math.sqrt(len(z))) if len(z) > 2 else float("nan")
    return {"ok": True, "buckets": means, "n": ns, "se": ses, "edges": edges,
            "monotone_up": mono_up, "monotone_down": mono_dn,
            "spearman": rs, "slope_per_sd": slope, "slope_se": se,
            "slope_t": slope / se if se and se == se and se > 0 else float("nan"),
            "spread_top_bottom": float(mm[-1] - mm[0])}


def _spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = math.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


# ---------------------------------------------------------------- gate 2: null of the search


def search_null(masks, pnl, n_rep=200, min_n=30, seed=0, train=None):
    """What t does THIS SEARCH produce when there is no edge at all?

    `masks`: (n_rules, n_periods) boolean array — every candidate rule you tried.
    Shuffling `pnl` destroys any association with the rules while preserving the return
    distribution exactly. The best-of-N t you get back is the bar your real t must clear.

    This is the single highest-yield function here. Without it, a searched t is a number
    with no interpretation: 968 rules put the median best-t at +3.02 in the measured case.
    """
    masks = np.asarray(masks, bool)
    pnl = np.asarray(pnl, float)
    sl = train if train is not None else slice(0, len(pnl))
    sub = masks[:, sl]
    cnt = sub.sum(1)
    ok = cnt >= min_n
    if not ok.any():
        return {"ok": False, "note": "nessuna regola raggiunge il minimo di osservazioni"}
    s = sub.astype(float)
    rng = np.random.default_rng(seed)
    best = []
    for _ in range(n_rep):
        p = rng.permutation(pnl)[sl]
        mean = (s @ p) / np.maximum(cnt, 1)
        sq = (s @ (p ** 2)) / np.maximum(cnt, 1)
        var = np.maximum(sq - mean ** 2, 0) * cnt / np.maximum(cnt - 1, 1)
        t = mean / np.maximum(np.sqrt(var / np.maximum(cnt, 1)), 1e-12)
        best.append(float(np.max(np.where(ok, t, -np.inf))))
    b = np.array(best)
    return {"ok": True, "n_rules": int(masks.shape[0]), "n_rep": n_rep,
            "median": float(np.median(b)), "q95": float(np.percentile(b, 95)),
            "max": float(b.max()), "dist": b}


# ---------------------------------------------------------------- gate 5: concentration


def jackknife_sign(pnl, blocks=None):
    """Drop one block at a time: does the sign of the mean survive every deletion?

    Also reports what fraction of the total comes from the best 3 observations — an effect
    that lives in three days is a story about three days, whatever its t says.
    """
    p = np.asarray(pnl, float)
    p = p[np.isfinite(p)]
    if len(p) < 10:
        return {"ok": False, "note": "troppe poche osservazioni"}
    n_b = blocks if blocks else min(20, len(p))
    idx = np.array_split(np.arange(len(p)), n_b)
    means = [float(np.delete(p, i).mean()) for i in idx]
    tot = p.sum()
    top3 = float(np.sort(p)[-3:].sum() / tot) if tot != 0 else float("nan")
    return {"ok": True, "sign_survives_all": bool(all((m > 0) == (p.mean() > 0) for m in means)),
            "worst_mean": float(min(means)), "best_mean": float(max(means)),
            "top3_share": top3, "n_blocks": n_b}


# ---------------------------------------------------------------- sizing


def sizing_from_slope(score, slope_per_sd, base_mean, w_min=0.5, w_max=2.0):
    """Weights from the FITTED SLOPE, capped. Never from ratios of bucket means.

    A ratio of bucket means is unbounded and explodes when a bucket mean is near zero,
    which is exactly where the estimate is least reliable. The slope is a single number
    with a standard error, and capping it bounds the damage when it is wrong.
    """
    z = (np.asarray(score, float) - np.mean(score)) / (np.std(score, ddof=1) or 1.0)
    if base_mean == 0:
        return np.ones_like(z)
    w = 1.0 + (slope_per_sd / abs(base_mean)) * z
    return np.clip(w, w_min, w_max)


def sizing_cost(weights):
    """If the confluence is PURE NOISE, what does sizing by it cost?

    Expected return per unit of capital is unchanged, but variance rises, so
        Sharpe_after / Sharpe_before = mean(w) / sqrt(mean(w^2)) <= 1
    Equality only if every weight is identical. 1x/3x on a coin flip → 0.894, an 11% loss.
    Run this BEFORE deciding how aggressive the sizing spread should be.
    """
    w = np.asarray(weights, float)
    w = w[np.isfinite(w)]
    if len(w) == 0:
        return float("nan")
    m2 = float(np.mean(w ** 2))
    return float(np.mean(w) / math.sqrt(m2)) if m2 > 0 else float("nan")


# ---------------------------------------------------------------- the whole gate


def confluence_gate(score, pnl, masks=None, split=0.6, n_rep=200, n_buckets=5,
                    replication=None, min_n=30, seed=0):
    """The five gates in order. Returns a verdict that DERIVES from them, not a summary.

    `score`       per-period confluence strength (this is what you would size on)
    `pnl`         per-period P&L of the base strategy
    `masks`       every rule the search considered — needed for gate 2. Omit only if the
                  confluence was pre-registered and NOT searched; say so explicitly.
    `replication` (score, pnl) from a DIFFERENT instrument, rule frozen — gate 4
    """
    score = np.asarray(score, float)
    pnl = np.asarray(pnl, float)
    n = len(pnl)
    cut = int(n * split)
    res = {"n": n, "gates": {}}

    dr = dose_response(score, pnl, n_buckets=n_buckets)
    res["dose_response"] = dr
    res["gates"]["1_dose_response"] = bool(dr.get("ok") and
                                           (dr.get("monotone_up") or dr.get("monotone_down")))

    if masks is not None:
        nl = search_null(masks, pnl, n_rep=n_rep, min_n=min_n, seed=seed,
                         train=slice(0, cut))
        res["search_null"] = {k: v for k, v in nl.items() if k != "dist"}
        real_t = dr.get("slope_t", float("nan"))
        res["gates"]["2_beats_own_null"] = bool(nl.get("ok") and real_t == real_t
                                                and real_t > nl["q95"])
    else:
        res["gates"]["2_beats_own_null"] = None
        res["search_null"] = {"note": "nessuna maschera fornita — vale SOLO se la confluenza "
                                      "era pre-registrata e non cercata. Dichiararlo."}

    dr_tr = dose_response(score[:cut], pnl[:cut], n_buckets=n_buckets)
    dr_te = dose_response(score[cut:], pnl[cut:], n_buckets=n_buckets)
    res["train"], res["test"] = dr_tr, dr_te
    st, se_ = dr_tr.get("slope_per_sd", float("nan")), dr_te.get("slope_per_sd", float("nan"))
    res["gates"]["3_out_of_sample"] = bool(st == st and se_ == se_ and (st > 0) == (se_ > 0)
                                           and abs(se_) > 0.25 * abs(st))

    if replication is not None:
        rs, rp = replication
        dr_r = dose_response(rs, rp, n_buckets=n_buckets)
        res["replication"] = dr_r
        sr = dr_r.get("slope_per_sd", float("nan"))
        res["gates"]["4_replicates"] = bool(sr == sr and st == st and (sr > 0) == (st > 0)
                                            and abs(dr_r.get("slope_t", 0)) > 2.0)
    else:
        res["gates"]["4_replicates"] = None

    top = score >= np.quantile(score, 1 - 1.0 / n_buckets)
    jk = jackknife_sign(pnl[top])
    res["concentration"] = jk
    res["gates"]["5_not_concentrated"] = bool(jk.get("ok") and jk["sign_survives_all"]
                                              and jk["top3_share"] < 0.5)

    passed = [v for v in res["gates"].values() if v is True]
    failed = [k for k, v in res["gates"].items() if v is False]
    untested = [k for k, v in res["gates"].items() if v is None]
    if failed:
        res["verdict"] = "NON USARE PER DIMENSIONARE — cancelli falliti: " + ", ".join(failed)
    elif untested:
        res["verdict"] = ("NON DECIDIBILE — cancelli non eseguiti: " + ", ".join(untested) +
                          ". Un cancello non eseguito non è un cancello superato.")
    else:
        res["verdict"] = "UTILIZZABILE PER DIMENSIONARE (con peso limitato)"
    res["n_passed"] = len(passed)
    if dr.get("ok"):
        w = sizing_from_slope(score, dr["slope_per_sd"], abs(np.mean(pnl)) or 1.0)
        res["sizing_weights_range"] = (float(w.min()), float(w.max()))
        res["sharpe_cost_if_noise"] = sizing_cost(w)
    return res
