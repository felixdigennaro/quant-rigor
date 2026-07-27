#!/usr/bin/env python3
"""qmontecarlo — null distributions by resampling and simulation. STDLIB ONLY (no numpy/scipy).

Answers the question a t-statistic cannot: how often does a process with NO edge produce a result
as good as mine? Four families, each attacking a different null:

  permutation      — destroy the timing, keep the distribution. Tests whether the ORDER matters.
  sign flip        — keep the magnitudes, randomise direction. Tests whether the MEAN is real.
  surrogate paths  — simulate under an explicit no-edge process (block bootstrap, GBM, OU).
  Reality Check /  — the formal multiple-comparison test when many strategies are screened
  SPA                against a benchmark (White 2000, Hansen 2005).

Plus `stress_moments`: does the verdict survive the estimation error in its own inputs?

All routines are deterministic given `seed`. Reference notes: references/07-montecarlo.md (45-49).
"""
import math
import random
import statistics as st

# ---------------------------------------------------------------- helpers


def _mean(x):
    return st.mean(x) if x else float('nan')


def _sharpe(x):
    if len(x) < 2:
        return float('nan')
    sd = st.stdev(x)
    return (st.mean(x) / sd) if sd > 0 else float('nan')


def _cum(x):
    out, a = [], 0.0
    for v in x:
        a += v
        out.append(a)
    return out


def _max_dd(x):
    peak, worst = 0.0, 0.0
    for c in _cum(x):
        peak = max(peak, c)
        worst = min(worst, c - peak)
    return worst


STATS = {"total": sum, "mean": _mean, "sharpe": _sharpe, "max_dd": _max_dd}


def _pval(observed, null_samples, tail="greater"):
    """Empirical p-value with the +1 correction (never returns exactly 0, which would be a lie)."""
    n = len(null_samples)
    if n == 0:
        return float('nan')
    if tail == "greater":
        k = sum(1 for v in null_samples if v >= observed)
    elif tail == "less":
        k = sum(1 for v in null_samples if v <= observed)
    else:
        k = sum(1 for v in null_samples if abs(v) >= abs(observed))
    return (k + 1) / (n + 1)


def _summary(observed, null_samples, tail, label, reps):
    ns = sorted(null_samples)
    q = lambda p: ns[min(len(ns) - 1, int(p * len(ns)))] if ns else float('nan')
    p = _pval(observed, null_samples, tail)
    return {"test": label, "observed": observed, "p_value": p, "reps": reps,
            "null_mean": _mean(null_samples), "null_p05": q(0.05), "null_p50": q(0.50),
            "null_p95": q(0.95), "null_p99": q(0.99),
            "verdict": ("distinguishable from the null (p<0.01)" if p < 0.01 else
                        "distinguishable from the null (p<0.05)" if p < 0.05 else
                        "borderline (0.05<=p<0.10)" if p < 0.10 else
                        "NOT distinguishable from the null")}

# ---------------------------------------------------------------- 45. permutation / sign flip


def permutation_test(returns, stat="sharpe", reps=5000, seed=7, block=1):
    """Shuffle the ORDER of the periods. Preserves the marginal distribution exactly and destroys
    any timing/sequencing structure.

    Use it when the claim depends on sequence: trend following, autocorrelation capture,
    stop placement, compounding. Note that for a plain sum or mean, reordering changes nothing —
    so `total`/`mean` are meaningless here and only `sharpe`/`max_dd` are informative.
    `block > 1` preserves local autocorrelation (a stricter null).
    """
    if stat in ("total", "mean"):
        return {"test": f"permutation({stat})", "p_value": float('nan'),
                "verdict": "not applicable — reordering cannot change a sum or a mean; "
                           "use sign_flip_test for the mean, or stat='sharpe'/'max_dd'"}
    f = STATS[stat]
    rnd = random.Random(seed)
    obs = f(returns)
    n = len(returns)
    null = []
    for _ in range(reps):
        if block <= 1:
            s = returns[:]
            rnd.shuffle(s)
        else:
            nb = max(1, n // block)
            blocks = [returns[i * block:(i + 1) * block] for i in range(nb + 1)]
            blocks = [b for b in blocks if b]
            rnd.shuffle(blocks)
            s = [v for b in blocks for v in b][:n]
        null.append(f(s))
    tail = "less" if stat == "max_dd" else "greater"
    return _summary(obs, null, tail, f"permutation({stat}, block={block})", reps)


def sign_flip_test(returns, stat="mean", reps=5000, seed=7):
    """Randomise the SIGN of each period, keeping magnitudes. The null is 'symmetric noise with
    this volatility'. This is the correct permutation-style test for the mean/total, and it is
    stricter than a t-test because it makes no distributional assumption."""
    f = STATS[stat]
    rnd = random.Random(seed)
    obs = f(returns)
    null = []
    for _ in range(reps):
        s = [v if rnd.random() < 0.5 else -v for v in returns]
        null.append(f(s))
    tail = "less" if stat == "max_dd" else "greater"
    return _summary(obs, null, tail, f"sign_flip({stat})", reps)

# ---------------------------------------------------------------- 46. surrogate paths


def surrogate_paths(returns, model="block_bootstrap", stat="sharpe", reps=5000, seed=7,
                    block=5, drift_zero=True):
    """Simulate paths under an EXPLICIT no-edge process calibrated to the observed moments, then
    ask how often the null reproduces the observed statistic.

    model:
      'block_bootstrap' — resample blocks of the actual returns (keeps fat tails and local
                          autocorrelation; the most honest default)
      'iid_bootstrap'   — resample single periods (destroys autocorrelation)
      'gbm'             — Normal returns with the observed sd (and zero drift if drift_zero)
      'ou'              — mean-reverting, theta estimated by lag-1 autoregression on the level

    drift_zero=True removes the observed mean, which is what makes it a NO-EDGE null.
    """
    n = len(returns)
    if n < 3:
        return {"test": f"surrogate({model})", "p_value": float('nan'),
                "verdict": "sample too small"}
    f = STATS[stat]
    obs = f(returns)
    rnd = random.Random(seed)
    mu, sd = st.mean(returns), st.stdev(returns)
    centred = [v - mu for v in returns] if drift_zero else returns[:]
    # OU parameters from a lag-1 regression on the cumulative level
    lvl = _cum(returns)
    theta = 0.0
    if model == "ou" and n > 3:
        lm = st.mean(lvl)
        num = sum((lvl[i] - lm) * (lvl[i + 1] - lvl[i]) for i in range(n - 1))
        den = sum((lvl[i] - lm) ** 2 for i in range(n - 1))
        theta = max(0.0, min(1.0, -num / den if den > 0 else 0.0))
    null = []
    for _ in range(reps):
        if model == "iid_bootstrap":
            s = [centred[rnd.randrange(n)] for _ in range(n)]
        elif model == "block_bootstrap":
            s = []
            while len(s) < n:
                i = rnd.randrange(n)
                s.extend(centred[i:i + block] if i + block <= n
                         else centred[i:] + centred[:(i + block) % n])
            s = s[:n]
        elif model == "gbm":
            s = [rnd.gauss(0.0, sd) for _ in range(n)]
        elif model == "ou":
            s, x = [], 0.0
            for _ in range(n):
                shock = rnd.gauss(0.0, sd)
                step = -theta * x + shock
                x += step
                s.append(step)
        else:
            raise ValueError(f"unknown model {model}")
        null.append(f(s))
    tail = "less" if stat == "max_dd" else "greater"
    out = _summary(obs, null, tail, f"surrogate({model}, {stat})", reps)
    out["calibration"] = {"observed_mean": mu, "observed_sd": sd,
                          "drift_removed": drift_zero, "ou_theta": theta if model == "ou" else None}
    return out

# ---------------------------------------------------------------- 47-48. Reality Check / SPA


def reality_check(loss_diffs, reps=2000, seed=7, block=5):
    """White's (2000) Reality Check for data snooping.

    loss_diffs: list of series, one per strategy, of the strategy's performance MINUS the
    benchmark's, period by period (positive = strategy beat the benchmark that period).

    Tests H0: the BEST of the k strategies does not beat the benchmark. This is the correct test
    when many candidates were screened — a per-strategy p-value is not.
    Returns the p-value for the best strategy and the per-strategy nominal means.
    """
    k = len(loss_diffs)
    if k == 0:
        return {"test": "reality_check", "p_value": float('nan'), "verdict": "no strategies"}
    T = min(len(s) for s in loss_diffs)
    if T < 4:
        return {"test": "reality_check", "p_value": float('nan'), "verdict": "sample too small"}
    means = [st.mean(s[:T]) for s in loss_diffs]
    V = max(math.sqrt(T) * m for m in means)
    best = max(range(k), key=lambda i: means[i])
    rnd = random.Random(seed)
    null = []
    for _ in range(reps):
        idx = []
        while len(idx) < T:                             # stationary/block bootstrap of the periods
            i = rnd.randrange(T)
            idx.extend(range(i, min(T, i + block)))
        idx = idx[:T]
        vstar = -float('inf')
        for j in range(k):
            mstar = st.mean([loss_diffs[j][t] for t in idx])
            vstar = max(vstar, math.sqrt(T) * (mstar - means[j]))   # re-centred under H0
        null.append(vstar)
    p = _pval(V, null, "greater")
    return {"test": "White reality_check", "n_strategies": k, "periods": T,
            "best_strategy_index": best, "best_mean_excess": means[best],
            "V_statistic": V, "p_value": p, "reps": reps,
            "per_strategy_mean_excess": means,
            "verdict": ("the best strategy beats the benchmark after correcting for the search "
                        f"(p={p:.3f})" if p < 0.05 else
                        "NO strategy survives the data-snooping correction "
                        f"(p={p:.3f}) — the winner is consistent with screening noise")}


def spa_test(loss_diffs, reps=2000, seed=7, block=5, threshold=None):
    """Hansen's (2005) Superior Predictive Ability test — a studentised, re-centred refinement of
    the Reality Check that is far less sensitive to including poor candidates in the comparison set.

    Same input as reality_check. Returns the SPA p-value (consistent variant).
    """
    k = len(loss_diffs)
    if k == 0:
        return {"test": "spa", "p_value": float('nan'), "verdict": "no strategies"}
    T = min(len(s) for s in loss_diffs)
    if T < 4:
        return {"test": "spa", "p_value": float('nan'), "verdict": "sample too small"}
    means = [st.mean(s[:T]) for s in loss_diffs]
    sds = [max(1e-12, st.stdev(s[:T])) if T > 1 else 1e-12 for s in loss_diffs]
    Tstat = max(0.0, max(math.sqrt(T) * means[j] / sds[j] for j in range(k)))
    if threshold is None:                                # Hansen's recommended rate
        threshold = -math.sqrt(2.0 * math.log(math.log(max(3, T))))
    keep = [j for j in range(k)
            if math.sqrt(T) * means[j] / sds[j] >= threshold]        # drop hopeless candidates
    rnd = random.Random(seed)
    null = []
    for _ in range(reps):
        idx = []
        while len(idx) < T:
            i = rnd.randrange(T)
            idx.extend(range(i, min(T, i + block)))
        idx = idx[:T]
        z = 0.0
        for j in keep:
            mstar = st.mean([loss_diffs[j][t] for t in idx])
            z = max(z, math.sqrt(T) * (mstar - means[j]) / sds[j])
        null.append(z)
    p = _pval(Tstat, null, "greater")
    return {"test": "Hansen SPA", "n_strategies": k, "n_retained": len(keep), "periods": T,
            "T_statistic": Tstat, "p_value": p, "reps": reps, "threshold": threshold,
            "verdict": (f"superior predictive ability established (p={p:.3f})" if p < 0.05 else
                        f"NO superior predictive ability (p={p:.3f})")}

# ---------------------------------------------------------------- 49. moment stress


def stress_moments(returns, cost_per_period=0.0, scenarios=None, stat="sharpe"):
    """Does the conclusion survive the estimation error in its own inputs?

    Applies explicit, adversarial perturbations and reports the statistic and t under each.
    Default scenarios: volatility understated by 30%, costs doubled, the best period removed,
    the two best removed, mean shrunk 25% (parameter uncertainty), and a fat-tail shock added.
    """
    n = len(returns)
    if n < 3:
        return {"test": "stress_moments", "rows": [], "verdict": "sample too small"}
    f = STATS[stat]
    base = f(returns)
    srt = sorted(range(n), key=lambda i: returns[i])
    worst_add = min(returns) * 2.0

    def _t(x):
        sd = st.stdev(x) if len(x) > 1 else 0.0
        return (st.mean(x) / (sd / math.sqrt(len(x)))) if sd > 0 else float('nan')

    sc = scenarios or [
        ("baseline", returns[:]),
        ("volatility +30% (understated risk)",
         [st.mean(returns) + 1.3 * (v - st.mean(returns)) for v in returns]),
        ("costs doubled", [v - cost_per_period for v in returns]),
        ("best period removed", [returns[i] for i in range(n) if i != srt[-1]]),
        ("two best periods removed", [returns[i] for i in range(n) if i not in (srt[-1], srt[-2])]),
        ("mean shrunk 25% (parameter error)",
         [0.75 * st.mean(returns) + (v - st.mean(returns)) for v in returns]),
        ("one tail event added (2x the worst period)", returns + [worst_add]),
    ]
    rows = []
    for name, s in sc:
        rows.append({"scenario": name, "n": len(s), stat: f(s), "t": _t(s),
                     "total": sum(s), "max_dd": _max_dd(s)})
    survives = [r for r in rows[1:] if (r[stat] > 0 if stat != "max_dd" else True)
                and abs(r["t"]) >= 1.6]
    return {"test": "stress_moments", "stat": stat, "baseline": base, "rows": rows,
            "n_scenarios_surviving_t1.6": len(survives), "n_scenarios": len(rows) - 1,
            "verdict": (f"conclusion survives {len(survives)}/{len(rows)-1} adversarial scenarios"
                        if survives else
                        "conclusion does NOT survive any adversarial scenario at |t|>=1.6")}

# ---------------------------------------------------------------- combined report


def report(returns, cost_per_period=0.0, label="strategy", reps=4000, seed=7):
    """Run the whole Monte Carlo battery and print the honest summary."""
    L = [f"=== monte carlo · {label} (n={len(returns)}, reps={reps}) ==="]
    for r in (sign_flip_test(returns, "mean", reps, seed),
              permutation_test(returns, "sharpe", reps, seed),
              permutation_test(returns, "max_dd", reps, seed),
              surrogate_paths(returns, "block_bootstrap", "sharpe", reps, seed),
              surrogate_paths(returns, "gbm", "sharpe", reps, seed)):
        if r.get("p_value") is None or (isinstance(r.get("p_value"), float)
                                        and math.isnan(r["p_value"])):
            L.append(f"  {r['test']:<38} n/a — {r['verdict']}")
            continue
        L.append(f"  {r['test']:<38} observed {r['observed']:+.4f} · null p50 "
                 f"{r['null_p50']:+.4f} p95 {r['null_p95']:+.4f} · p={r['p_value']:.4f}")
    sm = stress_moments(returns, cost_per_period)
    L.append(f"  {'stress across scenarios':<38} {sm['verdict']}")
    for row in sm["rows"]:
        L.append(f"      {row['scenario']:<44} sharpe {row['sharpe']:+.3f} · t {row['t']:+.2f}")
    return "\n".join(L)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", nargs="+", type=float, required=True)
    ap.add_argument("--cost", type=float, default=0.0)
    ap.add_argument("--reps", type=int, default=4000)
    a = ap.parse_args()
    print(report(a.series, a.cost, reps=a.reps))
