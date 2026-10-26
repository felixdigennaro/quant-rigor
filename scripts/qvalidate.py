#!/usr/bin/env python3
"""qvalidate — overfitting detection and honest validation. STDLIB ONLY (no numpy/scipy).

The tools that separate a backtest from a discovery: walk-forward and purged CV splits,
Deflated Sharpe, PBO, HAC/Newey-West errors, minimum track record length, FDR control.

Reference notes: references/05-validation.md (formulas 31-37); RAS, formula 68, in
references/08-microstructure-costs.md.

    import sys; sys.path.insert(0,'<skill>/scripts')
    from qvalidate import deflated_sharpe, pbo, newey_west_se, purged_kfold_splits
"""
import random
import math
import statistics as st
from itertools import combinations

EULER_GAMMA = 0.5772156649015329

# ---------------------------------------------------------------- normal helpers


def norm_cdf(z):
    """Φ(z) via erf — exact to double precision."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def norm_ppf(p):
    """Inverse normal CDF (Acklam's rational approximation, |error| < 1.15e-9)."""
    if not 0.0 < p < 1.0:
        if p == 0:
            return float('-inf')
        if p == 1:
            return float('inf')
        return float('nan')
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    elif p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    else:
        q = p - 0.5
        r = q * q
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    # one Halley refinement step
    e = norm_cdf(x) - p
    u = e * math.sqrt(2 * math.pi) * math.exp(x * x / 2)
    return x - u / (1 + x * u / 2)

# ---------------------------------------------------------------- 35. HAC / Newey-West


def newey_west_se(x, lags=None):
    """HAC standard error of the MEAN (Bartlett kernel).

    A naive SE assumes independence. With autocorrelated or overlapping observations the naive
    SE is too small and every t-stat is inflated. Returns (se, lags_used, inflation_vs_naive).
    """
    n = len(x)
    if n < 3:
        return float('nan'), 0, float('nan')
    if lags is None:
        lags = max(1, int(round(4 * (n / 100.0) ** (2.0 / 9.0))))
    lags = min(lags, n - 2)
    m = st.mean(x)
    e = [v - m for v in x]
    g0 = sum(v * v for v in e) / n
    s = g0
    for k in range(1, lags + 1):
        gk = sum(e[i] * e[i + k] for i in range(n - k)) / n
        w = 1.0 - k / (lags + 1.0)
        s += 2.0 * w * gk
    if s <= 0 or g0 <= 0:                     # zero/degenerate variance: no usable SE
        return float('nan'), lags, float('nan')
    se_hac = math.sqrt(s / n)
    se_naive = math.sqrt(g0 * n / (n - 1) / n)
    return se_hac, lags, (se_hac / se_naive if se_naive else float('nan'))


def hac_t(x, lags=None):
    """t-statistic of the mean using HAC errors, plus the naive one for comparison."""
    se, L, infl = newey_west_se(x, lags)
    m = st.mean(x)
    naive_se = st.stdev(x) / math.sqrt(len(x)) if len(x) > 1 else float('nan')
    fin = lambda a, b: (a / b) if (b and math.isfinite(b) and b > 0) else float('nan')
    return {"mean": m, "t_hac": fin(m, se), "t_naive": fin(m, naive_se),
            "se_hac": se, "se_naive": naive_se, "lags": L, "se_inflation": infl}


def effective_sample_size(x, max_lag=None):
    """n_eff = n / (1 + 2 Σ ρ_k). Overlapping/serial observations count for less than 1."""
    n = len(x)
    if n < 4:
        return float(n)
    if max_lag is None:
        max_lag = max(1, int(round(4 * (n / 100.0) ** (2.0 / 9.0))))
    m = st.mean(x)
    var = sum((v - m) ** 2 for v in x)
    if var == 0:
        return float(n)
    s = 0.0
    for k in range(1, min(max_lag, n - 2) + 1):
        s += sum((x[i] - m) * (x[i + k] - m) for i in range(n - k)) / var
    d = 1 + 2 * s
    return min(float(n), n / d) if d > 0 else float(n)   # never report n_eff > n

# ---------------------------------------------------------------- 33. Deflated Sharpe


def expected_max_sharpe(n_trials, sd_of_trial_sharpes):
    """SR0 — the Sharpe you expect from the BEST of n_trials independent random strategies.
    Anything below this is what searching produces on noise alone."""
    N = max(2, int(n_trials))
    a = (1 - EULER_GAMMA) * norm_ppf(1 - 1.0 / N)
    b = EULER_GAMMA * norm_ppf(1 - 1.0 / (N * math.e))
    return sd_of_trial_sharpes * (a + b)


def deflated_sharpe(sr_observed, n_obs, n_trials, sd_of_trial_sharpes,
                    skew=0.0, kurtosis=3.0):
    """DSR (Bailey & López de Prado): probability the observed Sharpe is real, after correcting
    for (a) the number of configurations tried, (b) non-normality, (c) sample length.

    sr_observed / sd_of_trial_sharpes must be on the SAME per-period basis (do not annualise one).
    Returns dict with sr0 (the hurdle) and dsr (a probability). DSR < 0.95 → not established.
    """
    if not all(math.isfinite(z) for z in (sr_observed, sd_of_trial_sharpes)):
        return {"sr0": float('nan'), "dsr": float('nan'),
                "verdict": "undefined — Sharpe not computable (zero variance?)"}
    sr0 = expected_max_sharpe(n_trials, sd_of_trial_sharpes)
    denom = 1.0 - skew * sr_observed + ((kurtosis - 1.0) / 4.0) * sr_observed ** 2
    if denom <= 0 or n_obs < 2:
        return {"sr0": sr0, "dsr": float('nan'), "note": "variance term non-positive"}
    z = (sr_observed - sr0) * math.sqrt(n_obs - 1) / math.sqrt(denom)
    return {"sr0": sr0, "z": z, "dsr": norm_cdf(z),
            "verdict": ("established (DSR>=0.95)" if norm_cdf(z) >= 0.95 else
                        "NOT established — consistent with search on noise")}

# ---------------------------------------------------------------- 36. Minimum track record


def min_track_record_length(sr_observed, skew=0.0, kurtosis=3.0, sr_benchmark=0.0, alpha=0.05):
    """MinTRL — how many periods you need before an observed Sharpe is credible.
    Same per-period basis for sr_observed and sr_benchmark."""
    if not math.isfinite(sr_observed) or sr_observed <= sr_benchmark:
        return float('inf')
    z = norm_ppf(1 - alpha)
    var_term = 1.0 - skew * sr_observed + ((kurtosis - 1.0) / 4.0) * sr_observed ** 2
    return 1.0 + max(0.0, var_term) * (z / (sr_observed - sr_benchmark)) ** 2

# ---------------------------------------------------------------- 31-32. splits


def walk_forward_splits(n, n_folds=5, anchored=True, min_train=None):
    """Rolling (anchored=False) or expanding/anchored (True) train/test splits.
    Returns [(train_idx_range, test_idx_range), …] as (start, stop) half-open tuples."""
    if n_folds < 1 or n < n_folds + 1:
        return []
    test_len = n // (n_folds + 1)
    if test_len < 1:
        return []
    if min_train is None:
        min_train = test_len
    out = []
    for k in range(n_folds):
        test_start = n - (n_folds - k) * test_len
        test_stop = test_start + test_len
        train_start = 0 if anchored else max(0, test_start - min_train)
        if test_start - train_start < 1:
            continue
        out.append(((train_start, test_start), (test_start, test_stop)))
    return out


def purged_kfold_splits(n, n_folds=5, label_span=1, embargo=0):
    """Purged K-fold (López de Prado). Removes from TRAIN any observation whose label window
    overlaps the test window, plus an embargo of `embargo` bars after the test block.

    Without purging, overlapping labels leak the test into the train and accuracy is inflated.
    Returns [(train_indices, test_indices), …] as explicit index lists.
    """
    if n_folds < 2 or n < n_folds:
        return []
    fold = n // n_folds
    out = []
    for k in range(n_folds):
        t0 = k * fold
        t1 = n if k == n_folds - 1 else (k + 1) * fold
        test = list(range(t0, t1))
        lo = t0 - (label_span - 1)          # a label starting here still touches the test
        hi = t1 + embargo
        train = [i for i in range(n) if (i + label_span - 1) < lo or i >= hi]
        out.append((train, test))
    return out

# ---------------------------------------------------------------- 34. PBO


def pbo(perf_matrix, n_splits=8, max_combos=200):
    """Probability of Backtest Overfitting via combinatorially symmetric cross-validation.

    perf_matrix: list of per-period performance lists, one per CONFIGURATION
                 (perf_matrix[c][t] = config c's result in period t).
    Splits periods into n_splits blocks; for every balanced train/test partition, picks the
    best config in-sample and records its out-of-sample rank. PBO = P(chosen config lands
    below the OOS median).

    PBO > 0.5 → the selection procedure is worse than random. Anything above ~0.2 means the
    'best' configuration is mostly a fitting artefact.
    """
    C = len(perf_matrix)
    if C < 2:
        return {"pbo": float('nan'), "note": "need >= 2 configurations"}
    T = min(len(r) for r in perf_matrix)
    n_splits = max(2, min(n_splits, T))
    if n_splits % 2:
        n_splits -= 1
    if n_splits < 2:
        return {"pbo": float('nan'), "note": "not enough periods"}
    edges = [round(i * T / n_splits) for i in range(n_splits + 1)]
    blocks = [list(range(edges[i], edges[i + 1])) for i in range(n_splits)]
    half = n_splits // 2
    combos = list(combinations(range(n_splits), half))
    if len(combos) > max_combos:                      # deterministic thinning, no RNG
        step = len(combos) / max_combos
        combos = [combos[int(i * step)] for i in range(max_combos)]
    logits, below = [], 0
    for tr in combos:
        te = [i for i in range(n_splits) if i not in tr]
        tr_idx = [i for b in tr for i in blocks[b]]
        te_idx = [i for b in te for i in blocks[b]]
        if not tr_idx or not te_idx:
            continue
        is_perf = [sum(perf_matrix[c][i] for i in tr_idx) for c in range(C)]
        best = max(range(C), key=lambda c: is_perf[c])
        oos = [sum(perf_matrix[c][i] for i in te_idx) for c in range(C)]
        order = sorted(range(C), key=lambda c: oos[c])
        rank = order.index(best) + 1                   # 1 = worst, C = best
        w = rank / (C + 1.0)
        if w <= 0.5:
            below += 1
        w = min(max(w, 1e-9), 1 - 1e-9)
        logits.append(math.log(w / (1 - w)))
    if not logits:
        return {"pbo": float('nan'), "note": "no valid partitions"}
    p = below / len(logits)
    return {"pbo": p, "n_partitions": len(logits),
            "median_logit": st.median(logits),
            "verdict": ("selection is worse than random" if p > 0.5 else
                        "acceptable (<0.2)" if p < 0.2 else "substantial overfitting risk")}

# ---------------------------------------------------------------- 37. FDR


def benjamini_hochberg(pvals, q=0.10):
    """BH step-up FDR control. Returns indices declared significant and the cutoff used.
    Far less brutal than Bonferroni when testing many variants."""
    m = len(pvals)
    if m == 0:
        return {"significant": [], "cutoff": float('nan')}
    order = sorted(range(m), key=lambda i: pvals[i])
    k = 0
    for rank, i in enumerate(order, start=1):
        if pvals[i] <= rank / m * q:
            k = rank
    sig = [order[j] for j in range(k)]
    return {"significant": sorted(sig), "cutoff": (k / m * q if k else 0.0),
            "n_tested": m, "q": q,
            "bonferroni_cutoff": q / m}


def bonferroni(pvals, alpha=0.05):
    return {"significant": [i for i, p in enumerate(pvals) if p <= alpha / max(1, len(pvals))],
            "cutoff": alpha / max(1, len(pvals))}

# ---------------------------------------------------------------- report


def report(returns, n_trials=1, sd_of_trial_sharpes=None, label="strategy"):
    """One-shot honest validation summary for a per-period return series."""
    n = len(returns)
    m = st.mean(returns)
    sd = st.stdev(returns) if n > 1 else 0.0
    sr = m / sd if (sd and math.isfinite(sd)) else float('nan')
    sk = 0.0
    ku = 3.0
    if n > 3 and sd:
        sk = (n / ((n - 1) * (n - 2))) * sum(((v - m) / sd) ** 3 for v in returns)
        ku = sum(((v - m) / sd) ** 4 for v in returns) / n
    h = hac_t(returns)
    L = [f"=== validation · {label} ===",
         f"  n {n} · n_eff {effective_sample_size(returns):.1f} · per-period SR {sr:+.3f} "
         f"(skew {sk:+.2f}, kurt {ku:.2f})",
         f"  t naive {h['t_naive']:+.2f} vs t HAC {h['t_hac']:+.2f} "
         f"(lags {h['lags']}, SE inflation {h['se_inflation']:.2f}×)"]
    mintrl = min_track_record_length(sr, sk, ku)
    L.append(f"  MinTRL for SR>0 at 95%: {mintrl:.0f} periods "
             f"({'REACHED' if n >= mintrl else 'NOT reached — ' + str(int(mintrl - n)) + ' more needed'})")
    if n_trials > 1:
        sdt = sd_of_trial_sharpes if sd_of_trial_sharpes else abs(sr) * 0.5
        d = deflated_sharpe(sr, n, n_trials, sdt, sk, ku)
        L.append(f"  best-of-{n_trials}: hurdle SR0 {d['sr0']:+.3f} · DSR {d['dsr']:.3f} → {d['verdict']}")
    return "\n".join(L)


# ---------------------------------------------------------------- 68. Rademacher / RAS
# Added after studying stefan-jansen/machine-learning-for-trading ch.16 §16.7.
# Reference: Paleologo, "Elements of Quantitative Investing" (2024) ch.8;
#            Mohri, Rostamizadeh & Talwalkar, "Foundations of Machine Learning" (2018).

def rademacher_complexity(series, reps=2000, seed=7, normalise=True):
    """Empirical Rademacher complexity of a FAMILY of strategy return series.

    R_hat = E_sigma[ max_i (1/T) * sum_t sigma_t * x_t^(i) ],  sigma_t drawn uniform in {-1,+1}

    WHY THIS BEATS THE DEFLATED SHARPE FOR A REAL SEARCH: DSR corrects for the NUMBER of
    trials assuming they are independent. A real parameter sweep produces highly CORRELATED
    variants (lookback 11 vs 12, threshold 0.95 vs 1.00). Over a correlated family the
    supremum is much smaller than over an independent one, so DSR over-penalises and can
    reject an edge that survives. Rademacher complexity measures the family you ACTUALLY
    searched, whatever its correlation structure, and is distribution-free.

    `series`: list of equal-length return series (one per configuration tried).
    With normalise=True each series is scaled to unit sd, so R_hat is in per-period
    Sharpe units — multiply by sqrt(periods_per_year) to annualise.

    FAILURE MODE: it bounds the OPTIMISM OF THE SEARCH, not the truth of any single
    strategy. A small R_hat does not make a bad strategy good; it only means your
    selection procedure did not buy you much by looking.
    """
    if not series:
        return float('nan')
    T = min(len(s) for s in series)
    if T < 4:
        return float('nan')
    X = []
    for s in series:
        s = list(s[:T])
        if normalise:
            m = st.mean(s)
            sd = st.stdev(s) if len(s) > 1 else 0.0
            s = [(v - m) / sd for v in s] if sd > 0 else [0.0] * T
        X.append(s)
    rnd = random.Random(seed)
    total = 0.0
    for _ in range(reps):
        sg = [1.0 if rnd.random() < 0.5 else -1.0 for _ in range(T)]
        best = max(sum(sg[t] * x[t] for t in range(T)) / T for x in X)
        total += best
    return total / reps


def massart_bound(n_strategies, n_obs):
    """Massart's finite-class upper bound on Rademacher complexity: sqrt(2*ln(N)/T).

    The worst case, attained only when the N strategies are mutually uncorrelated. Compare
    it with the empirical `rademacher_complexity`: the ratio tells you how much your search
    was really worth. A ratio far below 1 means your variants were near-duplicates and the
    effective number of independent trials is far smaller than N.
    """
    if n_strategies < 1 or n_obs < 1:
        return float('nan')
    return math.sqrt(2.0 * math.log(max(n_strategies, 1)) / n_obs)


def effective_trials(n_strategies, n_obs, r_hat):
    """The number of INDEPENDENT trials your correlated search was equivalent to.

    Inverts Massart: N_eff = exp(T * r_hat^2 / 2). Report this instead of the raw count of
    configurations — quoting "we tried 217 variants" overstates the search when 200 of them
    were neighbours on the same grid.
    """
    if not math.isfinite(r_hat) or r_hat <= 0 or n_obs < 1:
        return float('nan')
    return math.exp(min(700.0, n_obs * r_hat * r_hat / 2.0))


def ras_sharpe_adjustment(observed_sharpe, series, periods_per_year=252, reps=2000,
                          seed=7, delta=0.05):
    """Rademacher Anti-Serum: a distribution-free floor on the true Sharpe after a search.

        theta_true >= theta_observed - 2*R_hat - sqrt(log(1/delta) / (2T))

    All terms are converted to ANNUALISED Sharpe. If the floor is <= 0 the search has not
    established anything, however good the winner looked.

    Report this ALONGSIDE the Deflated Sharpe, never instead of it: they answer slightly
    different questions and disagreeing is informative.
    """
    T = min(len(s) for s in series) if series else 0
    if T < 4:
        return {"error": "need at least 4 observations"}
    r_hat = rademacher_complexity(series, reps=reps, seed=seed)
    ann = math.sqrt(periods_per_year)
    est_err = math.sqrt(math.log(1.0 / delta) / (2.0 * T))
    floor = observed_sharpe - 2.0 * r_hat * ann - est_err * ann
    mb = massart_bound(len(series), T)
    return {"observed_sharpe_ann": observed_sharpe,
            "rademacher_per_period": r_hat, "rademacher_ann": r_hat * ann,
            "massart_bound_per_period": mb, "massart_ann": mb * ann,
            "correlation_discount": (r_hat / mb) if mb > 0 else float('nan'),
            "effective_trials": effective_trials(len(series), T, r_hat),
            "n_configurations": len(series), "n_obs": T,
            "estimation_error_ann": est_err * ann,
            "sharpe_floor_ann": floor,
            "verdict": ("established: floor above zero after the search penalty" if floor > 0
                        else "NOT established — the search penalty consumes the entire edge")}


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", nargs="+", type=float, required=True)
    ap.add_argument("--trials", type=int, default=1)
    ap.add_argument("--sd-trials", type=float, default=None)
    a = ap.parse_args()
    print(report(a.series, a.trials, a.sd_trials))
