#!/usr/bin/env python3
"""qportfolio — portfolio risk, decomposition and weighting. STDLIB ONLY (no numpy).

Makes formulas 9-15 executable: covariance and correlation matrices, portfolio volatility,
marginal risk contributions, risk parity, eigendecomposition (Jacobi), PC1 concentration,
condition number, Ledoit-Wolf-style shrinkage, the Marchenko-Pastur noise edge, and
tail-conditional correlation.

Reference notes: references/02-linear-algebra.md (formulas 9-15).

Matrices are plain lists of lists; asset returns are {name: [r1, r2, …]} or a list of lists.
"""
import math
import statistics as st

# ---------------------------------------------------------------- 12. covariance / correlation


def _as_matrix(returns):
    if isinstance(returns, dict):
        names = list(returns.keys())
        rows = [list(returns[k]) for k in names]
    else:
        names = [f"a{i}" for i in range(len(returns))]
        rows = [list(r) for r in returns]
    n = min(len(r) for r in rows)
    return names, [r[:n] for r in rows], n


def cov_matrix(returns, ddof=1):
    """Sample covariance. Note: with T observations and N assets this needs T >> N to be usable;
    below ~10N it is badly conditioned and must not be inverted (see condition_number)."""
    names, rows, T = _as_matrix(returns)
    N = len(rows)
    means = [st.mean(r) for r in rows]
    d = max(1, T - ddof)
    C = [[sum((rows[i][t] - means[i]) * (rows[j][t] - means[j]) for t in range(T)) / d
          for j in range(N)] for i in range(N)]
    return names, C


def corr_matrix(returns):
    names, C = cov_matrix(returns)
    N = len(C)
    sd = [math.sqrt(C[i][i]) if C[i][i] > 0 else 0.0 for i in range(N)]
    R = [[(1.0 if i == j else
           (C[i][j] / (sd[i] * sd[j]) if sd[i] > 0 and sd[j] > 0 else 0.0))
          for j in range(N)] for i in range(N)]
    return names, R


def tail_conditional_corr(returns, quantile=0.10, reference=None):
    """Pearson correlation computed only over the worst `quantile` of periods for a reference
    series (default: the equal-weight portfolio).

    WARNING — this statistic is biased and routinely misread. Conditioning on a tail of the MEAN
    truncates the joint distribution, which mechanically *lowers* the sample correlation inside the
    subset even when true tail dependence is higher (Boyer-Gibson-Loretan). A falling number here is
    therefore NOT evidence that diversification improves in stress. Use `tail_dependence()` for the
    question people actually mean: how often do both series land in their own tails together.
    Kept because it is what most reports show — and so it can be compared against the correct one."""
    names, rows, T = _as_matrix(returns)
    ref = reference if reference is not None else [st.mean([rows[i][t] for i in range(len(rows))])
                                                  for t in range(T)]
    k = max(2, int(quantile * T))
    worst = sorted(range(T), key=lambda t: ref[t])[:k]
    sub = {names[i]: [rows[i][t] for t in worst] for i in range(len(rows))}
    return corr_matrix(sub)


def tail_dependence(returns, quantile=0.10):
    """Lower tail dependence: λ = P(X below its own q-quantile | Y below its own q-quantile),
    estimated non-parametrically for every pair. Free of the conditioning bias that distorts
    `tail_conditional_corr`.

    Benchmarks: independent series → λ ≈ quantile. λ >> quantile means the assets fail together.
    Returns (names, matrix of λ with 1.0 on the diagonal, mean off-diagonal λ, independence
    benchmark). Diversification should be sized on this, not on average correlation."""
    names, rows, T = _as_matrix(returns)
    N = len(rows)
    k = max(1, int(round(quantile * T)))
    tails = []
    for i in range(N):
        thr = sorted(rows[i])[k - 1]
        tails.append({t for t in range(T) if rows[i][t] <= thr})
    L = [[1.0] * N for _ in range(N)]
    off = []
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            both = len(tails[i] & tails[j])
            L[i][j] = both / len(tails[j]) if tails[j] else float('nan')
            if i < j:
                off.append(L[i][j])
    return names, L, (st.mean(off) if off else float('nan')), quantile


# ---------------------------------------------------------------- 9. portfolio algebra


def portfolio_vol(weights, cov):
    """σ_p = √(wᵀΣw)."""
    N = len(weights)
    v = sum(weights[i] * cov[i][j] * weights[j] for i in range(N) for j in range(N))
    return math.sqrt(max(0.0, v))


def marginal_risk_contribution(weights, cov):
    """Risk contribution of each position: wᵢ·(Σw)ᵢ/σ_p. Sums to σ_p.
    This answers 'who is actually carrying the risk' — usually not who you think."""
    N = len(weights)
    sp = portfolio_vol(weights, cov)
    if sp == 0:
        return [0.0] * N, 0.0
    sw = [sum(cov[i][j] * weights[j] for j in range(N)) for i in range(N)]
    return [weights[i] * sw[i] / sp for i in range(N)], sp


def risk_parity_weights(cov, iters=500, tol=1e-10):
    """Weights equalising marginal risk contributions. No matrix inversion, no expected returns —
    which is why it survives out-of-sample far better than mean-variance optimisation."""
    N = len(cov)
    w = [1.0 / N] * N
    for _ in range(iters):
        rc, sp = marginal_risk_contribution(w, cov)
        if sp == 0:
            break
        target = sp / N
        new = [max(1e-12, w[i] * (target / rc[i]) if rc[i] > 0 else w[i]) for i in range(N)]
        s = sum(new)
        new = [x / s for x in new]
        if max(abs(new[i] - w[i]) for i in range(N)) < tol:
            w = new
            break
        w = new
    return w


def equal_weights(n):
    return [1.0 / n] * n

# ---------------------------------------------------------------- 11. eigen (Jacobi)


def eigen_symmetric(A, iters=100, tol=1e-12):
    """Jacobi rotation eigendecomposition for a symmetric matrix.
    Returns (eigenvalues desc, eigenvectors as columns aligned to those eigenvalues)."""
    N = len(A)
    a = [row[:] for row in A]
    v = [[1.0 if i == j else 0.0 for j in range(N)] for i in range(N)]
    for _ in range(iters):
        off = max((abs(a[i][j]), i, j) for i in range(N) for j in range(N) if i != j) \
            if N > 1 else (0.0, 0, 0)
        if off[0] < tol:
            break
        _, p, q = off
        if a[p][p] == a[q][q]:
            theta = math.pi / 4
        else:
            theta = 0.5 * math.atan2(2 * a[p][q], a[p][p] - a[q][q])
        c, s = math.cos(theta), math.sin(theta)
        for k in range(N):
            akp, akq = a[k][p], a[k][q]
            a[k][p], a[k][q] = c * akp + s * akq, -s * akp + c * akq
        for k in range(N):
            apk, aqk = a[p][k], a[q][k]
            a[p][k], a[q][k] = c * apk + s * aqk, -s * apk + c * aqk
        for k in range(N):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p], v[k][q] = c * vkp + s * vkq, -s * vkp + c * vkq
    vals = [a[i][i] for i in range(N)]
    order = sorted(range(N), key=lambda i: -vals[i])
    return [vals[i] for i in order], [[v[r][i] for i in order] for r in range(N)]


def pc1_share(cov):
    """λ₁/Σλ — the fraction of variance in the first principal component.
    A 'diversified' book loading mostly on PC1 is one position with extra fees."""
    vals, _ = eigen_symmetric(cov)
    tot = sum(abs(x) for x in vals)
    return (abs(vals[0]) / tot) if tot > 0 else float('nan')


def condition_number(cov):
    """κ = λ_max/λ_min. Inversion multiplies estimation error by this. Financial covariance
    matrices are near-singular, which is the whole reason 'optimal' weights are unstable."""
    vals, _ = eigen_symmetric(cov)
    vals = [abs(x) for x in vals]
    lo = min(vals)
    return (max(vals) / lo) if lo > 1e-15 else float('inf')


def marchenko_pastur_edge(n_assets, n_obs, sigma2=1.0):
    """Eigenvalues below (1+√q)²σ² with q=N/T are indistinguishable from noise.
    With small T that is nearly all of them — clip before using the matrix for anything."""
    if n_obs <= 0:
        return float('inf')
    q = n_assets / n_obs
    return sigma2 * (1 + math.sqrt(q)) ** 2


def noise_eigenvalue_share(cov, n_obs):
    """Fraction of eigenvalues that fall below the Marchenko-Pastur noise edge."""
    vals, _ = eigen_symmetric(cov)
    avg_var = sum(cov[i][i] for i in range(len(cov))) / len(cov)
    edge = marchenko_pastur_edge(len(cov), n_obs, avg_var)
    return sum(1 for v in vals if abs(v) < edge) / len(vals), edge

# ---------------------------------------------------------------- 10. shrinkage


def shrink_cov(cov, delta=None, n_obs=None):
    """Shrink towards a scaled identity (constant-variance target):
    Σ_shrunk = (1−δ)Σ + δ·μ·I, μ = mean(diag Σ).
    If δ is None, use a Ledoit-Wolf-flavoured heuristic δ = min(1, N/T). Shrinkage is what makes
    a covariance matrix safe to invert; raw sample matrices are not."""
    N = len(cov)
    mu = sum(cov[i][i] for i in range(N)) / N
    if delta is None:
        delta = min(1.0, N / max(1, n_obs)) if n_obs else 0.2
    return [[(1 - delta) * cov[i][j] + (delta * mu if i == j else 0.0) for j in range(N)]
            for i in range(N)], delta

# ---------------------------------------------------------------- 15. constrained weights


def solve_min_variance(cov, cap=None, iters=4000, step=0.02):
    """Long-only minimum-variance weights with an optional per-asset cap, by projected gradient.
    Deliberately avoids inverting Σ (see condition_number). Returns (weights, σ_p).
    The Lagrange/KKT view: the cap that binds is information about the strategy's appetite,
    not a bug to relax."""
    N = len(cov)
    w = [1.0 / N] * N
    for _ in range(iters):
        g = [2 * sum(cov[i][j] * w[j] for j in range(N)) for i in range(N)]
        gm = sum(g) / N
        w = [w[i] - step * (g[i] - gm) for i in range(N)]
        w = [max(0.0, x) for x in w]
        w = _project_simplex_cap(w, cap, N)
    return w, portfolio_vol(w, cov)


def _project_simplex_cap(w, cap, N):
    """Normalise to sum 1 while never exceeding `cap`. Naively capping AFTER normalising lets the
    renormalisation push weights back above the cap — so redistribute the excess by water-filling
    onto the uncapped positions until nothing exceeds it."""
    s = sum(w)
    w = [1.0 / N] * N if s <= 0 else [x / s for x in w]
    if cap is None:
        return w
    if cap * N < 1.0 - 1e-12:                     # cap too tight to sum to 1
        return [cap] * N
    for _ in range(N + 2):
        over = [i for i, x in enumerate(w) if x > cap + 1e-15]
        if not over:
            break
        excess = sum(w[i] - cap for i in over)
        for i in over:
            w[i] = cap
        free = [i for i in range(N) if w[i] < cap - 1e-15]
        if not free:
            break
        room = sum(cap - w[i] for i in free)
        for i in free:                            # distribute in proportion to remaining room
            w[i] += excess * ((cap - w[i]) / room) if room > 0 else excess / len(free)
    return w


def binding_caps(weights, cap, tol=1e-6):
    """Which caps bind at the solution — read these before relaxing any risk limit."""
    return [i for i, x in enumerate(weights) if cap is not None and x >= cap - tol]

# ---------------------------------------------------------------- reporting


def report(returns, weights=None, cap=None, label="portfolio"):
    names, rows, T = _as_matrix(returns)
    N = len(rows)
    _, C = cov_matrix(returns)
    w = weights if weights else equal_weights(N)
    rc, sp = marginal_risk_contribution(w, C)
    _, R = corr_matrix(returns)
    _, L, lam, q = tail_dependence(returns)
    off = [R[i][j] for i in range(N) for j in range(N) if i != j]
    noise_share, edge = noise_eigenvalue_share(C, T)
    rp = risk_parity_weights(C)
    L = [f"=== {label} · N={N} assets, T={T} observations ===",
         f"  portfolio vol (given weights) {sp:.5f} · sum of risk contributions {sum(rc):.5f}",
         f"  risk contributions: " + " · ".join(f"{names[i]} {100*rc[i]/sp:.0f}%" for i in range(N))
         if sp else "  risk contributions: undefined",
         f"  mean off-diagonal correlation {st.mean(off):+.3f}",
         f"  tail dependence lambda {lam:.3f} vs {q:.3f} if independent → "
         f"{'assets FAIL TOGETHER ' + f'{lam/q:.1f}x more than independence' if lam > q * 1.3 else 'close to independent'}",
         f"  PC1 share of variance {100*pc1_share(C):.0f}% "
         f"({'concentrated — effectively one bet' if pc1_share(C) > 0.6 else 'diversified'})",
         f"  condition number {condition_number(C):.1f} "
         f"({'do NOT invert this matrix' if condition_number(C) > 100 else 'invertible with care'})",
         f"  T/N ratio {T/N:.1f} (want >= 10) · eigenvalues below the MP noise edge "
         f"({edge:.2e}): {100*noise_share:.0f}%",
         f"  risk-parity weights: " + " · ".join(f"{names[i]} {100*rp[i]:.0f}%" for i in range(N))]
    if cap:
        mv, mvvol = solve_min_variance(C, cap=cap)
        b = binding_caps(mv, cap)
        L.append(f"  min-variance (long-only, cap {cap:.0%}) vol {mvvol:.5f} · "
                 f"binding caps: {[names[i] for i in b] or 'none'}")
    L.append("  note: covariance is second-order only — it cannot represent default, depeg, "
             "oracle failure or a venue going to zero. Pair with a scenario table.")
    return "\n".join(L)


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(description="qportfolio — pass a JSON dict {asset: [returns…]}")
    ap.add_argument("--json", required=True, help="path to a JSON file or an inline JSON string")
    ap.add_argument("--cap", type=float, default=None)
    a = ap.parse_args()
    try:
        data = json.load(open(a.json))
    except (OSError, ValueError):
        data = json.loads(a.json)
    print(report(data, cap=a.cap))

# ---------------------------------------------------------------- 61. breadth con guardia


def effective_breadth(returns, weights=None):
    """Scommesse EFFICACI, con la guardia che serve quando rho medio e' NEGATIVA.

    ERRORE COMMESSO (2026-07-30): ho applicato N/(1+(N-1)*rho) a 7 strategie sullo STESSO
    strumento con rho medio −0.067 e ottenuto **11.71 su 7**. Numero privo di significato
    come "scommesse indipendenti", e la mia stessa prosa lo contraddiceva due righe sotto.

    La formula misura RIDUZIONE DI VARIANZA, non "scommesse indipendenti con rendimento
    atteso". Con rho negativa da > N, ed e' corretta come riduzione di varianza — ma leggerla
    come breadth e' un errore, perche' due strategie che si ANNULLANO hanno rho = −1 e
    "breadth infinita" mentre pagano due commissioni per tenere zero.

    ⚠ NON TRASPORTARLA AL MULTIPLE TESTING (aggiunto 2026-08-01, correzione presa da una
    sessione parallela). `N/(1+(N-1)rho)` e' la dimensione efficace per la VARIANZA DI UNA
    MEDIA di variabili equicorrelate. Usarla per correggere la soglia sul MASSIMO di N test
    e' un'euristica, non un risultato: la direzione e' giusta, la magnitudine no. Per la
    penalita' da ricerca si usa RAS (`qvalidate.rademacher_complexity` +
    `effective_trials`), che la MISURA dalle serie invece di dedurla da un rho ipotizzato.

    LA FIRMA DELL'ANNULLAMENTO, da riportare sempre insieme al numero:
        rho medio <= 0  E  PC1 alta  →  non e' decorrelazione, e' cancellazione reciproca.
    In quel caso la domanda giusta non e' la breadth ma: quanto dell'esposizione LORDA
    sopravvive nel netto, e i costi si sommano mentre i rendimenti si annullano?
    """
    names, rows, T = _as_matrix(returns)
    N = len(rows)
    if N < 2:
        return {"n": N, "usable": False, "note": "servono >= 2 serie"}
    _, C = corr_matrix(returns)
    off = [C[i][j] for i in range(N) for j in range(N) if i != j]
    rho = st.mean(off)
    _, cov = cov_matrix(returns)
    pc1 = pc1_share(cov)
    w = [1.0 / N] * N if weights is None else [float(x) for x in weights]
    sw = sum(abs(x) for x in w) or 1.0
    w = [x / sw for x in w]
    sd = [math.sqrt(cov[i][i]) if cov[i][i] > 0 else 0.0 for i in range(N)]
    num = sum(abs(w[i]) * sd[i] for i in range(N)) ** 2
    den = sum(w[i] * w[j] * cov[i][j] for i in range(N) for j in range(N))
    var_ratio = (num / den) if den > 0 else float("nan")
    naive = N / (1.0 + (N - 1) * rho) if (1.0 + (N - 1) * rho) > 0 else float("inf")
    cancelling = (rho <= 0.0) and (pc1 > 0.35)
    out = {"n": N, "mean_corr": rho, "pc1_share": pc1,
           "variance_reduction_ratio": var_ratio, "naive_formula": naive,
           "min_corr": min(off), "max_corr": max(off),
           "valid_as_breadth": bool(rho > 0), "cancellation_signature": bool(cancelling)}
    if rho <= 0:
        out["breadth"] = None
        out["note"] = (f"rho medio {rho:+.4f} <= 0: la formula da {naive:.2f} su {N} e NON va "
                       "letta come breadth. E' riduzione di varianza. "
                       + ("FIRMA DI CANCELLAZIONE (rho<=0 con PC1 al "
                          f"{100*pc1:.0f}%): misura invece quanta esposizione lorda "
                          "sopravvive nel netto — i costi si sommano, i rendimenti si "
                          "annullano." if cancelling else
                          "Verifica se le serie si annullano a vicenda."))
    else:
        out["breadth"] = naive
        out["note"] = f"breadth {naive:.2f} su {N} ({100*naive/N:.0f}% di N)"
    return out


def diversification_sign_check(mean_per_period, label="portafoglio"):
    """Diversificare con mu < 0 PEGGIORA lo Sharpe. Controllo del segno prima di applicarla.

    MISURATO (2026-07-30): 7 strategie con Sharpe medio −0.87 aggregate danno −2.74, cioe'
    PEGGIO della media delle parti. Non e' un paradosso: la diversificazione riduce il
    DENOMINATORE, e con numeratore negativo amplifica il segno. Se la sd fosse rimasta a
    livello di gamba lo Sharpe sarebbe stato −0.77 invece di −2.74.
    """
    ok = mean_per_period > 0
    return {"mean": mean_per_period, "safe_to_diversify": bool(ok),
            "note": ("ok: mu > 0, ridurre la varianza migliora lo Sharpe" if ok else
                     "⚠ mu <= 0: diversificare RIDUCE il denominatore e PEGGIORA lo Sharpe. "
                     "Non applicare tecniche di riduzione della varianza a un aggregato con "
                     "expectancy negativa — prima si risolve il segno.")}
