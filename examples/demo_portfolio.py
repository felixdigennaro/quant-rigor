#!/usr/bin/env python3
"""quant-rigor portfolio demo — STDLIB ONLY. Run: python3 examples/demo_portfolio.py

Third domain: a multi-sleeve book rather than a single strategy. Exercises formulas 9-15
(covariance, eigenvalues, risk contributions, shrinkage, constrained weights) and shows the three
ways a "diversified" portfolio is not diversified:

  1. risk contributions that bear no relation to capital weights
  2. correlations that rise in the worst decile — exactly when diversification is needed
  3. a covariance matrix too noisy to optimise on (T/N too small, condition number too large)

Deterministic: same numbers on any machine.
"""
import math
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import qportfolio as P    # noqa: E402
import qstats as S        # noqa: E402
import qsize as Z         # noqa: E402


class LCG:
    def __init__(self, seed=424242):
        self.s = seed

    def u(self):
        self.s = (1103515245 * self.s + 12345) % (2 ** 31)
        return self.s / (2 ** 31)

    def n(self):
        return math.sqrt(-2 * math.log(max(1e-12, self.u()))) * math.cos(2 * math.pi * self.u())


def build_book(T=250, seed=424242):
    """Four sleeves with deliberately different shapes:
      base      — low vol, steady, tiny positive drift (a lending-like carry sleeve)
      carry     — moderate vol, correlated with `base` through a shared funding factor
      momentum  — high vol, uncorrelated in calm periods
      convex    — small positive drift most of the time, rare large loss (short-gamma shape)
    A common stress factor hits everything in the worst ~10% of periods. This is DELIBERATELY
    constructed so the tail-correlation effect is visible: real books behave this way, but here it
    is built in rather than discovered, which is the honest way to use a synthetic example.
    """
    r = LCG(seed)
    out = {"base": [], "carry": [], "momentum": [], "convex": []}
    for t in range(T):
        funding = r.n()                              # shared factor, always on
        stress = 1.0 if r.u() < 0.10 else 0.0        # ~10% of periods are common-stress periods
        common = -(1.0 + abs(r.n())) if stress else 0.0
        out["base"].append(0.00095 + 0.0015 * r.n() + 0.0008 * funding + 0.004 * common)
        out["carry"].append(0.00170 + 0.0040 * r.n() + 0.0030 * funding + 0.009 * common)
        out["momentum"].append(0.00230 + 0.0110 * r.n() + 0.008 * common)
        conv = 0.00160 + 0.0020 * r.n()
        if r.u() < 0.02:                             # rare idiosyncratic loss (short-gamma shape)
            conv -= 0.03
        out["convex"].append(conv + 0.007 * common)
    return out


def main():
    print("=" * 74)
    print("quant-rigor — portfolio demo (formulas 9-15)")
    print("=" * 74)
    book = build_book()
    names = list(book)
    T = len(book[names[0]])
    _, C = P.cov_matrix(book)
    _, R = P.corr_matrix(book)

    print(f"\n[1] SLEEVES ({T} periods)")
    print(f"  {'sleeve':<10}{'mean':>10}{'sd':>9}{'skew':>8}{'worst':>9}{'ann. SR':>9}")
    for k in names:
        d = S.describe(book[k])
        sr = (d["mean"] / d["sd"] * math.sqrt(252)) if d["sd"] else float('nan')
        print(f"  {k:<10}{d['mean']:>10.5f}{d['sd']:>9.4f}{d['skew']:>+8.2f}"
              f"{d['worst']:>+9.3f}{sr:>9.2f}")

    print("\n[2] CAPITAL WEIGHTS vs RISK CONTRIBUTIONS — equal weight")
    ew = P.equal_weights(len(names))
    rc, sp = P.marginal_risk_contribution(ew, C)
    print(f"  {'sleeve':<10}{'capital':>10}{'risk':>10}{'ratio':>9}")
    for i, k in enumerate(names):
        share = 100 * rc[i] / sp if sp else float('nan')
        print(f"  {k:<10}{100*ew[i]:>9.0f}%{share:>9.0f}%{share/(100*ew[i]):>9.1f}x")
    print(f"  portfolio vol {sp:.5f} · equal capital is NOT equal risk")

    print("\n[3] RISK PARITY vs MINIMUM VARIANCE (long-only, 40% cap)")
    rp = P.risk_parity_weights(C)
    mv, mvvol = P.solve_min_variance(C, cap=0.40)
    bind = P.binding_caps(mv, 0.40)
    print(f"  {'sleeve':<10}{'equal':>9}{'risk-par':>10}{'min-var':>10}")
    for i, k in enumerate(names):
        print(f"  {k:<10}{100*ew[i]:>8.0f}%{100*rp[i]:>9.0f}%{100*mv[i]:>9.0f}%")
    print(f"  vol: equal {sp:.5f} · risk-parity {P.portfolio_vol(rp,C):.5f} · min-var {mvvol:.5f}")
    print(f"  binding caps at 40%: {[names[i] for i in bind] or 'none'}"
          f" — a binding cap is information about appetite, not a bug to relax")

    print("\n[4] DO THE SLEEVES FAIL TOGETHER?")
    _, L, lam, q = P.tail_dependence(book, quantile=0.10)
    _, Rt = P.tail_conditional_corr(book, quantile=0.10)
    off = [(names[i], names[j], R[i][j], L[i][j], Rt[i][j])
           for i in range(len(names)) for j in range(i + 1, len(names))]
    print(f"  {'pair':<22}{'correlation':>12}{'tail dep λ':>12}{'vs indep':>10}"
          f"{'[biased ρ|tail]':>17}")
    for a, b, c, l, ct in off:
        print(f"  {a+' / '+b:<22}{c:>+12.3f}{l:>12.3f}{l/q:>9.1f}x{ct:>+17.3f}")
    m1 = st.mean([x[2] for x in off])
    print(f"  mean correlation {m1:+.3f} · mean tail dependence {lam:.3f} vs {q:.2f} if independent")
    if lam > q * 1.3:
        print(f"  → the sleeves fail together {lam/q:.1f}x more often than independence implies;")
        print(f"    average correlation understates it. Size on the tail measure.")
    print(f"  the last column is Pearson correlation inside the worst decile. It is BIASED")
    print(f"  downward by conditioning on the mean's tail and must not be read as evidence that")
    print(f"  diversification improves in stress — that is why lambda is the number reported.")

    print("\n[5] IS THIS MATRIX USABLE?")
    N = len(names)
    kappa = P.condition_number(C)
    pc1 = P.pc1_share(C)
    share, edge = P.noise_eigenvalue_share(C, T)
    vals, _ = P.eigen_symmetric(C)
    print(f"  T/N = {T/N:.1f} (want >= 10) · condition number {kappa:.1f}")
    print(f"  PC1 carries {100*pc1:.0f}% of variance"
          f"{' — effectively one bet' if pc1 > 0.6 else ''}")
    print(f"  eigenvalues: {[f'{v:.2e}' for v in vals]}")
    print(f"  Marchenko-Pastur noise edge {edge:.2e} → {100*share:.0f}% of eigenvalues are noise")
    sh, delta = P.shrink_cov(C, n_obs=T)
    print(f"  shrinkage δ={delta:.2f} lowers condition number {kappa:.1f} → {P.condition_number(sh):.1f}")
    print(f"  {'SAFE to invert' if P.condition_number(sh) < 100 else 'still do NOT invert — use risk parity / 1-N'}")

    print("\n[6] SIZING THE BOOK")
    port = [sum(w * book[k][t] for w, k in zip(rp, names)) for t in range(T)]
    d = S.describe(port)
    print(f"  risk-parity book: mean {d['mean']:+.5f}/period · sd {d['sd']:.4f} · "
          f"skew {d['skew']:+.2f} · t {d['t']:+.2f}")
    dd = S.drawdown(S.equity_from_returns(port))
    vc = S.var_cvar(port)
    print(f"  max DD {dd['max_dd']:+.4f} · VaR95 {vc['VaR_empirical']:+.4f} · "
          f"CVaR {vc['CVaR_expected_shortfall']:+.4f} · worst {vc['worst']:+.4f}")
    print("  " + Z.report(d["mean"], d["sd"], dd_tolerated=0.15, label="book").replace("\n", "\n  "))

    print("\n[7] WHAT COVARIANCE CANNOT SEE")
    print("  Second-order risk only. Absent from every number above: counterparty default,")
    print("  a venue going to zero, oracle failure, a depeg, regulatory action, and liquidity")
    print("  evaporating in the periods where the tail correlation just rose. Those require a")
    print("  scenario table and position caps — not a risk number.")
    print("=" * 74)


if __name__ == "__main__":
    main()
