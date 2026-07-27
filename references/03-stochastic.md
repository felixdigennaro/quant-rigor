# Formulas 16-22 — Calculus & Stochastic Processes: why paths behave the way they do

This module contains the **theoretical reason** for the most expensive practical mistake in
strategy testing: believing a path-dependent rule that your data cannot resolve. Read §21 before
ever testing a trailing stop, a breakeven-move, or anything defined on the running max/min.

---

## 16. Taylor Expansion

`f(x) = f(a) + f'(a)(x−a) + ½f''(a)(x−a)² + … + Rₙ`

**What it does:** approximates a function locally by a polynomial. The organising principle of all
sensitivity analysis: first order = direction, second order = curvature.

**Use in practice:**
- The greeks **are** a Taylor expansion of option value:
  `ΔV ≈ Δ·dS + ½Γ·dS² + Θ·dt + 𝜈·dσ`. Delta-hedging removes the first term; **gamma is what's
  left**, and that's why a hedged book still bleeds or profits on big moves.
- Bond risk: `ΔP/P ≈ −D·Δy + ½C·Δy²` (duration + convexity).
- Log-return approximation `ln(1+r) ≈ r − r²/2`: fine for `r = 1%`, wrong for `r = 30%`.
  **Where the `−σ²/2` in GBM comes from** (§20).

**How it lies:** the remainder is not small when moves are large — exactly in the scenarios that
matter. A delta-gamma approximation that's excellent for ±1% is badly wrong for a −34% crash.
Truncating at second order **systematically understates tail P&L**.

**Do this:** use Taylor for intuition and hedging in normal times; for tail risk use **full
revaluation** (recompute the actual P&L at the stressed price), never the expansion.

---

## 17. Partial Derivatives

`∂f/∂xᵢ = lim_{h→0} [f(…,xᵢ+h,…) − f(…,xᵢ,…)]/h` — one variable moves, others held fixed

**What it does:** isolates one driver's effect.

**Use in practice:** every greek is a partial derivative. Every "what if we change one parameter"
question is a partial derivative. Compute them numerically on the backtest:
`∂metric/∂param` at the operating point.

**How it lies — "others held constant" is usually false in markets.** Vol rises when price falls;
funding spikes when everyone is on one side; liquidity evaporates exactly when you need the exit.
A partial derivative computed as if these are independent will misprice the joint move. The
correct object is often a **total** derivative (§18) including the induced changes.

**Do this:** after computing a one-parameter sensitivity, ask what else moves with it *in reality*,
and test the joint change. For a stop: widening the stop also changes hold time, which changes
funding cost and the probability of hitting the time-flat — three effects, not one.

---

## 18. Chain Rule

`dy/dx = (dy/du)(du/dx)`

**What it does:** composes sensitivities through a pipeline.

**Use in practice:** risk flows through chains, and the chain rule is how you find the real
exposure. Example: a yield-funded overlay — `∂(fund value)/∂(asset price)` =
`∂fund/∂overlay_position × ∂overlay_position/∂gold`. Since the overlay is built only from accrued
coupon, the first factor is small and bounded — **that's the structural protection, expressed as a
derivative**. Another: strategy P&L → account equity → position size (if size is recomputed from
equity) → next period's P&L. That last chain is a **feedback loop**, and it's why compounding
sizing changes the return distribution's shape, not just its scale.

**How it lies:** people stop the chain too early and miss a link. A classic: sizing computed from `withdrawable` (which is net of margin already in use) creates a hidden link
`position → margin used → withdrawable → next size`, so intraday re-arms shrink mechanically.
The chain was there in the code; nobody had differentiated it.

**Do this:** write the full dependency chain for any feedback path, especially anything where
**position size depends on equity or free margin**. Then decide deliberately whether you want that
feedback (Kelly-like) or not (fixed notional).

---

## 19. Itô's Lemma

For `dX = μdt + σdW`:
`df(X,t) = (∂f/∂t + μ∂f/∂x + ½σ²∂²f/∂x²)dt + σ(∂f/∂x)dW`

**What it does:** the chain rule for stochastic processes. The extra `½σ²f''` term appears because
`(dW)² = dt` — randomness is not negligible at second order.

**Use in practice:**
- Derives Black-Scholes (§25) and every no-arbitrage pricing result.
- **Explains volatility drag directly:** apply Itô to `f = ln S` under GBM and you get
  `d(ln S) = (μ − σ²/2)dt + σdW`. The `−σ²/2` is not a fee, it's a mathematical consequence of
  compounding a random process. **This single term is why leverage can raise expected return and
  lower realised wealth** — the mechanism behind most retracted leveraged theses.
- Any non-linear function of a random price (a payoff, a rebalanced portfolio, a leveraged ETF)
  has a drift you didn't put there. Leveraged ETF decay is the `½σ²` term.

**How it lies:** Itô assumes continuous paths and finite variance — no jumps. With gaps
(weekend opens, halts, oracle updates) the correct calculus is jump-diffusion, and the Itô-derived
hedge under-hedges precisely at the gap.

**Do this:** whenever a strategy is non-linear in price *and* rebalanced/compounded, compute the
induced drift term explicitly instead of assuming the average of the parts. Then check the gap case
separately.

---

## 20. Geometric Brownian Motion

`dS = μS dt + σS dW` → `S_t = S₀·exp[(μ − ½σ²)t + σW_t]`

**What it does:** the baseline model of a positive price: log-normal, proportional moves,
constant vol.

**Use in practice — read the two drifts and never confuse them:**
- **Mean** of `S_t`: `S₀e^{μt}` (arithmetic).
- **Median** of `S_t`: `S₀e^{(μ−σ²/2)t}` (what a typical path does).
- With `μ = 10%`, `σ = 40%`: mean grows 10%/yr, **median grows 2%/yr**. Report both in every
  compounding projection, and never present the mean path as "expected outcome".
- Sizing: max growth of the median is fractional-Kelly; beyond `f* = μ/σ²`, more leverage
  **lowers** median wealth while raising the mean. That is the entire quantitative case against
  leverage that only looks good in-sample.

**How it lies:** constant `μ`, `σ`; no jumps; independent increments; log-normal tails
(too thin). Real returns have vol clustering (GARCH), fat tails, and regime shifts. GBM is a
**null model to reason against**, not a forecast.

**Do this:** use GBM to compute the drag and the Kelly bound analytically, then verify empirically
with a **block bootstrap** of real returns (preserves autocorrelation and fat tails) instead of
simulating GBM paths.

---

## 21. Wiener Process — ★ the resolution limit (read before testing any path-dependent rule)

`W₀ = 0`; `W_t − W_s ~ N(0, t−s)`; independent increments; continuous paths a.s.

**What it does:** the canonical continuous random walk. Two properties that decide what you can
measure:

1. **Quadratic variation is finite:** `Σ(ΔW)² → t`. This is why variance/vol is a well-defined,
   measurable, resolution-**stable** quantity.
2. **Total variation is infinite:** `Σ|ΔW| → ∞` as `Δt → 0`. The measured path length grows
   like `√(1/Δt)` and **never converges**.

**Use in practice:** it is the null model for "is this move meaningful?". Under a driftless walk,
the expected absolute move over `t` is `σ√(2t/π)` and the expected range grows like `√t` — so any
observed excursion should be compared against `√t`, not against a fixed number. It also sets the
**vol scaling** (`σ_T = σ_1·√T`) used everywhere, and gives the measurability test below.

**How it lies:** the paths are continuous and non-differentiable, which real prices are not — real
markets **jump** (gaps, halts, liquidation cascades, oracle updates). So the Wiener null
*understates* tail moves while *overstating* the number of small oscillations relative to a market
that trades in discrete ticks with a spread. It also assumes independent increments: real returns
show vol clustering, so quiet periods under-predict and violent periods over-predict the path length.
Treat it as the reference process, never as the market.

### The consequence that costs money
- A rule that depends on **fixed levels** (limit TP, fixed stop, "did price touch X") depends on
  the *extremes* of the path, which converge fast. Measure it at 5m, you get the same answer at 1m.
  **Resolution-safe.**
- A rule that depends on the **path itself** — trailing stop, running max/min, breakeven-move,
  "highest close since entry", most volatility-scaled exits — depends on total variation, which
  **diverges**. Every refinement of the data reveals more oscillation, so the measured performance
  keeps changing. **Not resolution-safe.**

### The practical test (this is §2 of SKILL.md, here's why it works)
- Refinement factor from resolution `Δ` to `Δ/n`: for a pure random walk,
  `Σ|ΔW|` grows by `√n`. Measure it: `sum(fine ranges)/coarse range`. If ≈ `√n`, the path has
  **no intrinsic scale** and there is no "fine enough" data — only rules whose distance `D` is
  large versus the *finest available* bar range are testable.
- Rule of thumb: need `D ≥ 3-5 × median bar range` at the finest resolution you have.
- Corollary: making the trail wider until it's resolvable usually makes it ≈ the existing stop,
  i.e. the idea collapses into what you already do. **That is a legitimate finding, not a failure.**

**Do this:** classify every exit rule as fixed-level or path-dependent **before** testing.
For path-dependent rules, either (a) prove `D ≥ 3-5×` finest bar range, or (b) declare it not
measurable and offer fixed-level alternatives that address the same goal.

---

## 22. Stochastic Differential Equation (SDE)

`dX_t = μ(X_t,t)dt + σ(X_t,t)dW_t`

**What it does:** the general template: a deterministic drift plus a scaled random shock. Choosing
`μ(·)` and `σ(·)` *is* choosing your market model.

**Use in practice — the model choice is the thesis:**
- `μ = μX`, `σ = σX` → GBM: trends persist, no reversion. **A momentum thesis.**
- `μ = θ(m − X)` → Ornstein-Uhlenbeck: mean reversion to `m` at speed `θ`, stationary.
  **A fade/mean-reversion thesis** — this is the correct null for any fade: an intraday zone fade is literally an OU bet on reversion. The OU half-life `ln2/θ` gives the natural holding period, and therefore whether a fixed time-flat is early, right, or late.
- Add jumps → jump-diffusion (gaps, liquidations). Add `√X` → CIR (rates, funding, positive).
- Vol as its own SDE → Heston (vol clustering).

**Use it as a discipline:** write down the SDE your strategy implicitly assumes. Then check whether
your *modifications* are consistent with it. Concrete: extending a fade's take-profit **beyond**
the reversion target is no longer an OU bet — it requires overshoot, i.e. it silently switches the
thesis from mean-reversion to momentum. That's a different SDE and needs its own evidence.

**How it lies:** the parameters are estimated, unstable, and often not identifiable from short
samples (`θ` for OU needs many mean-reversion cycles). Fitting an SDE to 15 sessions gives
parameters with enormous error bars, and a fitted model with tight-looking output is the most
seductive way to be wrong.

**Do this:** state the implicit SDE, estimate its parameters **with bands**, and check that every
proposed change is consistent with the same process. If a change requires a different SDE, say so
explicitly — it needs its own out-of-sample evidence, not the parent strategy's.
