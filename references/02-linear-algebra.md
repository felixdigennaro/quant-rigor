# Formulas 9-15 — Linear Algebra & Optimisation: portfolios, risk decomposition, parameter fitting

These turn "a set of positions" into "a risk object you can reason about", and "a set of knobs"
into "an optimisation you can distrust properly". The recurring danger: **every one of them
amplifies estimation error**, and the amplification is invisible in the output.

---

## 9. Matrix Multiplication

`C = AB`, `cᵢⱼ = Σₖ aᵢₖbₖⱼ`

**What it does:** composes linear transformations. In finance it's the machinery of aggregation.

**Use in practice:**
- Portfolio return: `r_p = wᵀr`. Portfolio variance: `σ²_p = wᵀΣw` — one expression that
  contains every pairwise interaction (this is the whole point of the matrix form).
- Factor model: `r = Bf + ε`, so `Σ = BΣ_f Bᵀ + D`. This decomposition is how you separate
  "market risk I'm paid for" from "idiosyncratic risk I should diversify".
- Scenario matrix: positions × scenarios → P&L vector, in one multiply. **This is the cheapest
  and most useful risk tool that exists** and it needs no distributional assumption.

**How it lies:** it doesn't lie, it **propagates**. Garbage in `Σ` becomes confident-looking
`σ_p`. And `wᵀΣw` is a sum of `n²` terms: with `n` sleeves you are trusting `n(n+1)/2`
estimated numbers, most of which you never look at.

**Do this:** build the scenario matrix by hand for the 5-10 scenarios that would actually hurt
(venue failure, funding spike, index −10% overnight, stablecoin depeg, oracle break) and read the
P&L column. No estimation, no normality, no covariance. Do this **before** any optimiser.

---

## 10. Matrix Inverse

`A⁻¹A = I`, `A⁻¹ = adj(A)/det(A)`

**What it does:** solves linear systems; appears in every optimal-weight formula
(Markowitz: `w* ∝ Σ⁻¹μ`; OLS: `β = (XᵀX)⁻¹Xᵀy`).

**Use in practice:** mean-variance weights, regression betas, Kalman updates.

**How it lies — the most dangerous item in this file:**
- **Inversion amplifies noise by the condition number** `κ(Σ) = λ_max/λ_min`. Financial
  covariance matrices are near-singular (assets are highly correlated), so `κ` is huge and
  `Σ⁻¹` multiplies your estimation error by it. Markowitz's famous instability is not a
  Markowitz flaw, it's this.
- **`n` too small relative to dimension.** With `T` observations and `N` assets, the sample
  covariance is not even invertible for `T < N`, and is badly conditioned for `T < 10N`.
  15 sessions and 3 rungs is already marginal.
- **Result: "optimal" weights that are extreme, unstable, and worse out-of-sample than 1/N.**
  Equal weighting beats optimised weights in most honest out-of-sample tests.

**Do this:** never invert a raw sample covariance for a decision. Options, in order of preference:
(1) **1/N or risk-parity** — no inversion, no estimation of μ; (2) **shrinkage**
(Ledoit-Wolf: `Σ_shrunk = (1−δ)Σ_sample + δ·target`); (3) constrained optimisation with
position caps that bound the damage. Always report `κ(Σ)` if you invert anything.
Check stability: perturb inputs by their standard error and re-solve. If weights move a lot,
the solution is noise.

---

## 11. Eigenvalues & Eigenvectors

`Av = λv`

**What it does:** finds the directions a transformation only stretches. For a covariance matrix,
eigenvectors are **uncorrelated risk factors** (PCA) and eigenvalues are the variance in each.

**Use in practice:**
- **PCA on returns:** in crypto/equities PC1 is usually "the market" and typically carries
  60-90% of variance. If your "diversified" book loads mostly on PC1, you have one position.
- **Risk concentration diagnostic:** the fraction `λ₁/Σλ` is a one-number answer to
  "am I actually diversified?" — far more honest than counting positions.
- **Conditioning:** `κ = λ_max/λ_min` (see §10). Tiny eigenvalues = directions where the
  optimiser will happily take enormous, meaningless bets.

**How it lies:**
- **Marchenko-Pastur**: with finite `T`, most small eigenvalues of a sample covariance are
  **pure noise**. The MP edge for `q = N/T` is `λ_max_noise = (1+√q)²·σ²` — eigenvalues below
  it carry no information. With `T` small, that's nearly all of them.
- **Eigenvectors rotate across regimes.** PC1 in a calm month is not PC1 in a crash.
- **Sign/order are not stable** across samples; don't attach a story to PC3.

**Do this:** report `λ₁/Σλ` (concentration) and clip eigenvalues below the MP edge before
using `Σ` for anything. If a factor's eigenvalue is below the noise edge, it does not exist.

---

## 12. Covariance Matrix

`Σ = E[(X−μ)(X−μ)ᵀ]`, `Σᵢⱼ = Cov(Xᵢ,Xⱼ)`, variances on the diagonal

**What it does:** the full second-order risk description of a book.

**Use in practice:**
- `σ_p = √(wᵀΣw)`; **marginal risk contribution** of position i: `(Σw)ᵢ·wᵢ/σ_p`
  — the number that tells you *who* is actually carrying your risk (usually not who you think).
- Risk parity: choose `w` so all marginal contributions are equal.
- Build `Σ` over your sleeves (e.g. base yield, overlay, opportunistic) and check they are not secretly one bet on the same liquidity factor.

**How it lies:**
- **Estimation error dominates** for anything but tiny `N` (see §10, §11).
- **Correlations rise in stress** — the matrix you need is the conditional one. But measure it
  correctly: computing Pearson correlation *inside the worst decile of the portfolio mean* is
  **biased downward** by the conditioning itself (Boyer-Gibson-Loretan). Truncating the joint
  distribution mechanically shrinks the sample correlation, so a falling number there is not
  evidence that diversification improves in stress — a trap that reads as good news.
  Use the **tail dependence coefficient** instead: `λ = P(X ≤ q_α | Y ≤ q_α)`, estimated
  non-parametrically. Independent series give `λ ≈ α`; `λ >> α` means the assets fail together.
  `qportfolio.tail_dependence()` implements it; `tail_conditional_corr()` is kept only to show
  the biased number most reports quote.
- **Second-order only**: `Σ` cannot see skew, jumps, or a venue going to zero. A lending sleeve
  has near-zero variance right up until an exploit takes 100% — variance is the *wrong*
  risk measure for that shape (use scenario P&L, §9).
- Non-synchronous / different-frequency series bias covariance toward zero.

**Do this:** always pair `Σ`-based risk with a **scenario table**. State explicitly which risks
`Σ` cannot represent (default, depeg, oracle, regulatory) — for a mostly-lending fund that's most
of the real risk.

---

## 13. Gradient

`∇f = [∂f/∂x₁, …, ∂f/∂xₙ]ᵀ` — points in the direction of steepest increase

**What it does:** the local direction to improve; the engine of every fitting procedure.

**Use in practice:** parameter fitting, and more usefully **sensitivity analysis**: `∂metric/∂param`
tells you which knobs matter. If the gradient of P&L w.r.t. a parameter is tiny, stop arguing
about that parameter. In options, the greeks *are* gradients (§26, §27).

**How it lies:**
- **Gradient ascent on a backtest metric is literally overfitting.** The surface you're climbing
  is `true_signal + noise`, and with 15 observations the noise gradient dominates. You will
  reliably find a peak that is a noise artefact.
- Discontinuous/path-dependent objectives (stops, triggers) have **no useful gradient** — the
  surface is a step function; finite differences give garbage that depends on step size.
- A local optimum is not global, and on a noisy surface "global" is meaningless anyway.

**Do this:** use gradients for **sensitivity**, not for choosing parameters. Choose parameters
from the **middle of a plateau** on a coarse grid (see charts §4), never from an optimiser's peak.
If a parameter's plateau is narrow, the parameter is not a parameter — it's a fit.

---

## 14. Hessian Matrix

`H = ∇²f`, `Hᵢⱼ = ∂²f/∂xᵢ∂xⱼ` — curvature

**What it does:** classifies critical points (positive-definite = min, negative-definite = max,
mixed = saddle) and quantifies **how sharp** an optimum is.

**Use in practice — this is the overfit detector:**
- The Hessian's eigenvalues at your chosen parameters *are* the sharpness of the peak.
  **Large curvature = fragile fit.** Flat = robust.
- Second-order risk: Gamma is the Hessian of option value w.r.t. price (§27).
- Newton steps `−H⁻¹∇f` converge fast and overfit faster.

**How it lies:** on a noisy objective the Hessian is mostly measuring **noise curvature**;
a beautiful sharp maximum is exactly what fitting noise produces. Also, near-zero eigenvalues mean
the parameter is unidentified — the optimiser's chosen value there is arbitrary, and reporting it
to 3 decimals is false precision.

**Do this:** at the chosen parameter set, estimate curvature numerically along each axis and
report it as a **stability statement**: "metric drops X% if the parameter moves ±1 step". If X is
large, do not adopt. Prefer a parameter where the metric is flat within ±1-2 steps.

---

## 15. Lagrange Multipliers

`ℒ(x,λ) = f(x) + Σλᵢgᵢ(x)`, solve `∇ₓℒ = 0`, `gᵢ(x) = 0`
(inequalities → KKT: `λᵢ ≥ 0`, complementary slackness)

**What it does:** optimise subject to constraints. The multiplier `λᵢ` is the **shadow price** of
constraint i — how much the objective would improve per unit of relaxation.

**Use in practice:**
- Portfolio with budget/exposure limits: `max wᵀμ − (γ/2)wᵀΣw s.t. Σw = 1, w ≥ 0, wᵢ ≤ cap`.
- **Read the multipliers, not just the weights.** Hard risk caps (e.g. a ceiling on the strategy sleeve, a per-venue limit) are constraints; `λ` on each tells you which
  constraint is actually costing return. If a cap's `λ ≈ 0`, it's free — keep it. That is a
  quantitative argument for keeping a risk limit, which is otherwise a matter of taste.
- Risk budgeting: allocate to equalise marginal risk-adjusted return.

**How it lies:**
- Multipliers inherit all the estimation error of `μ` and `Σ` — a shadow price computed from a
  noisy `μ` is noise. Report `λ` with a band or not at all.
- **Constraints are load-bearing, not obstacles.** A constraint that binds is often the only thing
  standing between you and the optimiser's fantasy. "Relax it, `λ` says we'd earn more" is how
  blow-ups get authorised. Relax a risk constraint only with an out-of-sample argument.
- Equality constraints in code that should be inequalities (forcing full investment when cash is
  a legitimate position) silently removes the best answer.

**Do this:** write constraints explicitly, solve, then report **which bind** and their `λ`.
Treat a binding risk cap as information about the strategy's appetite, not as a bug to fix.
