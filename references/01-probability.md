# Formulas 1-8 — Probability & Statistics: what a sample can and cannot support

These eight decide **whether you know anything at all**. Every wrong strategy conclusion I have
seen traces back to one of these being skipped, not to a pricing model being wrong.

---

## 1. Bayes' Theorem

`P(A|B) = P(B|A)·P(A) / P(B)`

**What it does:** updates a belief when evidence arrives. The posterior is the prior times the
likelihood ratio.

**Use in practice — the base-rate discipline.** The question is never "does this backtest look
good?" but "given the base rate of strategies that look good and then fail, what's the posterior
that this one is real?" Write it explicitly:

- Prior: fraction of proposed edges that survive out-of-sample. Empirically brutal — for
  self-generated, single-regime backtests, well under 20%. If your own track record includes a retracted thesis that once looked strong, that IS your base rate — use it.
- Likelihood ratio: how much more likely is this result under "real edge" than under "noise"?
  A `t = 1.5` result is only ~2-3× more likely under a real effect than under noise.
- Posterior: `0.2 prior × 2.5 LR` → odds 0.25:0.8 → **~24%**. Not a finding.

**How it lies:** people apply it to the *result* and forget the *search*. If you tested 60 variants,
the evidence isn't "the winner beat baseline" — it's "the best of 60 beat baseline", which is a
much weaker likelihood ratio (see SKILL §5 best-of-N). Also: reusing the same data to form the
prior and to update makes the posterior meaningless.

**Do this:** before reporting any positive backtest, state the prior, the LR, and the posterior in
one line. If the posterior is under ~50%, the honest label is "hypothesis", not "edge".

---

## 2. Conditional Probability

`P(A|B) = P(A ∩ B) / P(B)`

**What it does:** the probability of A restricted to the world where B happened.

**Use in practice:** almost every trading rule is a conditional statement. "After a 100pt extension
beyond the zone, the price returns to the edge" is `P(return | extension ≥ 100)`. The whole edge
lives in the conditioning event, so **define it exactly** and count `P(B)` — how often the condition
even occurs — because that sets your trade frequency and therefore your `n`.

**How it lies — the two classics:**
1. **Conditioning on the future** (look-ahead). If B is knowable only after the fact, the
   probability is fake. Test: could I have evaluated B in real time, with data available then?
2. **Survivorship / selection in the conditioning set.** Computing `P(win | trade taken)` where
   "taken" was decided with hindsight inflates everything.
3. **Base-rate neglect**: a rule with `P(A|B) = 70%` but `P(B) = 2%` of days contributes almost
   nothing; a 55% rule firing daily dominates it.

**Do this:** for every rule, write both `P(A|B)` and `P(B)`, and derive expected trades per period.
`n` follows from `P(B)`, and `n` decides whether you can measure anything.

---

## 3. Expected Value

`E[X] = Σ xᵢ·P(xᵢ)`  (continuous: `∫ x·f(x) dx`)

**What it does:** the long-run average. The only quantity that adds across independent bets.

**Use in practice:** per-trade EV in the unit you actually care about:
`EV = P(win)·avg_win − P(loss)·avg_loss − costs`.
Always net of measured costs. Then `EV × trades_per_period` = expected period return, which is
what feeds sizing.

**How it lies:**
- **EV is not what you experience.** With multiplicative returns, the median path is far below the
  mean path (see §20 GBM: median grows at `μ − σ²/2`). A positive-EV bet can have a median
  outcome of ruin. Report both mean and median terminal wealth in any compounding projection.
- **EV of a truncated distribution.** Stops truncate the left tail, targets the right. Estimating
  EV from a sample where those truncations differ from the live rule is a category error.
- **EV under a changing regime is a mixture**, and the mixture weights are unknown out-of-sample.

**Do this:** report EV **and** the probability of the bad branch, **and** the median compounded
outcome. "Positive EV" alone is never a sufficient reason to size up.

---

## 4. Variance

`Var(X) = E[(X−μ)²] = E[X²] − μ²`;  `σ = √Var`

**What it does:** dispersion around the mean. The denominator of everything risk-adjusted, and
the driver of the compounding drag.

**Use in practice:**
- Per-period σ, then `SE = σ/√n` for the mean — this is the number that decides significance.
- Scaling: variance adds over independent periods, so **σ scales with √t**, not `t`.
  Annualising: `σ_annual = σ_period·√(periods per year)`.
- Compounding drag: geometric return ≈ `μ − σ²/2`. Doubling σ costs 4× the drag —
  this is why leverage that raises σ can destroy the compounded return even with unchanged mean
  (**the classic mechanism behind retracted leveraged theses**).

**How it lies:**
- **Non-stationarity.** σ estimated in calm regimes underestimates the risk you'll meet.
- **Autocorrelation.** If returns are serially correlated, `σ·√t` is wrong (use Newey-West or
  block bootstrap). Overlapping windows create this silently.
- **Fat tails.** With infinite/undefined higher moments, sample variance is unstable and
  Sharpe/VaR built on it are false precision.
- Small-sample bias: use `n−1` (Bessel). With `n = 15`, σ itself has ~±19% standard error —
  so your `t` is also uncertain.

**Do this:** never report a mean without σ and n. When comparing variants, use the **paired**
difference's σ, not each variant's σ.

---

## 5. Covariance

`Cov(X,Y) = E[(X−μₓ)(Y−μᵧ)]`

**What it does:** co-movement in raw units. Foundation of portfolio risk (§12 covariance matrix).

**Use in practice:** portfolio variance for two sleeves:
`σ²_p = w₁²σ₁² + w₂²σ₂² + 2w₁w₂·Cov(1,2)`.
This is the whole quantitative argument for a barbell: if the safe sleeve and the strategy sleeve
have low or negative covariance, total risk is much less than the weighted sum. **Verify the
covariance, don't assume it** — sleeves that look unrelated often share a funding or liquidity factor.

**How it lies:**
- **Unit-dependence** — magnitude is meaningless without scaling; use correlation (§6) to compare.
- **Covariance goes to +1 in crises.** The diversification you measured in calm data is precisely
  what disappears when you need it. Estimate a **stressed** covariance (conditional on the worst
  decile of the market) and size on that.
- **Fills decorrelate, tails cluster.** Concrete case: laddered entries at increasing distances fill on different days (so *frequency* is decorrelated), but a trend day stops **all** of them together
  (so *losses* are perfectly correlated). Covariance of returns computed over all days hides this.
  Always compute covariance **conditional on the loss tail**.

**Do this:** report both unconditional and tail-conditional covariance. Size on the latter.

---

## 6. Correlation

`ρ = Cov(X,Y)/(σₓσᵧ)`, `−1 ≤ ρ ≤ 1`

**What it does:** scale-free linear co-movement.

**Use in practice:** de-duplicating signals and sleeves. Two rungs/strategies with `ρ > 0.9` are
one position with extra fees. Diversification benefit ≈ `√((1+ρ)/2)` for two equal-weight assets:
`ρ = 0.9` gives almost nothing; `ρ = 0.3` gives a real cut.

**How it lies:**
- **Linear only.** `ρ ≈ 0` with a strong non-linear relation is routine. Plot it (Anscombe).
- **Unstable**: with `n = 15`, the 95% CI on `ρ = 0.5` is roughly `[0.0, 0.8]` — useless for
  a decision. Needed `n` for ±0.1 precision is in the hundreds.
- **Correlation of *what*?** Returns, P&L, or signed exposure give different answers.
- **Spurious via common trend.** Two rising series correlate. Difference/detrend first.

**Do this:** report `ρ` with its CI and `n`. If the CI spans 0, do not build an argument on it.

---

## 7. Standard Normal / Z-score

`Z = (X − μ)/σ`; `P(|Z| ≤ 1) ≈ 68%`, `≤2 ≈ 95%`, `≤3 ≈ 99.7%`

**What it does:** puts an observation in units of its own noise. The bridge from raw numbers to
"is this surprising?"

**Use in practice:** this *is* the significance test for a strategy result:
`Z = (mean return per period) / (σ/√n)` = the `t`-statistic (t for small n).
Also: normalising features before any optimisation, and computing VaR quantiles (§30).

**How it lies — this is the big one for trading:**
- **Financial returns are not normal.** Excess kurtosis is the norm; the 5σ events that "can't
  happen" happen. A −34% index month is ~−8σ under a normal fitted to calm data.
- Therefore: **any risk number derived from Z assuming normality understates the tail** — VaR,
  Sharpe-based probability statements, "3σ stop".
- With `n < 30`, use **t** (heavier tails) and remember σ is itself estimated.

**Do this:** use Z/t for *inference about the mean* (that's where CLT helps), but never for
*tail risk*. For tails use empirical quantiles or a fat-tailed model, and say which.

---

## 8. Law of Large Numbers

`X̄ₙ → μ` in probability as `n → ∞`  (and CLT: `X̄ₙ ~ N(μ, σ²/n)` for large n)

**What it does:** guarantees averages converge — **eventually**, and at rate `1/√n`.

**Use in practice — the sample-size calculator you must run before collecting:**
`n ≈ (t_target · σ / effect)²`

Worked example: `σ_daily ≈ 250-290pt`, effect to detect `≈ 52pt/day`, `t_target = 2`:
`n ≈ (2·270/52)² ≈ 108 sessions ≈ 5 months of weekdays`. **State that number.** It converts
"we need more data" into a decision-grade timeline, and it usually reveals that the question
being argued over is unanswerable in the available horizon.

**How it lies:**
- **`1/√n` is slow.** Halving the error needs 4× the data. 15 → 60 sessions only doubles precision.
- **Independence is required.** Overlapping windows, same-day rungs, autocorrelated regimes →
  effective `n` far below row count. Estimate effective `n` (e.g. `n/(1+2Σρₖ)`) and use *that*.
- **Stationarity is required.** If the process changes faster than `n` grows, the average converges
  to nothing meaningful. This is the fatal one for markets: more history ≠ more relevant history.
- **LLN says nothing about the path.** You can be right in the limit and bankrupt in the interim —
  which is a sizing question (Kelly/fractional), not an inference question.

**Do this:** always report `n`, **effective n**, and the `n` required for the claim.
If required `n` >> available `n`, the verdict is "not decidable" — write it and move on.
