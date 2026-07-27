# Formulas 31-37 — Validation, overfitting, and honest inference

These are the tools that separate a backtest from a discovery. The first thirty formulas describe the
*data*; none of them says anything about the **research process** that produced the result — and with
financial data that process is almost always where the return came from.

---

## 31. Walk-Forward Validation

```
rolling:   fold k trains on [k·S, k·S + W),  tests on [k·S + W, k·S + W + T)
anchored:  fold k trains on [0,   k·S + W),  tests on [k·S + W, k·S + W + T)
folds:     K = floor((N − W − T)/S) + 1     with S = T for non-overlapping tests
```

**What it does:** simulates the only thing you actually get to do in production — fit on the past,
trade the future — repeated K times. Rolling forgets old data (adapts, higher variance); anchored
keeps everything (stable, but later folds are fit on a different amount of data than earlier ones).

**Use in practice:**
- Size the **train** window by free parameters, not by calendar: at minimum ~30 independent events
  per fitted parameter, and prefer 100+. A 4-parameter rule fit on 60 trades is not fit, it is
  interpolated.
- Size the **test** window by the standard error of what you are measuring. For a Sharpe,
  `SE(SR) ≈ sqrt((1 + SR²/2)/n)`. With `n = 40` trades and `SR = 1` (per-trade units), `SE ≈ 0.19` —
  so a fold Sharpe of 0.6 and one of 1.4 are the same fold. Report the SE next to every fold number.
- With `N = 1000` days, `W = 500`, `T = 125`, `S = 125` you get `K = 4`. Four numbers. Say "four".

**How it lies:**
- **The fold count is not the evidence count.** K non-overlapping test windows inside one regime are
  one observation of one regime repeated K times, not K independent tests. Label each fold with its
  regime (vol tercile, trend sign, rate direction) and count *distinct regimes*, not folds.
- **The configurations tried across folds are invisible in the output.** If you ran a 6×8 grid inside
  every fold, the reported walk-forward equity curve is the best-of-48 path, not a live path. The
  number that matters — total configurations evaluated — appears nowhere on the chart. Log it.
- **Choosing anything by looking at the aggregated test curve makes it in-sample.** The moment you
  compare two walk-forward runs and keep the better one, the out-of-sample set has become a tuning
  set. There is no statistical repair for this, only a fresh holdout you have never queried.
- Anchored walk-forward makes early folds look bad for a boring reason (less training data), which
  invites you to discard them — reintroducing selection.
- Parameters re-fit each fold hide instability: a rule whose optimum lookback jumps 5 → 60 → 12 across
  folds has no optimum. Plot the chosen parameter per fold; a wandering optimum is a null result.

**Do this:** freeze the configuration in writing before the final fold, report a per-fold table with
`n`, SE, and regime tag, plot the selected parameter across folds, and record the total number of
configurations evaluated to feed into §33 below.

---

## 32. Purged K-Fold Cross-Validation with embargo

```
purge:    drop train obs j if  [t_{j,0}, t_{j,1}] ∩ [t_{i,0}, t_{i,1}] ≠ ∅  for any test obs i
embargo:  drop train obs j if  t_{j,0} ∈ ( max_i t_{i,1},  max_i t_{i,1} + h ]
h ≈ max(label horizon, 0.01·N)
```

**What it does:** makes cross-validation legal on time series by removing the train observations whose
information overlaps the test set. `t_{i,0}` is when observation *i*'s features are known, `t_{i,1}`
when its label resolves.

**Use in practice:**
- Standard shuffled K-fold fails on financial data for three separate reasons: **label overlap** (a
  5-day forward return labelled every day means 5 consecutive labels share 4 days of the same
  outcome), **serial correlation** (neighbouring rows are near-duplicates, so a shuffled test set is
  full of rows the model effectively saw), and **feature lookback** (a 60-day moving average at a
  test timestamp encodes 60 days that may sit in train).
- Purging removes the overlap; the embargo removes the residual leakage from features and slow
  information decay *after* the test block. With daily bars and a 10-day label, `h = 10` days
  minimum.
- Report the fraction of the training set destroyed by purging. With overlap ratio `o = h/S_fold`,
  you lose roughly `2·h` observations per fold boundary; at 10 folds on 1000 daily bars with
  `h = 10` that is ~20% of the data.

**How it lies:**
- **Without purging, accuracy is inflated almost every time, and never deflated.** The bias is
  one-directional and grows with the label horizon. An AUC of 0.52 (marginal) routinely reads 0.60+
  unpurged; a "highly significant" classifier on overlapping labels is the single most common false
  discovery in this field.
- **Purged CV fixes leakage, not selection.** It says nothing about how many feature sets, horizons,
  and hyperparameters you tried. A purged CV score is still a best-of-N score (§37).
- The performance drop after purging is often blamed on "purging being too aggressive" and then
  quietly reverted. The drop *is* the leakage you were harvesting.
- CV shuffles time, so it cannot validate anything **path-dependent** (compounding, drawdown limits,
  position carry-over, stops). Those require walk-forward (§31), not CV, regardless of purging.
- Purging assumes you know `t_{i,1}`. For labels with a variable resolution time (first touch of a
  barrier), `t_{i,1}` is itself a random variable, and using its expectation under-purges.

**Do this:** never run CV on financial labels without both purge and embargo; store `t_{i,0}` and
`t_{i,1}` per observation as columns, not as an implicit convention; and report the purged score, the
unpurged score, and the gap — the gap is your leakage measurement.

---

## 33. Deflated Sharpe Ratio (DSR)

```
DSR = Z[ (SR − SR0)·sqrt(n − 1) / sqrt(1 − g1·SR + ((g2 − 1)/4)·SR²) ]

SR0 ≈ sqrt(V) · [ (1 − γ)·Z⁻¹(1 − 1/N) + γ·Z⁻¹(1 − 1/(N·e)) ]
γ ≈ 0.5772 (Euler-Mascheroni),  e ≈ 2.71828,  V = variance of the N trial Sharpes
```

**What it does:** the probability that the true Sharpe exceeds zero *after* correcting for three
things a raw Sharpe ignores: the number of configurations searched (via `SR0`), non-normality of
returns (skew `g1`, kurtosis `g2`), and sample length (`n`). `Z[·]` is the standard normal CDF,
`Z⁻¹[·]` its inverse (`01-probability.md` §7).

**Use in practice:**
- **Units must match or the answer is nonsense.** `SR` and `SR0` are expressed at the observation
  frequency, and `n` is the number of those observations. Convert an annualised Sharpe with
  `SR_obs = SR_ann/sqrt(f)`, `f` = 252 daily, 12 monthly. Using an annualised SR with a daily `n` is
  the standard implementation bug and inflates DSR to ~1.00 for everything.
- **`SR0` is the number to internalise.** With `N = 100` trials and a spread of trial Sharpes of
  `sqrt(V) = 0.5` (annualised): `Z⁻¹(0.99) = 2.326`, `Z⁻¹(1 − 1/271.8) = 2.680`, so
  `SR0 = 0.5·(0.4228·2.326 + 0.5772·2.680) = 0.5·2.530 ≈ 1.27`. **The best of 100 variants is
  expected to show an annualised Sharpe of ~1.27 with zero real edge.** A reported 1.4 is not a
  finding; it is the median outcome of the search.
- `SR0` grows roughly with `sqrt(ln N)`, so it is brutal at first and then slow: `N = 10 → 1.58·sqrt(V)`,
  `N = 100 → 2.53·sqrt(V)`, `N = 1000 → 3.26·sqrt(V)` (same `sqrt(V)` units). Going from 10 trials to
  1000 only doubles the bar — the first handful of trials does most of the damage.
- Threshold convention: `DSR > 0.95` to call something established, `0.90-0.95` provisional,
  `< 0.90` not established. Report the value, not the verdict alone.

**How it lies:**
- **`N` is self-reported and always undercounted.** Every grid point, every discarded idea on the
  same data, every re-run after a bug fix, every "quick look" at a different exit is a trial. If you
  log 12 and actually ran 400, `SR0` is understated by ~40% and DSR is fiction. `N` must come from a
  trial log written *during* the research, never reconstructed afterwards.
- **`V` computed over the surviving variants understates the true dispersion**, because the losers
  you deleted were the wide tail. Small `V` → small `SR0` → generously high DSR. Compute `V` over
  *all* trials including the failures.
- `g1` and `g2` are estimated from the same short sample and are themselves extremely noisy —
  kurtosis needs a huge `n` to estimate at all. They sit in the denominator, so a bad kurtosis
  estimate moves DSR directly.
- Trials are correlated (a 6×8 grid of neighbouring parameters is not 48 independent draws), so
  nominal `N` overstates independence. In practice this over-corrects far less often than
  under-reporting `N` under-corrects, so it is not an excuse to shrink `N`.
- DSR assumes iid-ish returns. Autocorrelated or overlapping returns need the SE fixed first (§35),
  otherwise `sqrt(n − 1)` is counting observations you do not have.

**Do this:** keep a trial log with one row per configuration evaluated and its Sharpe; compute `V`
and `N` from that file, never from memory; report `SR`, `SR0`, `N`, `n`, `g1`, `g2`, and `DSR`
together, in observation-frequency units, or report none of them.

---

## 34. Probability of Backtest Overfitting (PBO)

```
1. build a T×N matrix M of per-period returns: one column per configuration
2. split rows into S disjoint contiguous blocks of equal length
3. for each of the C(S, S/2) combinations c: half the blocks = IS, complement = OOS
4. n*(c) = argmax_j  perf_IS(j)             ← the config the procedure would have picked
5. r(c)  = OOS rank of n*(c) among N,  ω_c = r/(N + 1),  λ_c = ln(ω_c/(1 − ω_c))
6. PBO = #{c : λ_c < 0} / C(S, S/2)         ← fraction of splits where the IS winner is
                                              below the OOS median
```

**What it does:** measures the **selection procedure**, not the strategy. It answers: if I pick the
best configuration on one half of history, how often is it a below-median performer on the other
half? `S = 16` gives `C(16,8) = 12,870` splits, enough to make the estimate stable.

**Use in practice:**
- Interpretation: `PBO < 0.10` the selection carries information; `0.10-0.25` borderline, size small
  and re-measure live; `0.25-0.50` the procedure is barely better than picking at random;
  `> 0.50` the IS winner is *systematically* an OOS loser — your search is anti-predictive and the
  honest action is to stop, not to re-tune.
- Two by-products are as useful as PBO itself: the **IS-vs-OOS scatter with its regression slope**
  (how much performance is expected to degrade — a slope near zero means IS ranking predicts
  nothing) and the **probability of loss** of the selected configuration OOS.
- Run it on the *raw* variant set from the trial log, before any human filtering.

**How it lies:**
- **PBO only sees the configurations you hand it.** Discard 300 variants by eye, feed it the
  surviving 4, and PBO comes back low because the hard selection already happened outside the
  measurement. This is the dominant failure and it produces exactly the reassuring number you wanted.
- **A degenerate variant set breaks it in both directions.** If all N columns are near-identical,
  the OOS ranks are pure noise and PBO → 0.5 even for a real edge. If you pad the set with obviously
  broken variants, the winner beats them everywhere and PBO → 0 for free. Report the pairwise
  correlation distribution of the columns alongside PBO.
- Recombining non-adjacent blocks destroys the time ordering, so returns that depend on **carry-over
  state** (open positions, compounding, vol-scaling memory) are misrepresented. PBO is valid on a
  matrix of per-period returns of *self-contained* runs, not on stitched fragments.
- `perf_IS` choice matters: ranking on Sharpe versus on total return versus on Calmar can move PBO by
  0.1 or more. Pick the metric you would actually select on, state it, and do not shop.
- PBO says nothing about whether the *best* configuration is good in absolute terms — a low PBO on a
  family of unprofitable variants is a reliable way to pick a loser.

**Do this:** run PBO on the complete trial set, report `S`, `N`, the metric used, PBO, the IS-OOS
slope, and the column-correlation spread; treat any PBO above 0.25 as a verdict of "not decidable"
rather than an invitation to reduce the variant set.

---

## 35. HAC / Newey-West standard errors

```
S    = γ0 + 2·Σ_{k=1..L} w_k·γ_k        γ_k = autocovariance at lag k
w_k  = 1 − k/(L + 1)                    (Bartlett kernel)
L    = floor(4·(T/100)^(2/9))            (rule of thumb)
SE_HAC(mean) = sqrt(S/T)                n_eff = T·γ0/S = T / (1 + 2·Σ w_k·ρ_k)
```

**What it does:** corrects the standard error of a mean or a regression coefficient for
autocorrelation and heteroskedasticity, without modelling either. Everything downstream — t-stats,
p-values, confidence bands, Sharpe significance — is wrong by the same factor if you skip it.

**Use in practice:**
- **Overlapping returns are the worst case and the most common.** Sampling 5-day returns daily gives
  `ρ_k ≈ 1 − k/5`. With `L = 4`, `ρ = (0.8, 0.6, 0.4, 0.2)` and `w = (0.8, 0.6, 0.4, 0.2)`:
  `S/γ0 = 1 + 2(0.64 + 0.36 + 0.16 + 0.04) = 3.40`. So `SE` is `sqrt(3.40) = 1.84×` the naive value,
  a naive `t = 2.6` is really `t = 1.4`, and `n_eff = T/3.4` — you have 29% of the sample you
  thought you had.
- Mild single-lag persistence still matters: `ρ1 = 0.3`, `L = 1` → `S/γ0 = 1.30`, `SE` up 14%,
  `t` down 12%. Enough to move a result across the 2.0 line.
- With `T = 1000`, the rule gives `L = floor(4·10^(2/9)) = floor(4·1.67) = 6`.
- Apply to Sharpe significance too: the Sharpe t-stat is `SR·sqrt(n)` only under iid; with
  autocorrelation use `SR·sqrt(n_eff)`.

**How it lies:**
- **In small samples HAC standard errors are biased downward** — still too optimistic, just less so.
  A HAC-corrected `t = 2.1` on 200 observations is not a comfortable result.
- **`L` is a free parameter and the answer depends on it.** Too small under-corrects; too large makes
  `S` noisy and can drive it toward zero or negative in finite samples. Report the `t` at `L`, `2L`,
  and `L/2`; if the conclusion changes, you have no conclusion.
- HAC fixes the *uncertainty*, never the *point estimate*. A mean return biased by survivorship,
  look-ahead, or mid-price fills stays biased with a wider band around it.
- It does nothing about selection. A HAC-corrected p-value on the best of 200 variants is still
  meaningless (§37) and a HAC-corrected Sharpe still needs deflating (§33).
- **Autocorrelation from overlapping windows is a design flaw, not a statistical inconvenience.**
  Patching the SE keeps a sample whose effective size you have shrunk by 3× or more. Non-overlapping
  sampling costs you rows and buys you real degrees of freedom.
- Non-stationarity masquerades as autocorrelation. If `γ_k` decays slowly and never dies, the process
  may have a changing mean — in which case no SE is meaningful and the correct answer is a regime
  split.

**Do this:** report every t-stat with its HAC lag `L`, the variance-inflation factor `S/γ0`, and
`n_eff` beside the nominal `n`; and prefer redesigning to non-overlapping observations over
correcting for overlap.

---

## 36. Minimum Backtest Length / Minimum Track Record Length

```
MinTRL = 1 + [ 1 − g1·SR + ((g2 − 1)/4)·SR² ] · ( Z_α / (SR − SR*) )²
```

**What it does:** the number of observations required before an observed Sharpe `SR` can be
distinguished from a benchmark `SR*` at confidence `α`. It is the inverse of the Sharpe significance
test, and it is the honest way to answer "is this track record long enough yet?".

**Use in practice:**
- Units again: `SR` at the observation frequency, `MinTRL` in that same unit. `Z_0.95 = 1.645`.
- **Monthly, `SR_ann = 1.0`, normal returns.** `SR_m = 1/sqrt(12) = 0.2887`, `g1 = 0`, `g2 = 3`:
  bracket `= 1 + 0.5·0.0833 = 1.042`; `(1.645/0.2887)² = 32.5`;
  `MinTRL = 1 + 1.042·32.5 ≈ 35 months ≈ 2.9 years`. An annual Sharpe of 1 needs about three years of
  monthly data to be distinguishable from zero. Most proposals have a fraction of that.
- **Same Sharpe, ugly distribution** (`g1 = −1.5`, `g2 = 8`): bracket `= 1 + 0.433 + 0.146 = 1.579`;
  `MinTRL ≈ 52 months ≈ 4.4 years`. Negative skew and fat tails cost you eighteen extra months for
  the same headline number.
- **Daily, `SR_ann = 2.0`** (a plausible intraday claim): `SR_d = 0.126`, `(1.645/0.126)² = 170`;
  normal → `~173 days ≈ 8 months`; with `g1 = −2`, `g2 = 12` → bracket `1.296` → `~222 days ≈ 11
  months`. This is the useful direction of the formula: it tells you the *minimum* observation window
  before a live result means anything, which is usually longer than the pilot you planned.
- Sensitivity is quadratic in the gap: halving `SR − SR*` quadruples `MinTRL`. Benchmarking against
  a risk-free rate instead of zero can double the requirement.

**How it lies:**
- **It assumes stationarity** — that the process generating the observed `SR` continues unchanged for
  the whole `MinTRL`. If the regime turns over faster than `MinTRL`, the requirement can never be
  met and the correct verdict is permanent "not decidable", not "wait longer". Compare `MinTRL` to
  your estimated regime length; if `MinTRL > ` regime length, say so explicitly.
- **`SR* = 0` is almost always the wrong benchmark for a searched strategy.** If `SR` is the best of
  `N` variants, the null is not zero but `SR0` from §33. Substituting `SR* = SR0` typically multiplies
  `MinTRL` by 3-10×. Using zero is how a 3-year requirement is reported as 8 months.
- `g1` and `g2` come from the same short sample whose adequacy you are testing — circular, and both
  are badly estimated at small `n`. Sensitivity-test with the worst plausible skew/kurtosis, not the
  measured one.
- It is a **minimum, necessary, not sufficient**: reaching `MinTRL` clears one hurdle (sample length)
  and none of the others (leakage §32, selection §37, costs `06-sizing-execution.md` §43).
- The formula rewards a strategy for having a high in-sample Sharpe by *shortening* the required
  history — precisely the wrong incentive if that Sharpe is itself the product of overfitting.

**Do this:** compute `MinTRL` with `SR* = SR0` and with pessimistic `g1`/`g2`, compare it to both the
available history and the plausible regime length, and state the ratio. If you have 20% of `MinTRL`,
the deliverable is a hypothesis and a measurement plan, not a Sharpe.

---

## 37. Multiple Testing & False Discovery Rate

```
Bonferroni (FWER):  reject H_i if  p_i ≤ α/m
Holm (FWER, stronger): sort p ascending; reject while  p_(i) ≤ α/(m − i + 1)
Benjamini-Hochberg (FDR): sort p ascending; find largest i with  p_(i) ≤ (i/m)·q;
                          reject H_(1)…H_(i)
Benjamini-Yekutieli (arbitrary dependence): replace q by  q / Σ_{k=1..m} (1/k) ≈ q/(ln m + 0.577)
```

**What it does:** rescales the significance bar for the fact that you tested `m` hypotheses.
Bonferroni/Holm control the probability of **any** false positive (FWER); BH controls the **expected
proportion** of false positives among the rejections (FDR) and is much more powerful when `m` is
large.

**Use in practice:**
- **The bar moves a long way.** With `m = 200` variants and `α = 0.05`, Bonferroni requires
  `p ≤ 0.00025`, i.e. two-sided `|t| ≥ 3.66`. The habitual `t > 2` bar carries a ~64% chance of at
  least one false positive at just `m = 20` (`1 − 0.95²⁰`), and near-certainty at `m = 200`. This is why a serious bar for
  a newly claimed effect is `t > 3`, not `t > 2`.
- BH at `q = 0.10` with `m = 200`: the smallest p-value needs `≤ 0.0005`, the tenth needs
  `≤ 0.005` — far more permissive than Bonferroni's flat `0.00025`, which is why BH is the right
  default when screening many candidate signals rather than defending one.
- **Search headroom.** Treat your test budget as a resource: `m ≈ (product of grid sizes) × (number
  of re-runs) × (number of variants examined by hand)`. Decide `m` and the bar *before* searching.
  Once spent, the only way to buy new headroom is new data — a fresh instrument, a fresh period, a
  live period — not a new metric on the same data.
- This is the same quantity as `N` in §33 and the same discipline as the base-rate calculation in
  `01-probability.md` §1. Consistency check: if you report `m = 5` here and `N = 300` there, one of
  them is a lie.

**How it lies:**
- **`m` is self-reported and structurally undercounted.** The variants you abandoned after a glance,
  the "data checks" that were really tests, the parameter you nudged three times, the rerun after
  fixing a bug — all are tests. Undercounting `m` by 10× moves the Bonferroni bar by roughly `0.6` in
  t-units, which is exactly the range where marginal results live.
- **The p-values being corrected are usually already wrong.** Autocorrelated or overlapping data
  inflate every nominal `t` (§35). Correcting inflated p-values for multiplicity gives you a
  precise-looking number built on two errors.
- **Tests are correlated,** so Bonferroni is genuinely conservative and BH's guarantee needs a
  positive-dependence condition. This is a real objection and also the most abused one: it is used to
  justify no correction at all. Use BY if you are worried about dependence; never use dependence as
  a reason to keep `m = 1`.
- **FDR is a portfolio-level statement, not a per-decision one.** "10% of my discoveries are false"
  does not tell you whether *this* one, the one you are about to fund, is real. For a single deploy
  decision, the posterior in `01-probability.md` §1 is the right object.
- Passing a multiplicity correction says nothing about economics: a survivor at `t = 4` with a gross
  edge below transaction costs (`06-sizing-execution.md` §43) is statistically real and financially
  worthless.

**Do this:** write `m` down before you start and log every test as you run it; report the corrected
threshold in t-units next to every claim; and when `m` is genuinely unknown, say the result is not
decidable rather than quoting a nominal p-value.

---

Next: `06-sizing-execution.md` — how much of a validated edge survives sizing, capacity, and fills.
