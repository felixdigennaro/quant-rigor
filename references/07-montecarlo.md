# Formulas 45-49 — Monte Carlo: building the null instead of assuming it

A t-statistic assumes a distribution. These build the distribution from your own data, which is the
only honest way to answer "how often does no-edge produce a result this good?" — and the only way to
handle statistics with no closed form (max drawdown, hit rate, path-dependent metrics).

---

## 45. Permutation and sign-flip tests

```
permutation:  reorder the periods            -> null keeps the marginal distribution, destroys timing
sign flip:    r_t -> ±r_t with p = 1/2       -> null keeps magnitudes, destroys direction
p_hat = (1 + #{null >= observed}) / (1 + reps)
```

**What it does:** constructs an exact finite-sample null by re-randomising the one structure the
claim depends on. Distribution-free: no normality, no independence assumption beyond the shuffle.

**Use in practice:** choose the shuffle that matches the claim.
- The claim is about the **mean** (does the strategy make money?) → **sign flip**. Reordering cannot
  change a sum, so a permutation test on the mean is vacuous — a mistake that appears in real papers.
- The claim is about **timing/sequence** (trend capture, stop placement, compounding, drawdown) →
  **permutation** of the order, with `sharpe` or `max_dd` as the statistic.
- Autocorrelated data → shuffle **blocks**, not single periods: a stricter and more honest null.

The `+1` in the numerator matters: it prevents reporting `p = 0`, which would claim infinite evidence
from a finite number of draws.

**How it lies:** a permutation test answers only the question its shuffle encodes. Shuffling order and
then concluding "the mean is significant" is a category error. It also assumes exchangeability under
the null — with strong regime structure, plain shuffling produces a null that is *too easy* to beat,
inflating significance. And it says nothing about out-of-sample: it is a statement about this sample.

**Do this:** state the shuffle, the statistic, and therefore the exact null being rejected. Use
sign-flip for the mean, block-permutation for anything sequential, and report `reps`.

---

## 46. Surrogate paths (bootstrap and parametric simulation)

```
block bootstrap: resample contiguous blocks of length b  -> keeps fat tails + local autocorrelation
iid bootstrap:   resample single periods                 -> destroys autocorrelation
GBM null:        r_t ~ N(0, sigma_hat)                   -> parametric, zero drift
OU null:         dx = -theta*x*dt + sigma*dW             -> mean-reverting null
```

**What it does:** simulates many histories from a process with **no edge**, calibrated to the
observed second moment, and reports where the real result falls in that distribution.

**Use in practice:** remove the observed mean (that is what makes it a no-edge null), keep the
volatility and — with block bootstrap — the tail shape and serial dependence. Report the observed
statistic against the null's median, 95th and 99th percentiles. For a fade strategy the OU null is
the informative one: it asks whether the result exceeds what mean reversion alone would give.
Block length should be at least the autocorrelation horizon; too short a block manufactures an
easy null.

**How it lies:** the null inherits your calibration. A GBM null has thin tails, so a strategy that
harvests fat tails looks miraculous against it. A bootstrap null cannot contain any regime the sample
does not contain. And `drift_zero=False` silently turns a no-edge test into a self-fulfilling one.

**Do this:** run at least two nulls of different character (block bootstrap and a parametric one).
If the verdict depends on which null you choose, that dependence *is* the result — report it.

---

## 47. White's Reality Check

```
V = max_k sqrt(T) * mean(d_k)            d_k,t = strategy k performance - benchmark, per period
V*_b = max_k sqrt(T) * (mean(d*_k) - mean(d_k))     re-centred bootstrap draw b
p = #{V*_b >= V} / B
```

**What it does:** tests `H0: the BEST of k screened strategies does not beat the benchmark`. This is
the correct test after a search; a per-strategy p-value is not, because the maximum of k statistics
has a different distribution from a single one.

**Use in practice:** feed one series per candidate of *excess* performance versus the benchmark, in
matched periods. Use a block/stationary bootstrap over the period index so all candidates are
resampled together — that preserves the cross-sectional dependence which is exactly what makes the
maximum non-trivial.

**How it lies:** it is **conservative**, and gets worse as you add poor candidates: a single hopeless
strategy inflates the null's maximum and can hide a genuinely good one. It also assumes the candidate
set is fixed in advance — adding strategies after seeing results breaks it.

**Do this:** report `k`, the identity of the winner, and its excess mean alongside the p-value. If
the set includes obviously weak candidates, use SPA (§48) as well and report both.

---

## 48. Hansen's SPA test

```
T_spa = max(0, max_k sqrt(T) * mean(d_k)/sd(d_k))        studentised
retain only k with sqrt(T)*mean(d_k)/sd(d_k) >= -sqrt(2 log log T)
p from the re-centred bootstrap of the retained set
```

**What it does:** the Reality Check's refinement. Studentising puts strategies with different
volatilities on the same scale, and the retention threshold drops candidates too poor to influence
the null — removing the RC's main weakness.

**Use in practice:** the default test when the candidate set is large or heterogeneous. Report both
SPA and RC: SPA more powerful, RC more conservative. If they disagree, the honest reading is that the
result depends on how the comparison set was assembled.

**How it lies:** the threshold is asymptotic and arbitrary at small `T`; results can move with it.
Studentisation punishes low-volatility strategies less, which can flip the ranking versus RC. Neither
test protects against a candidate set built by looking at the data first.

**Do this:** pre-register the candidate set. Report `k`, how many were retained, both p-values, and
the threshold used.

---

## 49. Moment stress (parameter uncertainty)

```
for each adversarial perturbation: recompute the statistic and its t
  sigma understated by 30%        costs doubled          best period removed
  two best periods removed        mean shrunk 25%        one extra tail event
```

**What it does:** asks whether the conclusion survives the estimation error in its own inputs. Every
number in a backtest is an estimate; the verdict should not hinge on the third decimal of one of them.

**Use in practice:** run it before adopting anything. A conclusion that survives 6/6 perturbations is
robust; one that survives 1/6 is a point estimate wearing a verdict's clothes. The "best period
removed" rows double as a concentration check, and "costs doubled" is the cheapest test of whether an
edge is really a fee rebate.

**How it lies:** the scenario list is a choice, and a flattering list produces a flattering answer.
It also perturbs one dimension at a time, while real stress arrives jointly (vol up, liquidity down,
costs up, correlations to one).

**Do this:** publish the scenario table in full — including the rows where the conclusion fails —
and add at least one *joint* scenario that moves several inputs the wrong way at once.

---

Next: `references/backtest-protocol.md` for the checklist these tests plug into, and
`references/contrarian-protocol.md` for the attacks they are meant to survive.
