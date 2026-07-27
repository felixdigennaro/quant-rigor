# Charts — the five plots, and what each one is for

Tables average away shape. Every one of these plots exists to expose a specific failure that a
table hides. Produce them with `scripts/qcharts.py`; if a chart is going to be **published**
(site, investor doc), load the `dataviz` skill for styling first.

Rules for all charts: label units, state `n` and the period, mark costs as included/excluded,
and never plot a projection as a single smooth line.

---

## 1. Equity curve — baseline vs variant(s), same axes

**Exposes:** whether an improvement is a **slope** (persistent edge) or a **step** (one lucky period).
A variant that tracks the baseline and then jumps once is a single observation dressed as a strategy.

Read it for: monotonic vs step gains · when the curves diverge (and whether that date coincides
with a regime change) · flat stretches (edge absent, costs still paid).

---

## 2. Underwater / drawdown curve

`DD_t = equity_t / running_max(equity) − 1`

**Exposes:** depth **and duration**. Two strategies with equal return and equal max drawdown are
not equivalent if one recovers in a week and the other takes three months — the second is the one
you abandon at the worst moment.

Always report alongside: max DD, **time under water** (longest and median), and current DD.
For a leveraged or margin-financed strategy, mark the level at which a margin call/liquidation
would have occurred — a drawdown that crosses it is not a drawdown, it's a termination.

---

## 3. Per-period contribution bars (sorted)

One bar per period, sorted by magnitude; for a comparison, plot the **paired difference**.

**Exposes:** concentration. Instantly answers "does this live in 3 days?" and shows how many
periods are improved vs worsened. Annotate: top-3 share of total, and count of
improved/worsened/unchanged.

This is the chart that killed several plausible findings in this repo. It is the cheapest and most
informative of the five — produce it first.

---

## 4. Parameter surface (metric vs parameter, 1-D or 2-D heatmap)

**Exposes:** overfitting, visually and unambiguously.
- **Sharp isolated peak** → fit to noise. Do not adopt.
- **Broad plateau** → robust; adopt the **middle**, never the edge.
- **Monotone ramp** → the metric is tracking something else (often a hidden risk or resolution
  artefact); investigate before believing.

Overlay the **search headroom** (`σ_period·√n·√(2 ln N)`) as a band: any variation inside it is
noise, and drawing that band on the surface usually ends the parameter debate immediately.

---

## 5. Return distribution + QQ vs normal

Histogram of per-period returns with the normal fit overlaid, plus a QQ plot.

**Exposes:** whether Sharpe and VaR are lying (both assume near-normality —
`04-finance.md` §28, §30). Look for negative skew and fat left tail: the signature of a
**short-gamma** strategy (§27). If the strategy is structurally short gamma but the sample shows
positive skew, the tail is **missing from the sample**, not absent from the world — flag it.

Annotate: skew, excess kurtosis, empirical 5% VaR, CVaR, worst period.

---

## Optional but high-value

- **Rolling metric** (return, Sharpe, or hit-rate over a moving window): shows decay. An edge that
  shrinks monotonically through the sample is being arbitraged, or was never there.
- **Trade-level scatter**: MAE (max adverse excursion) vs final P&L per trade. This is the honest
  way to reason about stops: it shows how much heat winners took, i.e. exactly what a tighter stop
  would have destroyed. **Use this instead of a trailing-stop backtest when resolution is
  insufficient** — MAE/MFE distributions are far more resolution-tolerant than a path-dependent
  simulation.
- **Scenario bar chart**: P&L under the 5-10 named scenarios (venue failure, gap, funding spike,
  index −10%). No distributional assumption; covers the risks VaR cannot see.
