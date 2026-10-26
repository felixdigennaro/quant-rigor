# Microstructure, costs and mean reversion — formulas 50-59

Added in v1.4 after studying four public repositories. This file states **what was taken,
from where, and what it is for** — provenance matters, because a formula whose origin you
cannot name is a formula whose assumptions you do not know.

## Provenance

| repo | what was taken |
|---|---|
| [stefan-jansen/machine-learning-for-trading](https://github.com/stefan-jansen/machine-learning-for-trading) | ch.18 cost taxonomy, spread estimation (Roll + Corwin-Schultz), the *cost cliff*, commission-vs-slippage decomposition; ch.16 §16.7 RAS protocol |
| [cantaro86/Financial-Models-Numerical-Methods](https://github.com/cantaro86/Financial-Models-Numerical-Methods) | Ornstein-Uhlenbeck fitting, Kalman-filter regression for a time-varying hedge ratio |
| [wilsonfreitas/awesome-quant](https://github.com/wilsonfreitas/awesome-quant) | data-source routing (below) |
| [LechGrzelak/Computational-Finance-Course](https://github.com/LechGrzelak/Computational-Finance-Course) | consulted; derivatives-pricing focus, nothing adopted — noted here so the omission is deliberate rather than an oversight |

Clones live in `~/quant-research-repos/`.

---

## 1. The gap this closes

Before v1.4 this package could measure **uncertainty** but not **cost**. That is a serious
hole, because an assumed cost model silently decides verdicts:

- a rule with a 40pt target and a 7pt round-trip is an 18% cost ratio;
- the same rule at 0.7pt is 1.8%;
- these are different businesses, and no downstream statistics can tell them apart.

`qmicro.py` estimates cost **from the price data you already have**, so the assumption becomes
a measurement. Use it before quoting any net figure.

---

## 2. Spread estimation (50-51)

### Roll (1984) — `qmicro.roll_spread`
With trades alternating bid/ask, consecutive price changes have covariance `-(s/2)²`, so
`s = 2·√(-cov)`.

**Guards that matter:** returns `None` when the covariance is non-negative (trending data —
the estimator is undefined), and also when the implied spread is a negligible fraction of
price. That second guard is not pedantry: floating-point rounding on a deterministic series
produces *alternating* differences, hence a strongly negative covariance and a spread of
order 1e-15. A physical guard catches what a sign test cannot.

### Corwin-Schultz (2012) — `qmicro.corwin_schultz_spread`
Uses the high-low range instead of serial covariance. Over a short interval the high is
predominantly at the ask and the low at the bid; volatility scales with time while the spread
does not, so comparing one-period and two-period ranges separates them.

```
beta  = ln(H_t/L_t)² + ln(H_{t-1}/L_{t-1})²
gamma = ln(max(H_t,H_{t-1}) / min(L_t,L_{t-1}))²
alpha = (√(2β) − √β)/(3 − 2√2) − √(γ/(3 − 2√2))
S     = 2(e^α − 1)/(1 + e^α)
```

**Prefer it on OHLC bars**: it does not break on trending data. Report *both* — when the two
disagree materially, that disagreement is information about the data, not a nuisance.

### `qmicro.bounce_check` — the test that settles short-horizon reversion
Measured negative autocorrelation at 1-5 minutes is nearly always the bid-ask bounce. The
test: estimate the implied spread from the same data and compare it to the real tick.

> If the implied spread is of the order of one tick, the reversion you measured **is** the
> cost you would pay trying to harvest it.

This is the single highest-yield check before building any short-horizon reversion system,
and it costs one function call.

---

## 3. The cost cliff (53)

`qmicro.cost_cliff`, `max_trade_frequency`, `required_edge_per_trade`.

Per-trade friction is fixed while the predictable component shrinks with horizon, so the
**same alpha** is valuable daily and worthless intraday. `max_trade_frequency` inverts the
relation: given a gross Sharpe, a cost and a target net Sharpe, it returns the highest
trade frequency that still clears the target.

**Run it before building an intraday system.** If the ceiling is below the frequency your
idea needs, the idea is dead and no backtest is required.

Rule of thumb worth saying out loud: if round-trip cost exceeds ~20% of the target, the
strategy is an execution problem wearing a strategy costume.

---

## 4. Where did the money go (54)

`qmicro.decompose_slippage` splits the gap between theoretical and realised P&L into
**entry slippage · exit slippage · commission · unexplained residual**.

The three components have three different fixes — order type, exit design, venue — and
confusing them wastes months. Run this the *first* time realised P&L disagrees with the
simulator, not the tenth.

The residual must close to zero. If it does not, the model of the trade is wrong, and that
is more urgent than the size of the gap.

---

## 5. Rademacher Anti-Serum (38, in `qvalidate.py`)

`rademacher_complexity`, `massart_bound`, `effective_trials`, `ras_sharpe_adjustment`.

### Why this beats the Deflated Sharpe for a real search
DSR corrects for the **number** of trials assuming they are independent. A real parameter
sweep produces highly **correlated** variants (lookback 11 vs 12, threshold 0.95 vs 1.00).
Over a correlated family the supremum is far smaller than over an independent one, so DSR
over-penalises and can reject an edge that would survive.

Empirical Rademacher complexity measures the family you *actually* searched, whatever its
correlation structure, and is distribution-free:

```
R̂ = E_σ[ max_i (1/T) Σ_t σ_t x_t^(i) ],   σ_t ∈ {−1,+1}
Massart bound:  R̂ ≤ √(2 ln N / T)        (attained only when the N are uncorrelated)
RAS floor:      θ_true ≥ θ̂ − 2R̂ − √(ln(1/δ)/(2T))
```

Measured behaviour on synthetic families of 40 variants over T=1000:

| family | R̂ | ratio to Massart | effective trials |
|---|---|---|---|
| correlated (shared base + noise) | 0.016 | 0.19 | **1.1 of 40** |
| independent | 0.068 | 0.80 | 10.4 of 40 |

`effective_trials` inverts Massart and is the number to **report instead of the raw count**.
Saying "we tried 217 configurations" overstates a search whose variants were neighbours on
one grid.

**Failure mode:** RAS bounds the optimism of the *search*, not the truth of any single
strategy. A small R̂ does not make a bad strategy good.

**Use it alongside DSR, never instead of it.** They answer slightly different questions, and
when they disagree that disagreement is the finding.

---

## 6. Mean reversion and pairs (55-59, `qpairs.py`)

Order of operations for any "this spread comes back" thesis:

1. **Does it revert?** — `hurst`, `adf_stat`
2. **How fast?** — `ou_half_life` ← decides the economics
3. **Is the relationship stable?** — `engle_granger`, `kalman_hedge_ratio`

### Half-life comes before the backtest
Half-life sets the holding period → the trade count → the cost hurdle. A 3-day half-life on
daily bars gives ~80 round trips a year; a 200-day half-life gives ~1 and can never pay for
a research programme. `pair_report` turns half-life straight into estimated trades/year and
annual cost, which is the number that kills most pairs before any simulation.

### Two traps encoded in the code
- **Engle-Granger critical values are more negative than plain ADF** (5%: −3.34 vs −2.86),
  because the residual comes from an *estimated* regression. Using the plain values is a
  classic way to "find" cointegration that is not there. Both tables are in `qpairs.py`.
- **`adf_stat` deliberately returns no p-value.** The null distribution is non-standard;
  a fabricated p-value would be exactly the false precision this package exists to prevent.
  Compare to the printed critical values.

### Hurst: definition matters
Uses the **second moment about zero**, per `E[|X(t+τ)−X(t)|²] ∝ τ^(2H)`. Using the variance
about the mean silently removes the drift and returns H≈0 for a pure trend — a real bug
caught by the test suite. A strongly drifting series returns H≈1, which is correct; de-trend
first if you wanted the Hurst of the fluctuations.

Validated: random walk 0.47 · OU(φ=0.9) 0.15 · trend 0.77.

### Kalman hedge ratio
A static OLS hedge ratio assumes the relationship never changes; it always does. But a filter
that adapts fast enough will track **any** relationship, including one that has genuinely
broken, and will keep promising reversion that never comes. Report the **path** of beta:
a beta wandering monotonically is a relationship dying, not a spread oscillating. The
function emits a warning when drift dominates range.

---

## 7. Testing statistical procedures — a lesson worth keeping

Four tests added in v1.4 failed at first. **The estimators were correct; the tests were wrong.**
They asserted a property of a *single realisation* of a 5%-level test, which fails 5% of the
time by construction.

The fix — and the rule going forward:

> To test a statistical procedure, test its **operating characteristics** over many
> replications (false-positive rate, median of the estimator), never a single draw.

Measured after the rewrite: ADF false-positive rate on random walks ≈3.3% at the 5%
threshold; half-life median over 200 replications 6.79 against a closed form of 6.58.
Those are the assertions now in the suite.

---

## 8. Data-source routing (from awesome-quant)

Extends `qdata.py`, which only speaks to Yahoo. **All free, no API key**, unless noted:

| need | source |
|---|---|
| **futures positioning (COT)** | `commitments-of-traders` — CFTC Legacy/Disaggregated/TFF via the official Socrata API. Weekly, decades of history. The most useful free dataset for a futures trader that this package did not know about. |
| US reference rates | `newyorkfed` — SOFR/EFFR/OBFR, SOFR averages, SOMA holdings |
| Treasury / debt | `treasurydirect`, `treasury-fiscaldata` — auctions, CUSIPs, debt, exchange rates |
| macro + multi-source | `pandas-datareader` (FRED, OECD, Fama/French, World Bank, Eurostat) |
| crypto OHLCV, multi-venue | `pricehub` (Binance, Bybit, Coinbase, OKX, Kraken, KuCoin, Bitget) |
| broad research terminal | OpenBB Terminal |
| CME intraday, paid | Databento (`qdata.db_cost` / `db_download` already wired) |

**Retention reality is unchanged**: free intraday history remains the binding constraint
(~60 days at 5m, ~7 days at 1m from Yahoo). None of the above fixes it. Only a vendor
(Databento) or a forward collector does.
