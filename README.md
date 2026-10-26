# quant-rigor

**A protocol that stops a plausible-but-wrong backtest from becoming a decision that moves money.**

Not a formula sheet. Formula sheets are free and everywhere. This is the operating procedure that
decides *whether you know anything at all* — and it is built around the failure modes that make
backtests lie: unmeasurable rules, contractions that win by algebra, best-of-N winners inside the
noise, autocorrelation-inflated t-stats, and costs quoted for the wrong variant.

Works as a Claude Code / Claude Agent skill, and the Python side runs standalone.

---

## Quickstart (30 seconds, no dependencies)

```bash
python3 tests/test_all.py            # checks against independent references; exits 1 if any fails
python3 examples/demo_strategy.py    # full protocol on one strategy + writes an HTML report
python3 examples/demo_portfolio.py   # portfolio risk: contributions, tail dependence, PCA
```

The demo needs no data, no network and no packages. It generates a reproducible market, runs a
strategy plus a trailing-stop variant through all eight steps, and prints the four traps being
caught — every number computed live.

Analyse your own series immediately:

```bash
# one variant
python3 scripts/qstats.py --series 410 -322 234 97 -46 233

# baseline vs variant, paired, with the best-of-N correction for 12 configurations searched
python3 scripts/qstats.py --a 410 -322 234 --b 710 -322 434 --variants 12

# can my data even measure a 40-point trailing stop? (pass your bar ranges)
python3 scripts/qstats.py --ranges 67 40 108 72 55 --dist 40

# validation: HAC t, Deflated Sharpe, minimum track record length
python3 scripts/qvalidate.py --series 0.01 -0.02 0.03 0.01 --trials 40

# sizing: Kelly ceiling vs drawdown-constrained fraction
python3 scripts/qsize.py --mean 0.004 --sd 0.02 --dd 0.15

# portfolio: risk contributions, tail dependence, PC1 concentration, condition number
python3 scripts/qportfolio.py --json '{"a":[0.01,-0.02,0.03],"b":[0.00,-0.01,0.02]}' --cap 0.4

# monte carlo: how often does no-edge produce this? (sign-flip, permutation, surrogates, stress)
python3 scripts/qmontecarlo.py --series 0.01 -0.02 0.03 0.01 -0.01 0.02 --reps 4000

# what history can I actually get for a symbol, right now?
python3 scripts/qdata.py --symbol NQ=F
```

---

## What's inside

| Path | What it is |
|---|---|
| `SKILL.md` | the workflow: 8 steps, in a fixed order, plus when it fires |
| `references/01-probability.md` | formulas 1-8 — inference, what a sample can support |
| `references/02-linear-algebra.md` | 9-15 — portfolios, risk decomposition, optimisation |
| `references/03-stochastic.md` | 16-22 — price dynamics; §21 is the resolution limit |
| `references/04-finance.md` | 23-30 — discounting, options, Sharpe, CAPM, VaR |
| `references/05-validation.md` | 31-37 — walk-forward, purged CV, DSR, PBO, HAC, MinTRL, FDR |
| `references/06-sizing-execution.md` | 38-44 — Kelly, drawdown sizing, impact, capacity, fills |
| `references/07-montecarlo.md` | 45-49 — permutation, surrogate paths, Reality Check, SPA, moment stress |
| `references/routing-table.md` | **which formula answers which question** — start here |
| `references/backtest-protocol.md` | the checklist, bias by bias |
| `references/contrarian-protocol.md` | 8 attacks whose mandate is to refute |
| `references/charts.md` | the five plots and what each one exposes |
| `references/worked-example.md` | a full case, end to end, with the four traps it catches |
| `scripts/qstats.py` | SE, t, bootstrap, jackknife, headroom, required n, measurability |
| `scripts/qvalidate.py` | HAC, Deflated Sharpe, PBO, purged K-fold, MinTRL, BH-FDR |
| `scripts/qsize.py` | Kelly, drawdown-constrained f, vol targeting, impact, capacity |
| `scripts/qportfolio.py` | covariance, eigen/PCA, risk contributions, risk parity, shrinkage, caps |
| `scripts/qcharts.py` | five charts as standalone SVG |
| `scripts/qmontecarlo.py` | permutation, sign-flip, surrogate paths, White RC, Hansen SPA, moment stress |
| `scripts/qdata.py` | the only networked module: public price history + your own CSV files |
| `scripts/qreport.py` | assembles everything into one self-contained HTML report |
| `tests/test_all.py` | assertions vs **independent** references; prints the passed and failed counts, exits 1 if any fails |
| `examples/demo_strategy.py` | runnable end-to-end demo (single strategy) |
| `examples/demo_portfolio.py` | runnable demo (multi-sleeve book, formulas 9-15) |

Every formula entry has the same four parts: the formula, **how it's used**, **how it lies**, and
the concrete action. The third part is the one you're paying for.

---

## Design decisions worth knowing

**Stdlib only, except `qconfluence` and `qregime`, which import numpy.** No scipy, pandas or
matplotlib; the other scripts, the demos and the test suite run without numpy. It runs inside
sandboxes, CI, restricted environments and agent runtimes where installing packages isn't an
option. The t-distribution, inverse normal, bootstrap and all charts are implemented from scratch.
Charts are emitted as self-contained SVG.

**Validated against independent references, not against itself.** The test suite checks published
Student-t critical values, known normal quantiles, closed-form results and analytic identities.
That discipline caught two real bugs during development: a 1% p-value error in the t-tail at df=1
(heavy Cauchy tails, truncated integration), and a `math.gamma` overflow above df≈340 that made any
multi-year daily series crash outright. Both are now regression-tested.

**Reports state what the evidence does not support.** `qreport.py` emits a self-contained HTML
document — assumptions, measurability, full metric set with standard errors, overfitting
diagnostics, sizing constraints, five charts, verdict, and an explicit "what this does not support"
section. No persuasive language, no omitted caveats, no external dependencies (publishable as an
artifact or printable to PDF).

**Opinionated defaults.** The protocol has a fixed order and a default action: *change nothing*
unless the evidence is decidable-positive. The burden of proof sits on the modification, never on
the status quo.

---

## Limitations — read before relying on it

- **Three domains exercised, not all.** Validated end to end on (a) intraday point-based strategies
  (60 sessions, ladder entries, fixed and path-dependent exits), (b) daily percentage returns with
  volatility clustering (756 observations, walk-forward, purged K-fold, HAC, sizing, charts) and
  (c) a four-sleeve portfolio (covariance, eigendecomposition, risk contributions, risk parity,
  constrained weights, tail dependence). Options books and intraday market-microstructure work use
  the same primitives but have not been exercised end to end — expect edge cases there.
- **The report generator produces a document, not a decision.** It refuses to hide caveats by
  design, which means its output is less flattering than most backtest reports.
- **It cannot manufacture data.** If your sample is too short or too coarse, the honest output is
  "not decidable" plus what would make it decidable. That is a feature and it will sometimes be
  the answer you didn't want.
- **No live trading and no execution.** Analysis only.
- **Data is the binding constraint, not the analysis.** `qdata.py` loads free public history, and
  what that yields is measured, not assumed: roughly **20 years of daily** bars, **2 years hourly**,
  but only **~60 calendar days of 5-minute** (~50 sessions) and **~7 days of 1-minute** data. Those
  intraday limits are hard — older intraday windows return empty. So a **daily** strategy can be
  evaluated on a sample long enough to matter, while an **intraday** strategy cannot be, from free
  sources: supply your own history via `from_csv`, or collect forward with `append_csv`.
  `qstats.required_n()` will usually quantify how far short the free sample falls; the protocol
  prints that instead of running the backtest anyway.
- **PBO and DSR need enough configurations and periods** to mean anything; with a handful of both,
  they return wide or undefined answers by design rather than false comfort.
- **Not a substitute for out-of-sample.** Nothing here rescues a strategy that has never been
  tested on data it wasn't built on.

---

## Disclaimer

This software and its documentation are provided for **research and educational purposes only**.
They do **not** constitute investment advice, a recommendation, or an offer to buy or sell any
financial instrument. Nothing here predicts returns. Trading and investing carry the risk of
substantial loss, including total loss of capital. Backtested and simulated results are
hypothetical, do not reflect actual trading, and are not indicative of future results. You are
solely responsible for your own decisions and for compliance with the laws and regulations that
apply to you. The authors accept no liability for any loss arising from use of this material.

## License

See `LICENSE`.
