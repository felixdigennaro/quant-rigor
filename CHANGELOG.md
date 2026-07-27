# Changelog

Versioning: MAJOR.MINOR.PATCH. MAJOR changes the protocol's steps or verdict semantics; MINOR adds
formulas, modules or functions; PATCH fixes correctness or presentation.

## 1.3.0 — 2026-07-26

### Added to `qdata.py`
- **Databento** support (paid, 16 years of CME/CBOT/NYMEX/COMEX plus ICE, equities, options):
  `db_cost()` estimates USD and record count BEFORE spending anything, and `db_download()`
  refuses to run without `confirm=True`. The API key is read from `DATABENTO_API_KEY` in the
  environment and is never accepted as an argument, printed, or written to disk.
- **Binance** free crypto: `binance_months()` / `binance_fetch()` pull the public monthly kline
  dumps — 1-minute bars from 2017-08, no API key, no rate limit. Verified reachable (HTTP 200,
  about 2 MB per month per symbol).
- **`CFD_TO_FUTURES` + `cfd_note()`**: a mapping from what CFD brokers call an instrument to the
  listed contract that actually has history. CFDs are OTC, so no vendor sells canonical CFD
  history — the honest method is to test on the underlying and add your broker's spread and
  financing as measured costs. 20 common CFD names mapped.

### Measured
- Databento historical pricing is $0.50/GB with $125 of signup credit, i.e. about 250 GB.
  NQ and GC both carry history from 2010-06-06 with OHLCV-1m available, so the datasets needed
  for an intraday study cost well under a dollar.

## 1.2.0 — 2026-07-26

### Added
- `references/07-montecarlo.md` — formulas 45-49.
- `qmontecarlo.py` — the null distributions a t-statistic assumes away: permutation and sign-flip
  tests (with the note that permuting order cannot test a mean), surrogate paths under four
  explicit no-edge processes (block/iid bootstrap, GBM, OU), White's Reality Check and Hansen's SPA
  for post-search comparison against a benchmark, and moment stress across six adversarial scenarios.
- `qdata.py` — the only networked module: public OHLCV history, your own CSV files, an idempotent
  forward collector, and `availability()` which MEASURES what a provider actually returns.
- 35 further assertions (255 total), including nulls validated against known cases: noise gives
  p > 0.10, a real edge p < 0.01, permutation ignores order on iid data but detects timing structure,
  Reality Check rejects a screen of eight noise strategies and finds a genuinely superior one.

### Documented (measured, not assumed)
- Free intraday history is the binding constraint on evaluating intraday strategies: ~20 years of
  daily bars and ~2 years hourly, but only ~60 calendar days of 5-minute (~50 sessions) and ~7 days
  of 1-minute. Older intraday windows return empty. A daily strategy is answerable from free data;
  an intraday one requires user-supplied history or forward collection. `SKILL.md` now checks data
  availability before agreeing to an analysis.

## 1.1.0 — 2026-07-26

### Added
- `qportfolio.py` — makes formulas 9-15 executable: covariance/correlation, portfolio volatility,
  marginal risk contributions, risk parity, Jacobi eigendecomposition, PC1 share, condition number,
  Marchenko-Pastur noise edge, shrinkage, capped long-only minimum variance.
- **`tail_dependence()`** — `λ = P(X ≤ q | Y ≤ q)`, the unbiased answer to "do these fail together".
- `qreport.py` — assembles an analysis into one self-contained HTML report: assumptions,
  measurability, full metric set with standard errors, paired difference, contraction check,
  overfitting diagnostics, sizing, five inlined SVG charts, verdict, and an explicit
  "what this evidence does not support" section. Objective register by construction.
- `examples/demo_portfolio.py` — third validation domain (four-sleeve book).
- 34 further assertions (220 total), including analytic portfolio identities (Euler risk
  decomposition, known eigenvalues, trace preservation, MP edge) and report-integrity checks.

### Fixed
- `solve_min_variance` violated the position cap: capping before renormalising let the
  renormalisation push weights back above it. Replaced with a water-filling projection.
- Correlation matrices could report a diagonal of 0 when an asset had zero variance in a subset;
  the diagonal is now 1 by convention.
- Documented a real statistical trap in `tail_conditional_corr`: Pearson correlation computed
  inside the worst decile of the portfolio mean is **biased downward** by the conditioning
  (Boyer-Gibson-Loretan), so a falling value is not evidence of improved diversification.
  `tail_dependence()` is now the reported measure; the biased one is kept only for comparison.

## 1.0.0 — 2026-07-26

First packaged release.

### Protocol
- Eight-step fixed order: restate, measurability gate, simulator validation, measure, uncertainty,
  contrarian pass, charts, verdict. Default action is "change nothing" absent decidable-positive evidence.
- Measurability gate promoted to step 2, before any backtest.

### Reference modules
- Formulas 1-30 across four modules (probability, linear algebra, stochastic calculus,
  pricing/risk), each entry with formula, practical use, **how it lies**, and a concrete action.
- Formulas 31-37 (`05-validation.md`): walk-forward, purged K-fold with embargo, Deflated Sharpe,
  PBO, HAC/Newey-West, MinTRL, multiple-testing/FDR.
- Formulas 38-44 (`06-sizing-execution.md`): Kelly and fractional Kelly, drawdown-constrained
  sizing, volatility targeting, market impact, capacity, turnover drag, latency and fill probability.
- `routing-table.md`: question-to-formula routing.
- `backtest-protocol.md`, `contrarian-protocol.md`, `charts.md`, `worked-example.md`.

### Code (stdlib only — no numpy, scipy, pandas or matplotlib)
- `qstats.py`: describe, paired difference, bootstrap (with block mode), jackknife, concentration,
  best-of-N headroom, required n, effective n, measurability gate, refinement factor, drawdown,
  empirical VaR/CVaR, Kelly ceiling.
- `qvalidate.py`: normal CDF/PPF, HAC standard errors, Deflated Sharpe, expected max Sharpe, PBO
  via CSCV, walk-forward and purged K-fold splits, MinTRL, Benjamini-Hochberg, Bonferroni.
- `qsize.py`: Kelly (continuous and discrete), growth rate, drawdown-constrained fraction,
  EWMA volatility, volatility targeting, square-root impact, capacity, turnover drag, cost ratio,
  effective rule distance under polling, adverse-fill adjustment.
- `qcharts.py`: equity, underwater, per-period contribution, parameter surface with noise band,
  return distribution with normal overlay — all as self-contained SVG. Plus ASCII sparkline.
- `tests/test_all.py`: 163 assertions against independent references (published t-tables, known
  normal quantiles, closed forms, analytic identities) and degenerate inputs.
- `examples/demo_strategy.py`: runnable end-to-end demo on a deterministic synthetic market.

### Fixed during development
- t-distribution tail was truncated, giving a 1% absolute p-value error at df=1 (Cauchy tails).
  Replaced with closed forms for df=1,2 and an analytic power-law tail remainder for df>=3;
  max error across 16 tabulated critical values now 3.3e-05.
- `dist_chart` ran off-canvas on a constant series (degenerate bin width).
- Zero-variance inputs produced ZeroDivisionError, 1e152 t-stats and a "strong evidence" verdict
  on NaN. All three now return explicit undefined states.
- `effective_n` could report more independent observations than the sample size under negative
  autocorrelation; now capped at n.
