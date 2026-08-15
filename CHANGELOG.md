# Changelog

Versioning: MAJOR.MINOR.PATCH. MAJOR changes the protocol's steps or verdict semantics; MINOR adds
formulas, modules or functions; PATCH fixes correctness or presentation.

## 1.5.0 — 2026-08-15

### Added
- `examples/demo_csp_leveraged_etf.py` — a second circulating retail claim through the protocol:
  wait for a red day on SOXL/TQQQ, confirm IV > 50%, sell a 30-day cash-secured put, close at
  50% of credit, rotate the full amount immediately, "3-4 times a month until a $5,000 week is
  the standard". Reuses the engine, validated 3x simulator and IV surface from
  `demo_wheel_leap_pmcc.py`.
- `to_tradeable()` in `demo_wheel_leap_pmcc.py` — rebuilds a split-back-adjusted series at the
  price level it was actually quoted at. Yahoo reads TQQQ at $0.22 in 2010; returns are right but
  level is not, and level decides option tick size, the $0.01 bid-ask floor and contract notional.
  Without this the real-TQQQ backtest is destroyed by a spread artefact rather than by the strategy.
- `Config.entry_red_day`, `Config.entry_min_iv`, `Config.do_covered_call`, plus per-trade P&L and
  holding-period logging in the engine.

### Measured
- **"IV above 50%" is not a filter on a 3x ETF.** TQQQ 30-day ATM IV is >= 50% on **84% of days**
  (mean 67%, p5 45%). The threshold binds only in the calmest quarter of the sample — the quarter
  with the smallest premiums. IV >= 80%, which binds on 21% of days, does select: +3.46% per trade
  (t = 3.90) against +1.04% unfiltered.
- **The red-day filter is refuted as an edge**: −0.526% per trade against no filter (t = −0.60),
  and it lowers full-system CAGR from 13.60% to 12.73%. Selling after a down day buys higher IV and
  a higher probability the underlying keeps falling; on a 3x ETF the two cancel.
- **"$5,000/week"**: the claimed variant earns $372/week per $100k of collateral on average, so the
  figure needs about **$1.35M of collateral** — and the same configuration draws down 69%.
  Measured rotation rate is 1.5/month, not 3-4; mean holding period is 14.1 sessions.
- **Assignment is the mechanism the claim omits.** On real TQQQ the position is assigned 22 times in
  16.5 years and holds 3x shares for **31% of all sessions**, always after a fall: CAGR 1.21% against
  42.82% for buy-and-hold TQQQ. Robust to the post-assignment rule (1.44% without covered calls).
- **400 surrogate markets**: beats buy-and-hold QQQ on CAGR 24.8% of the time, on drawdown **0.8%**,
  on **both 0.5%**. Median CAGR 10.1% vs 17.1%, median MaxDD −60.9% vs −33.4%.

## 1.4.0 — 2026-08-15

### Added
- `examples/demo_wheel_leap_pmcc.py` — the protocol applied end-to-end to a widely circulated
  retail options claim: sell 30-day cash-secured puts 5% OTM, close at 50%, skim 30-50% of every
  winning premium into a "LEAP fund", buy deep-ITM 18-24 month calls on TQQQ/SOXL, then sell
  weekly 25-delta calls against them. Self-contained options engine (Black-Scholes, a stressable
  IV surface anchored on VXN with term structure and skew, explicit commissions and half-spreads,
  assignment, wheel completion, LEAP rolls, share splits).

### Measured
- **3x-ETF simulator validated before use** (protocol step 3): `3*r_QQQ - 1.25 bp/day` reproduces
  real TQQQ over 16.5 years with residual sd 20.7 bp/day against a 386 bp/day underlying — R² 99.71%,
  cumulative divergence 0.48%/yr. The 2.96 IV multiplier for the leveraged ETF is the measured
  realized-vol ratio, not an assumption.
- **The ablation is the finding.** Over 20 years the wheel alone returns 8.45% CAGR against 16.65%
  for buy-and-hold QQQ. The full six-step system's headline 35.81% comes almost entirely from
  recycling step-6 premium into an uncapped number of LEAP contracts (peak 1,721 held at once);
  capped at 5 contracts it returns 11.58%, below the index.
- **Step 6 has no measurable per-trade edge.** Across 894 weekly cycles the poor-man's-covered-call
  earns $8.23 ± $7.26 per contract (t = 1.13, p = 0.26), skew −3.61, win rate 85% — the signature
  of selling volatility, and the total sits inside the best-of-N search headroom.
- **400 surrogate markets** (block bootstrap of (return, VXN) pairs, whole engine re-run per path):
  the system beats buy-and-hold on CAGR 42.2% of the time, on drawdown 22.2%, and on **both — the
  claim actually being made — 6.8%**.
- **Net delta goes negative in fast rallies**: a 0.90-delta LEAP against a 0.25-delta weekly call is
  net short the underlying after a ~30% move in the 3x ETF, i.e. short exactly the move the LEAP exists to capture.
- **The 1-2% premium claim requires 30%+ IV**; mean VXN over the sample was 22.6%. Rebuilding the
  underlying at 1.3x/1.6x/2.0x volatility raises the premium as promised and leaves wheel CAGR flat
  while drawdown roughly doubles.

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
