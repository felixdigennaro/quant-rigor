# Routing table — which formula answers which question

Start here. Find the question you actually have on the left; the formula and the one-line payoff
are on the right. Every entry lives in a `references/` file with the same number.

## The questions that come up most

| Question | Formula | What it gives you |
|---|---|---|
| "Is this backtest convincing?" | **1** Bayes | the posterior given the base rate of edges that die — usually sobering |
| "How many periods do I need?" | **8** LLN | `n ≈ (t·σ/effect)²`, converted into a calendar wait |
| "Is this difference real?" | **4** Variance + **7** Z | `SE = σ/√n`, then `t`; use the **paired** difference for two variants |
| "Can my data even measure this rule?" | **21** Wiener | fixed-level rules converge; path-dependent ones need `D ≥ 3-5×` the finest bar range |
| "I tried 60 variants and one won" | **8** + best-of-N | search headroom `σ√n·√(2 ln N)`; inside it = luck. Then **33** DSR, **34** PBO |
| "Is it alpha or just beta?" | **29** CAPM | regress on the factor; if `R²` is high, rename the strategy after the factor |
| "Which of these two is better?" | **28** Sharpe | `SE(SR) ≈ √((1+SR²/2)/n)` — gaps below ~2·SE are noise |
| "How much should I size?" | **20** GBM + **38** Kelly + **39** DD | Kelly is a ceiling; the drawdown constraint usually binds first |
| "Why does leverage hurt?" | **19** Itô | the `−σ²/2` term: mean rises, median falls |
| "How bad can it get?" | **30** VaR/CVaR + **9** scenarios | empirical CVaR, plus scenarios for what VaR structurally cannot see |
| "Is this yield worth it?" | **23** PV | discount at `r_f` + venue hazard; the juicier one often loses |
| "What will €10k become?" | **24** FV + **20** | geometric (median) rate, with a band — never the arithmetic mean |
| "Am I actually diversified?" | **12** Σ + **11** eigen + **6** ρ | `λ₁/Σλ` and **tail dependence λ** (not correlation inside the tail — that measure is biased downward) |
| "What weight per sleeve?" | **15** Lagrange + **10** inverse | read the shadow prices; and prefer 1/N to an inverted sample covariance |
| "Which parameter do I pick?" | **13** gradient + **14** Hessian | middle of a plateau; sharp curvature = fragile fit |
| "What am I actually assuming?" | **22** SDE | OU (reversion) vs GBM (momentum) — a change that switches process needs its own evidence |
| "What shape is my payoff?" | **26** Δ + **27** Γ + **25** BS | long or short gamma predicts the P&L shape before any backtest |
| "Do my costs matter?" | **16** Taylor + **43** turnover | cost as a **fraction of target**; edge needs ≥ 2-3× total cost |
| "Can I run this at size?" | **41** impact + **42** capacity | `impact ≈ Y·σ√(Q/V)`; capacity is where alpha = impact |
| "Is my t-stat honest?" | **35** HAC | overlapping/autocorrelated data inflates the naive t |
| "Do I have enough track record?" | **36** MinTRL | periods needed before an observed Sharpe is credible |
| "How do I validate without leaking?" | **31** walk-forward + **32** purged K-fold | purge overlapping labels, then embargo |
| "I tested 200 hypotheses" | **37** FDR | Benjamini-Hochberg instead of nominal p-values |
| "How often does no-edge do this?" | **45** permutation/sign-flip + **46** surrogates | a null built from your data, distribution-free |
| "I screened many strategies" | **47** Reality Check + **48** SPA | the only valid post-search tests against a benchmark |
| "Does the conclusion survive my own estimation error?" | **49** moment stress | vol +30%, costs x2, best period removed, mean shrunk |
| "Do I even have the data?" | `qdata.availability()` + **8** | ~20y daily, ~50 sessions of 5m — measure before promising |

## By stage of work

**Framing a new idea** → 2 (define the conditioning event *and* its frequency), 22 (write the SDE),
3 (EV per trade), 27 (which side of gamma you're on).

**Before backtesting** → 21 + measurability gate, 8 (required `n`), 43 (cost/target ratio).

**After a result** → 4, 7 (SE and t), 35 (HAC), paired difference, bootstrap, jackknife,
33/34 if you searched, 6 (correlation with what you already run).

**Before sizing** → 20 (median vs mean), 38 (Kelly ceiling), 39 (drawdown constraint),
5 (tail-conditional covariance), 30 (CVaR), 41/42 (impact and capacity).

**Before showing anyone** → 24 (geometric, banded), 29 (alpha vs beta), 28 (Sharpe with its SE),
30 (what VaR cannot see), plus the disclaimer.

## Executable support

| Module | Formulas it makes runnable |
|---|---|
| `qstats.py` | 3-8, plus the measurability gate (21) and best-of-N headroom |
| `qportfolio.py` | 9-15 |
| `qsize.py` | 38-44 (and 20's Kelly bound) |
| `qvalidate.py` | 31-37 |
| `qmontecarlo.py` | 45-49 |
| `qdata.py` | data loading (the only networked module) |
| `qcharts.py` / `qreport.py` | the five charts and the assembled report |

## The four that catch the most errors

1. **21 Wiener** — most "improved exit rules" are unmeasurable at the available resolution.
2. **Best-of-N headroom** (8) — most winners are the best of many, and inside the noise.
3. **35 HAC** — most t-stats on overlapping data are inflated.
4. **27 Gamma** — knowing you're short gamma predicts the tail before you observe it.
