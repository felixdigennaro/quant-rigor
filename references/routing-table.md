# Routing table — which formula answers which question

Start here. Find the question you actually have on the left; the formula and the one-line payoff
are on the right. Formulas 1-59 and 68 are explained in the `references/` file whose range covers
them (`01-probability.md` has 1-8, and so on); 60-67, 69 and 70 are explained outside `references/`,
and the last table below says where.

## The questions that come up most

| Question | Formula | What it gives you |
|---|---|---|
| "Is this backtest convincing?" | **1** Bayes | the posterior given the base rate of edges that die — usually sobering |
| "How many periods do I need?" | **8** LLN | `n ≈ (t·σ/effect)²`, converted into a calendar wait |
| "Can any backtest decide this Sharpe?" | **69** years to decide | `(t/Sharpe)²` years before it separates from zero, asked before writing the backtest; `literature_queue` sorts published strategies by decidability |
| "Is this difference real?" | **4** Variance + **7** Z | `SE = σ/√n`, then `t`; use the **paired** difference for two variants |
| "Can my data even measure this rule?" | **21** Wiener | fixed-level rules converge; path-dependent ones need `D ≥ 3-5×` the finest bar range |
| "I tried 60 variants and one won" | **8** + best-of-N | search headroom `σ√n·√(2 ln N)`; inside it = luck. Then **33** DSR, **34** PBO |
| "My variants are near-copies of each other" | **68** RAS | DSR assumes independent trials and over-penalises a correlated sweep; report `effective_trials` instead of the raw count, and use it alongside DSR, never in place of DSR |
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
| "Does this spread revert, fast enough to pay?" | **55-59** pairs | `hurst` and `adf_stat` (does it revert), `ou_half_life` (how fast, which sets trades per year and the cost hurdle), `engle_granger` and `kalman_hedge_ratio` (is it stable); `pair_report` runs `engle_granger`, with its ADF, half-life and Hurst, in one call, and `kalman_hedge_ratio` runs on its own |
| "What shape is my payoff?" | **26** Δ + **27** Γ + **25** BS | long or short gamma predicts the P&L shape before any backtest |
| "Do my costs matter?" | **16** Taylor + **43** turnover | cost as a **fraction of target**; edge needs ≥ 2-3× total cost |
| "What does a round trip actually cost here?" | **50** Roll + **51** Corwin-Schultz | the spread estimated from your own bars instead of assumed; report both, and read a disagreement as information about the data |
| "Is this short-horizon reversion real?" | **52** bounce check | if the implied spread is about one tick, the reversion **is** the cost of harvesting it |
| "How often can I trade before costs eat the edge?" | **53** cost cliff | `max_trade_frequency` returns the highest frequency at which the net Sharpe still reaches the target, and 0 when none does; if it is below what the idea needs, the idea is dead (with a rebate and a gross Sharpe below the target, the binding limit is a minimum frequency instead) |
| "Why is live P&L below the simulator?" | **54** slippage decomposition | splits the gap into entry slippage, exit slippage and commission, three different fixes; to check the model of the trade, compare its `realised_pnl` with the account statement |
| "Could I really have entered at that level?" | **70** fillability | `unfillable_share` counts entries at a level the price had already passed when the session opened; in the measured case they carried 106% and 96% of the profit |
| "Can I run this at size?" | **41** impact + **42** capacity | `impact ≈ Y·σ√(Q/V)`; capacity is where alpha = impact |
| "Is my t-stat honest?" | **35** HAC | overlapping/autocorrelated data inflates the naive t |
| "Do I have enough track record?" | **36** MinTRL | periods needed before an observed Sharpe is credible |
| "How do I validate without leaking?" | **31** walk-forward + **32** purged K-fold | purge overlapping labels, then embargo |
| "I tested 200 hypotheses" | **37** FDR | Benjamini-Hochberg instead of nominal p-values |
| "How often does no-edge do this?" | **45** permutation/sign-flip + **46** surrogates | a null built from your data, distribution-free |
| "I screened many strategies" | **47** Reality Check + **48** SPA | the only valid post-search tests against a benchmark |
| "Does the conclusion survive my own estimation error?" | **49** moment stress | vol +30%, costs x2, best period removed, mean shrunk |
| "It works better when X: should I size up on X?" | §5b `confluence_gate` | the five gates, starting from a monotone gradient across buckets instead of a binary filter, with a null that replays your search; `sizing_from_slope` turns the fitted slope into weights |
| "The edge stopped working: was it the regime?" | §5c `regime_report` | dilution lowers `t` but does not flip its sign, so a negative full sample needs losing years; strong windows are tested against a block-shuffled null |
| "Do I even have the data?" | `qdata.availability()` + **8** | ~20y daily, ~50 sessions of 5m — measure before promising |
| "Can this rule ever fire?" | **60** reachable threshold | the z-score of a point inside its own window `W` cannot exceed `√(W−1)` with population sd; a threshold of 2 on `W=4` (maximum 1.73) never fires |
| "Does the rule trade as often as I meant?" | **61** event rate | observed against intended rate, right after generating the signals and before the P&L; zero events is a rule that does not fire, not a negative result |
| "Is this drawdown in log terms?" | **62** log to simple | `−0.652` in log is −48%, not −65%; convert before publishing and put the scale in the label |
| "Is this spread estimate physically possible?" | **63** plausibility | two guards for disjoint cases (implied ticks, share of the bar range); no verdict when neither applies |
| "My verdict rests on someone else's number" | **63** borrowed parameter | `sign_flip_parameter` finds where the conclusion flips; if that falls inside the plausible range, deliver the threshold and go measure |
| "Is my breadth real?" | **64** effective breadth | refuses to call it breadth when `rho <= 0`; with a negative mean, diversifying makes the Sharpe worse (`diversification_sign_check`) |
| "Do my sentences match my numbers?" | **65-67** `qclaim` | `cmp_`, `rel`, `verdict` and `audit`: comparisons and verdicts generated from the values, and every untraceable number in a draft flagged |

## By stage of work

**Framing a new idea** → 2 (define the conditioning event *and* its frequency), 22 (write the SDE),
3 (EV per trade), 27 (which side of gamma you're on), 69 (can any backtest decide it?),
53 (can the trade frequency pay its costs?).

**Before backtesting** → 21 + measurability gate, 8 (required `n`), 43 (cost/target ratio),
50/51 (measure the spread), 52 (before any short-horizon reversion), 55-59 (half-life before any
pairs backtest), 60 (can the threshold fire?), 61 (once the signals exist, before any P&L).

**After a result** → 4, 7 (SE and t), 35 (HAC), paired difference, bootstrap, jackknife,
33/34 if you searched, 68 next to 33 when the variants are neighbours, 70 (could the entries at
a level really be filled?), 6 (correlation with what you already run).

**Before sizing** → 20 (median vs mean), 38 (Kelly ceiling), 39 (drawdown constraint),
5 (tail-conditional covariance), 30 (CVaR), 41/42 (impact and capacity), 64 (breadth, and the
sign of the mean first).

**Before showing anyone** → 24 (geometric, banded), 29 (alpha vs beta), 28 (Sharpe with its SE),
30 (what VaR cannot see), 62 (log to simple), 63 (no verdict on a borrowed number),
65-67 (sentences generated from the numbers), plus the disclaimer.

**When live P&L disagrees with the simulator** → 54 (where the money went), the first time it
happens, not the tenth.

## Executable support

| Module | Formulas it makes runnable |
|---|---|
| `qstats.py` | 3-8, plus the measurability gate (21) and best-of-N headroom; 60-62 and 69 |
| `qportfolio.py` | 9-15 and 64 |
| `qsize.py` | 38-44 (and 20's Kelly bound) |
| `qvalidate.py` | 31-37 and 68 |
| `qmontecarlo.py` | 45-49 |
| `qmicro.py` | 50-54, 63 and 70 |
| `qpairs.py` | 55-59 |
| `qclaim.py` | 65-67 |
| `qconfluence.py` | conditional edges and confluence sizing (SKILL.md §5b, no formula number) |
| `qregime.py` | "no edge works forever" (SKILL.md §5c, no formula number) |
| `qdata.py` | data loading (the only networked module) |
| `qcharts.py` / `qreport.py` | the five charts and the assembled report |

## The four that catch the most errors

1. **21 Wiener** — most "improved exit rules" are unmeasurable at the available resolution.
2. **Best-of-N headroom** (8) — most winners are the best of many, and inside the noise.
3. **35 HAC** — most t-stats on overlapping data are inflated.
4. **27 Gamma** — knowing you're short gamma predicts the tail before you observe it.

## Formulas explained outside `references/`

| Number | Function | Where it is explained |
|---|---|---|
| 60 | `qstats.reachable_threshold` | SKILL.md §0b, G1 |
| 61 | `qstats.event_rate_check` | SKILL.md §0b, G2 |
| 62 | `qstats.log_to_simple` | SKILL.md §0b, G3 |
| 63 | `qmicro.plausibility_flag`, `sign_flip_parameter`, `borrowed_parameter` | SKILL.md §0b, G4 and G5 |
| 64 | `qportfolio.effective_breadth`, `diversification_sign_check` | SKILL.md §0b, G6 and G7 |
| 65 | `qclaim.cmp_`, `rel`, `pct_change` | SKILL.md §0b, G10 (`pct_change` only in the CHANGELOG, v1.6) |
| 66 | `qclaim.audit` | SKILL.md §0b, G10 |
| 67 | `qclaim.verdict` | SKILL.md §0b, G10 |
| 69 | `qstats.years_to_decide`, `literature_queue` | SKILL.md §1b |
| 70 | `qmicro.unfillable_share` | SKILL.md §6, item 7b |

The numbers follow the CHANGELOG (60-64 in v1.5, 65-67 in v1.6, 68 in v2.1.1, 69 and 70 in
v2.1.2), and each module marks them with the same number in its section comments.
