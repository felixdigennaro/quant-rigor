# Worked example — a trailing stop that looks good and isn't

Every number below is produced by `examples/demo_strategy.py`, which generates a deterministic
synthetic market (no downloads, no dependencies). Run it and you get exactly these figures. The
point is not the strategy — it's watching the protocol catch four specific failures that a normal
backtest reports as success.

**Setup.** An intraday fade: a zone of ±40 points around the session open; if price extends 60
points beyond the zone edge, enter *against* the move; three rungs at 60/90/120 points; take-profit
back at the zone edge; fixed 60-point stop; flat at the end of the session. Costs measured, not
assumed: 2.0 fee + 2.5 slippage on market exits. 60 sessions, 20 coarse bars each, and the
generator retains the fine sub-path so resolution effects can be *measured* rather than argued about.

**The question.** Would a trailing stop improve this? Twelve variants were searched
(trail ∈ {20, 30, 40, 60} × activation ∈ {0, 20, 40}).

---

## Step 2 — Measurability gate (before any backtest)

| rule | `D / median bar range` | verdict |
|---|---|---|
| fixed stop 60 | R = 1.59 | not measurable — treat with care |
| trailing stop 40 | R = 1.06 | not measurable — artefacts |
| trailing stop 20 | **R = 0.53** | **fiction — the rule lives inside one bar** |

Measured path-refinement factor: **1.74** against a random-walk null of `√n = 2.24`. The path keeps
revealing more oscillation as resolution improves — it does not converge (`03-stochastic.md` §21).
So the trailing variants are, at best, weakly measurable with this data. **This is known before a
single P&L number is computed**, and it already tells you how much to trust what follows.

## Step 4 — Measure

Baseline: **+2288 points over 43 triggering sessions** (17 of 60 never triggered), outcomes
40 TP / 18 stop / 24 timeout. Cost ratio against the nearest 60-point target: **8%** — fine
(`06-sizing-execution.md` §43: an edge needs to clear 2-3× its costs).

## Step 5 — Uncertainty

| | baseline | best of 12 trailing variants |
|---|---|---|
| total | **+2288** | +1617 |
| mean/session | +53.2 ± 13.4 (SE) | +37.6 ± 10.8 |
| t | **+3.96** | +3.49 |
| max drawdown | −204.5 | −193.5 |
| win rate | 79% | 79% |
| skew | −0.38 | −0.01 |
| top-3 sessions' share | 15% | 15% |
| jackknife sign survives | 43/43 | 43/43 |

Paired difference: **−15.6 points/session ± 7.3, t = −2.13, p = 0.040**. The best trailing variant
is *worse*, and significantly so on a paired test.

Best-of-12 search headroom: **±1288 points**. Observed gap: **−671**. So even this negative result
sits **inside the noise band** — the honest reading is "no evidence the trailing helps", not
"proof it hurts by exactly 671".

## The four traps, caught

**1. Resolution.** The most attractive variants are the tightest trails — and those are precisely
the ones the data cannot resolve (R = 0.53). Widen the trail until it *is* resolvable and it
converges on the existing 60-point stop: the idea collapses into what the strategy already does.

**2. Contraction — algebra, not merit.** On the **9 sessions where the baseline lost**, the
baseline made **−586** and the variant **−235**. The variant "wins" there *by construction*: it
caps both tails, so it beats the original on any window where the original is negative. Selecting a
variant because it rescues bad periods is selecting on an identity, not on evidence.

**3. Best-of-N.** Twelve configurations were searched over 43 sessions with σ ≈ 88/session. The
headroom `σ·√n·√(2 ln N)` is ±1288 points — larger than any difference observed. And the
overfitting diagnostics agree: **PBO = 0.25** (substantial overfitting risk in the *selection
procedure itself*).

**4. Autocorrelation.** The baseline's naive t is +3.96; with HAC/Newey-West errors it is **+4.61**
(SE inflation 0.86×, 3 lags). Here the correction happens to *help* — the returns are mildly
negatively autocorrelated. The lesson is not the direction, it's that the naive t was the wrong
number and nobody would have known.

## Steps 5b, 6 — diagnostics and sizing

- **Deflated Sharpe** of the baseline: hurdle `SR0 = +0.138`, **DSR = 0.992** → established.
  Note the asymmetry this exposes: the *baseline* clears its hurdle; the *modification* does not.
- **MinTRL** for the baseline's Sharpe: **10 sessions**. It has 43 — comfortable.
- **Sizing:** Kelly ceiling `f* = 68.5`, quarter-Kelly `17.1`, but the drawdown constraint
  (≤15% with 10% probability) allows only **8.9** → the **drawdown constraint binds**, not Kelly.
  Sizing on Kelly alone would have been ~2× too aggressive.

## Step 8 — Verdict

> **NOT DECIDABLE** for the trailing stop: unmeasurable at this resolution, and its gap sits inside
> the best-of-N noise headroom. **DEFAULT ACTION: change nothing.**

The baseline, by contrast, *is* decidable-positive on this sample (t = +3.96, DSR = 0.992,
jackknife 43/43, only 15% concentrated in the top 3 sessions) — with the standing caveat that a
single synthetic regime is one observation, not 43 independent ones.

---

## Why this is the valuable output

A conventional report would have said: *"trailing stop tested, 12 variants, best one delivers +1617
with a smoother equity curve and a shallower drawdown — recommend adoption."* Every clause is true.
The conclusion is wrong.

Three habits produced the difference, and they generalise:

1. **Classify the rule before testing it.** Fixed-level exits (limit TP, fixed stop) are measurable
   on coarse data because they depend on the path's *extremes*, which converge. Path-dependent
   exits (trailing, running max, breakeven-move) depend on total variation, which diverges — they
   need `D ≥ 3-5×` the finest available bar range, and usually you do not have it.
2. **Check the algebra before the statistics.** Anything that caps both tails is a contraction and
   will beat the original wherever the original lost. That comparison carries zero information.
3. **Report the full metric set for the exact configuration you are recommending** — not another
   variant's number, and never a point estimate without its standard error and the search headroom.

A negative result that prevents a change to a working system is not a failed analysis. It is the
highest-value output this protocol produces.
