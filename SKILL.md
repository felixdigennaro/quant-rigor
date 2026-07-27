---
name: quant-rigor
description: Rigorous quantitative analysis of any trading or investment strategy — existing, proposed, or under evaluation. Applies 44 quant formulas (probability, linear algebra, stochastic calculus, pricing/risk, validation, sizing/execution) as working tools rather than decoration, each with its failure mode. Runs a measurability gate before any backtest, validates the simulator against reality, quantifies uncertainty (SE, t, bootstrap, jackknife, best-of-N headroom, Deflated Sharpe, PBO, HAC), runs an adversarial contrarian pass, and produces the charts that expose what tables hide. Use whenever a strategy, edge, backtest result, parameter change, sizing decision, or performance claim is on the table.
---

# Quant Rigor — analyse before you believe

> A number without its standard error is an opinion. A backtest without an out-of-sample is a story.
> A strategy that "works" on 15 days works on 15 days.

This skill turns the 30 core quant formulas into an operating procedure. It is not a formula sheet
— the formulas are in `references/`. This file is the **workflow**: the order of operations that
stops a plausible-but-wrong result from reaching a decision that moves money.

---

## 0. WHEN THIS FIRES

Fires on any of these, automatically, without being asked:
- A strategy is proposed, changed, or parameter-tuned ("what if we add a trailing stop / move the TP / raise size")
- A backtest result is reported or quoted (mine or someone else's)
- A performance claim is made ("this returns X%", "the edge is real", "win rate is Y")
- A sizing / leverage / risk-budget decision is on the table
- An existing live strategy is reviewed
- Someone asks whether an edge is real

Does NOT fire on: pure execution plumbing (order routing, API wiring) with no claim about returns,
or accounting of realised P&L already reconciled to an exchange balance (that's bookkeeping, not inference).

---

## 1. THE ORDER OF OPERATIONS (never reorder)

The sequence matters because each step can **kill the analysis**, and killing it early saves
the rest. Most wrong conclusions come from doing step 4 before steps 1-3.

```
1. RESTATE       what exactly is claimed, and what would falsify it
2. MEASURABILITY can the available data even measure this? (§2 — the step everyone skips)
3. VALIDATE      does the simulator reproduce reality where reality is known? (§3)
4. MEASURE       run it, with costs, on the full sample (§4)
5. UNCERTAINTY   band not point: SE, t, bootstrap, jackknife (§5)
6. CONTRARIAN    adversarial pass whose mandate is to destroy the finding (§6)
7. CHARTS        the four plots that expose what the table hides (§7)
8. VERDICT       decidable / not decidable, and what would make it decidable (§8)
```

Load `references/` files as needed — **start with `routing-table.md`** (question → formula):
- `routing-table.md` — which formula answers which question, by stage of work
- `01-probability.md` — formulas 1-8 (Bayes → LLN): inference, what a sample can support
- `02-linear-algebra.md` — 9-15 (matrices → Lagrange): portfolios, risk decomposition, optimisation
- `03-stochastic.md` — 16-22 (Taylor → SDE): price dynamics, why paths behave as they do
- `04-finance.md` — 23-30 (PV → VaR): discounting, options, Sharpe, CAPM, tail risk
- `05-validation.md` — 31-37: walk-forward, purged K-fold, Deflated Sharpe, PBO, HAC, MinTRL, FDR
- `06-sizing-execution.md` — 38-44: Kelly, drawdown-constrained sizing, vol targeting, impact,
  capacity, turnover drag, latency/fills
- `07-montecarlo.md` — 45-49: permutation & sign-flip, surrogate paths, White's Reality Check,
  Hansen's SPA, moment stress
- `backtest-protocol.md` — the full checklist, bias by bias
- `contrarian-protocol.md` — the attack playbook (run it, don't skim it)
- `charts.md` — the required plots and what each one is for
- `worked-example.md` — **a real case end to end**, with the four mistakes that shaped this
  checklist. Read this first if you want the skill in 5 minutes.
- `scripts/qstats.py` — describe / paired / bootstrap / jackknife / headroom / required_n /
  effective_n / measurability / refinement_factor / drawdown / VaR-CVaR / Kelly.
  `python3 qstats.py --a <base…> --b <variant…> --variants N` prints the whole decision set.
- `scripts/qvalidate.py` — HAC/Newey-West, Deflated Sharpe, expected-max Sharpe, PBO (CSCV),
  walk-forward & purged K-fold splits, MinTRL, Benjamini-Hochberg FDR.
- `scripts/qsize.py` — Kelly (continuous/discrete), drawdown-constrained fraction, EWMA vol &
  vol targeting, square-root impact, capacity, turnover drag, cost ratio, effective rule distance
  under polling, adverse-fill adjustment.
- `scripts/qportfolio.py` — covariance/correlation, **tail dependence** (the unbiased tail measure),
  portfolio vol, marginal risk contributions, risk parity, Jacobi eigendecomposition, PC1 share,
  condition number, Marchenko-Pastur noise edge, shrinkage, capped min-variance.
- `scripts/qreport.py` — assembles the whole analysis into ONE self-contained HTML report
  (assumptions · measurability · metrics with SE · difference · overfitting · sizing · charts ·
  verdict · what the evidence does not support). Objective register, no persuasive language.
- `scripts/qcharts.py` — stdlib-only, emits standalone **SVG** (no matplotlib):
  `equity_chart`, `drawdown_chart`, `contribution_chart`, `surface_chart(…, headroom=…)`,
  `dist_chart`, plus `sparkline()` for terminal output.

All four scripts ship with `tests/test_all.py` — 163 assertions against independent references
(published t-tables, known normal quantiles, closed forms) plus degenerate-input cases. Run it first.

---

## 2. MEASURABILITY GATE — can the data measure the claim at all?

**Run this BEFORE any backtest. It is the single highest-yield step and the most commonly skipped.**

An exit rule, a signal, or a stop is only measurable if the data's resolution is **much finer**
than the distance the rule cares about. Otherwise you are measuring your own blindness.

### The test
1. Compute the **typical bar range** for your data: `median(high - low)` over the relevant session.
2. Compute the rule's **characteristic distance** `D` (trail width, stop distance, target size).
3. Ratio `R = D / median_range`.
   - `R ≥ 3-5` → measurable
   - `1 ≤ R < 3` → **not measurable**; results are artefacts
   - `R < 1` → the rule lives entirely inside one bar; the backtest is fiction

### The convergence check (decisive)
Re-run the identical rule at **two or more data resolutions** (e.g. 5m and 1m):
- **Fixed-level rules** (limit TP, fixed stop) → the number **converges**: same result at both
  resolutions, because a fixed level is either touched in the bar or not. Resolution only
  reorders two levels inside one bar (1 bit of information, and only when both are touched).
- **Path-dependent rules** (trailing stop, running-max/min, breakeven-move, anything defined
  ON the path) → the number **does not converge**: it keeps changing as you refine, because
  measured path length grows like `√(1/Δt)` (Wiener process — see `03-stochastic.md` §21).
  Measure the refinement factor: `sum(fine bar ranges) / coarse bar range`. For a random walk
  this is `√(n)`. If you measure ≈`√n`, the path has no intrinsic scale and **no resolution
  you have is enough** — only rules with `D` far above the finest bar range are testable.

### Do NOT do this
**Aggregating fine bars into coarse bars is not a granularity test.** Coarsening *hides* adverse
excursions, so it monotonically improves everything path-dependent — and eventually starts hiding
the baseline's stop too. It proves nothing. The informative direction is **refinement**, never coarsening.

### Verdict language
If `R < 3` and refinement doesn't converge, the honest output is:
> "Not measurable with this data. To answer it I need <resolution>, collected forward
> with a pre-registered criterion. What IS measurable now: <fixed-level alternatives>."

Offer the non-path-dependent alternative that answers the same *business* question
(e.g. instead of a trailing stop: partial TP, time-based exit, daily loss limit — all fixed rules).

---

### Data first — check before promising an answer
Run `qdata.availability(symbol)` (or `python3 scripts/qdata.py --symbol X`) BEFORE agreeing to
evaluate a strategy. Measured, not assumed: ~20 years daily, ~2 years hourly, ~60 calendar days of
5-minute, ~7 days of 1-minute. Intraday retention is a hard limit — older intraday windows return
empty. Therefore: a daily strategy is answerable from free data; an **intraday** one is not, unless
the user supplies history (`qdata.from_csv`) or a forward collector has been running
(`qdata.append_csv`). Combine with `qstats.required_n()` and say the number of periods needed out
loud. Never substitute a short intraday sample for an answer.

## 3. VALIDATE THE SIMULATOR AGAINST REALITY

Never trust a simulator you haven't tied to a known outcome.

- If there are **live fills**: reproduce them. Entry levels must match **to the tick**; realised
  P&L must reconcile within costs. Report the comparison explicitly.
- If there are no live fills: reproduce a **published/known** result of the same rule, or
  hand-compute 2-3 periods and diff.
- Reconcile the **cash**: `start balance + Σ realised P&L + funding − fees = end balance`.
  If it doesn't tie to the cent, stop — something in the model is not what you think.
- Then, and only then, believe its output on the untested variants.

Also validate the **direction of conversion**: points → currency → % of capital. A 10× error here
is invisible in the ratios and fatal in the sizing.

---

## 4. MEASURE — with costs, on everything

- **Full sample, no cherry-picking of periods.** Report `n` (independent observations, not rows).
- **Costs from reality**, never assumed: fee per round-trip and slippage per exit type,
  measured from actual fills. Express them as a **fraction of the target**:
  `cost / target`. If a rule's target is 40pt and round-trip cost is 7pt, that's **18%** — a rule
  whose cost ratio is 3-5× the baseline's is a different economic proposition, not a tweak.
- **Adverse selection on costs**: exits that fire during momentum bursts (stops, trailing) slip
  worse than calm scheduled exits (time flat). Do not reuse one slippage number for both.
- **Unmodelled costs must be named**: funding/carry per hour held, borrow, taxes (never model
  taxes unless that is the task), overnight gaps. If a variant holds longer, its number is
  optimistic by an unquantified amount — say so.
- **Executability**: is the tested parameter the one that can actually run? A trail of 40pt
  polled every 15s behaves like ~50-55pt. A stop triggered on mark/oracle is not a stop on last.
  Test what you can execute, not what you can imagine.

---

## 5. UNCERTAINTY — the step that decides most questions

Use `scripts/qstats.py`. Formulas in `01-probability.md` (§3-8).

### Minimum reporting set
- `mean ± SE` where `SE = σ/√n`, and `t = mean/SE`
- For a **difference** between two variants, use the **paired** difference per period
  (`SE_diff = σ_diff/√n`) — it removes the common market factor and is far more powerful
- `n`, and the **effective** n (clustered/overlapping observations count for less)
- **Bootstrap** the metric (10k resamples) → 5/95 percentile band
- **Jackknife**: drop one period at a time; does the sign survive all `n` deletions?
- Contribution concentration: what fraction of the result comes from the top 3 periods?

### Thresholds, stated plainly
| t | reading |
|---|---|
| < 1.0 | noise; do not act |
| 1.0-1.6 | suggestive; a hypothesis, not a finding |
| 1.6-2.0 | weak evidence (p ≈ 0.05-0.11); act only with a structural argument |
| > 2.0 | evidence, still conditional on one regime |
| > 3.0 | strong, if the sample spans regimes |

### Formal overfitting diagnostics (use them, not just the headroom)
- **Deflated Sharpe** (`qvalidate.deflated_sharpe`) — corrects the Sharpe for the number of trials,
  non-normality and sample length. `DSR < 0.95` → not established.
- **PBO** (`qvalidate.pbo`) — probability that the *selection procedure* overfits. `> 0.2` is a warning,
  `> 0.5` means picking the in-sample best is worse than random.
- **HAC t** (`qvalidate.hac_t`) — the naive t is inflated on overlapping/autocorrelated data.
- **MinTRL** (`qvalidate.min_track_record_length`) — periods needed before a Sharpe is credible.
- **FDR** (`qvalidate.benjamini_hochberg`) — when testing many hypotheses, nominal p-values are void.

### Build the null, don't assume it (`qmontecarlo.py`)
- **sign_flip_test** — the distribution-free test for the MEAN. A permutation of the order cannot
  change a sum, so use sign flips here, not reordering.
- **permutation_test** — for claims about timing/sequence (`sharpe`, `max_dd`), block-shuffled when
  the data is autocorrelated.
- **surrogate_paths** — simulate a no-edge process calibrated to the observed volatility
  (block bootstrap / iid / GBM / OU) and locate the real result in that distribution. Run at least
  two nulls of different character; if the verdict depends on which, that dependence is the result.
- **reality_check / spa_test** — the correct tests after screening many candidates against a
  benchmark. A per-strategy p-value is not.
- **stress_moments** — does the conclusion survive a 30% volatility understatement, doubled costs,
  the best period removed, a shrunk mean, an extra tail event? Publish the failing rows too.

### The best-of-N trap (deflated Sharpe / multiple testing)
If you searched `N` variants and picked the best, the winner's edge is inflated by
`≈ σ_total · √(2 ln N)`. With `σ_total = σ_period·√n`:

> **Search headroom** = `σ_period · √n · √(2 ln N)`

Any result **below** the headroom is indistinguishable from luck. Report it:
> "Best of 60 variants, headroom ≈ ±2400pt; the winner beat baseline by 785pt → inside the noise."

Conversely — and this is the useful inverse — **if 0 of N variants beat the baseline while the
headroom was large, that is strong evidence the deficit is structural**, not sampling.

### Sample-size reality check (LLN, §8)
Needed periods for a given effect: `n ≈ (t_target · σ / effect)²`.
Run it *before* collecting, and say how long the wait is at the observed rate of new periods.
"15 sessions, ~1/weekday → 40-50 sessions is 5-7 weeks" is a decision-grade statement.

### The contraction trap (algebra masquerading as edge)
Any rule that **caps upside and downside** (trailing stop, partial TP, tighter stop, lower size)
is a **contraction**: `P&L_new ≈ α + β·P&L_base`, `0 < β < 1`.
A contraction **necessarily** beats the original on any window where the original is negative.
So "it wins on the recent bad stretch" carries **zero information** about merit. Test contractions
on the **full** sample, and never select one using a window chosen because the baseline lost there.

---

## 6. CONTRARIAN PASS — mandate: destroy the finding

Run `references/contrarian-protocol.md` in full. Spawn an independent adversarial agent for
anything that moves capital, changes strategy, or goes to a third party. Its brief is to **refute**,
not to review. Report whether the finding survived, and what it cost.

Non-negotiable attacks (details in the protocol):
1. **Measurability** — is the rule below the data's resolution? (§2)
2. **Algebra** — is the "win" mechanical (contraction, leverage identity, denominator change)?
3. **Selection** — was the winner chosen post-hoc, on a window used for another criterion?
4. **Concentration** — does the result live in 1-3 periods? Jackknife it.
5. **Regime** — one regime = one observation. Which regime is missing (trend, crash, chop, high-rate)?
6. **Executability** — poll rate, trigger basis (mark vs last), min size, margin, borrow.
7. **Costs** — measured or assumed? adverse selection? funding on longer holds?
8. **The missing number** — is the headline metric of the *recommended* variant actually reported,
   or was a different variant's number quoted? (This one catches real mistakes: always print the
   full metric set for the exact configuration being recommended.)

---

## 7. REQUIRED CHARTS

Tables hide shape; charts expose it. Produce these (see `charts.md`, `scripts/qcharts.py`),
and follow the `dataviz` skill for styling if a visual is being published:

1. **Equity curve** — baseline vs variant(s), same axes. Look for: is the gain a step (one period) or a slope?
2. **Underwater / drawdown curve** — depth AND duration. Two strategies with equal return and
   different time-under-water are not comparable.
3. **Per-period contribution bars** — sorted; instantly reveals concentration and which periods drive the delta.
4. **Parameter surface** — metric vs parameter. A **sharp isolated peak = overfit**; a broad plateau = robust.
   Never adopt a parameter at a lone spike; prefer the middle of a plateau.
5. **Return distribution / QQ vs normal** — where the tails actually are, and whether Sharpe/VaR
   (which assume normality) are lying (see `04-finance.md` §28, §30).

---

### Deliverable
Assemble the analysis with `qreport.Report` and hand over the HTML, not a terminal dump. The report
must contain the assumptions, the metric set for **every** variant (not just the recommended one),
the uncertainty, the charts, the verdict, and an explicit list of what the evidence does **not**
support. If a section would be empty, say so rather than omitting it silently.

## 8. VERDICT — say which of the three it is

1. **DECIDABLE, positive** — effect, band, t, why it survived the contrarian pass, what regime
   it has NOT seen. Give the sizing implication in currency and % of capital.
2. **DECIDABLE, negative** — the claim fails; state the mechanism that kills it (not just "not significant").
3. **NOT DECIDABLE** — say it plainly, then give the shortest path to decidability: what data,
   at what resolution, for how long, and the **pre-registered** criterion that will be applied to it.

Then the default action, explicitly: **change nothing** unless (1) is true. The burden of proof is
on the change, never on the status quo.

Close with the standard status line: **done · next · open risks**.

**Never present output as investment advice.** This protocol evaluates evidence; it does not
recommend positions to third parties, and backtested figures are hypothetical (see README disclaimer).

---

## 9. HOUSE RULES (adapt to your project)

- **Single source of numbers.** Every published figure must come from one designated script/query
  with its assumptions stated, or it doesn't get said. Check for drift by re-running it.
- **Bands, not points.** State the assumptions that move the band.
- **Base-rate discipline.** A previously retracted thesis is the base rate for the next one.
  Out-of-sample or nothing. Leverage that improves a backtest is a red flag until proven out-of-sample.
- **New strategies stay ring-fenced** from committed capital until a live, multi-regime track record exists.
- **Never quote an unvalidated edge** to investors, clients, or on a public page.
- **Never model taxes** unless that is explicitly the task.
- The analyst computes and verifies; a human signs and executes. No key handling, no order placement.
