# Backtest Protocol — the checklist, bias by bias

Run top to bottom. Any unchecked box must appear in the output as a stated limitation, not be
silently omitted. A backtest with named limitations is useful; one with hidden limitations is worse
than none.

---

## A. Before writing any code

- [ ] **State the claim** in one falsifiable sentence, and what result would kill it.
- [ ] **Classify every rule** as fixed-level or path-dependent (`03-stochastic.md` §21).
- [ ] **Measurability gate** (`SKILL.md` §2): `D / median_bar_range ≥ 3-5`? If not → stop, declare
      not measurable, propose fixed-level alternatives.
- [ ] **Write the implicit SDE / thesis** (§22). Check the proposed change doesn't silently switch
      it (e.g. reversion → overshoot).
- [ ] **Pre-register**: which variants, which metric decides, what threshold. Written before running.
      Without this, everything later is post-hoc.
- [ ] **Compute required `n`** (`01-probability.md` §8). If required `n` >> available `n`, say so now.

## B. Data hygiene

- [ ] Provenance and timezone of every series (ET vs UTC bugs are silent and fatal; check DST).
- [ ] **Session filter**: are non-trading days (weekends, holidays) excluded? If the instrument
      trades 24/7 but the strategy is tied to a cash open, weekend "sessions" are fictitious.
      *Verify* they're excluded rather than assuming.
- [ ] Survivorship: are dead venues/instruments in the sample? Their absence inflates results.
- [ ] Gaps/outages/bad prints. A single bad tick can create or destroy an entire edge.
- [ ] Resolution of each series, and whether it changes mid-sample.
- [ ] Same data for signal and execution, or explicitly modelled lag.

## C. Look-ahead & leakage (the silent killers)

- [ ] Every input available **at decision time**? Timestamp each field.
- [ ] Bar close used for a decision executed at that same close → 1-bar leak. Execute next bar.
- [ ] Parameters fitted on the full sample then tested on it → in-sample. Split or walk-forward.
- [ ] Session/zone definitions that use data from later in the session.
- [ ] Restatements/revisions (rates, indices) — use point-in-time values.
- [ ] Filters chosen after seeing the result ("exclude the crash day") → leakage.

## D. Execution realism

- [ ] **Fills**: limit orders assumed filled on touch? Real fills need the level *traded through*
      or queue position. Touch-and-reverse tails often don't fill — this **inflates** any
      target-based strategy.
- [ ] **Costs measured from real fills**, split by exit type; maker vs taker.
- [ ] **Adverse selection**: momentum exits (stops, trailing) slip more than scheduled exits.
- [ ] **Trigger basis**: mark/oracle vs last price. On a builder/HIP-3 venue these differ (a
      persistent mark-vs-oracle offset was measured here at ~+9-10pt) — a stop on mark is not a
      stop on last.
- [ ] **Poll/latency**: a rule re-emitted every `k` seconds behaves like a wider rule. Test the
      executable parameter, not the ideal one.
- [ ] Min size / lot rounding / min notional; margin sufficiency at the moment of arming.
- [ ] **Funding/carry per unit time held** — if a variant holds longer, its number is optimistic
      until this is modelled. Name it if unmodelled.
- [ ] Liquidity/impact at the intended size; the backtest is silent about size.
- [ ] Never model taxes unless that is explicitly the task.

## E. Statistics

- [ ] `n`, effective `n`, mean ± SE, `t` (`01-probability.md`).
- [ ] **Paired** differences for variant comparisons.
- [ ] Bootstrap band (block bootstrap if autocorrelated) and jackknife sign-stability.
- [ ] Contribution concentration: top-3 periods' share of the result.
- [ ] **Multiple-testing correction**: report `N` searched and the `σ√n·√(2 ln N)` headroom.
- [ ] Out-of-sample or walk-forward. If impossible, say the result is in-sample — always.
- [ ] Regime labelling: which regimes are in the sample, which are missing.

## F. The traps that produce real errors

- [ ] **Coarsening ≠ granularity test.** Aggregating bars hides adverse excursions and improves
      everything monotonically. Only *refinement* is informative.
- [ ] **Contraction algebra.** Any upside-and-downside-capping variant beats the baseline on any
      window where the baseline is negative — mechanically. Never select on such a window.
- [ ] **The missing headline number.** Always print the full metric set for the *exact* variant
      being recommended. Quoting another variant's cost is how a −47% gets presented as −11%.
- [ ] **Zero-in-`n` bounds, not proves.** 0 bad events in 38 trades bounds the rate at ~7.6%
      (95%), i.e. up to ~2.9 expected per 38 — not "it never happens".
- [ ] **Sample-period luck.** Check whether the profit lives in a handful of exceptional periods,
      and report the result with and without them.
- [ ] **Wrong clock.** Verify the current date/time from the environment before reasoning about
      elapsed time or interest accrual.
- [ ] **Wrong endpoint/scope.** Query the exact venue/dex/subaccount the strategy uses; a default
      endpoint can return an empty or unrelated state and look like a catastrophe.
- [ ] **Unit conversion**: points → currency → % of capital, verified against a known trade.

## G. Sizing (separate from edge — never conflate)

- [ ] Kelly bound `f* ≈ μ/σ²` (or the discrete-payoff version) as a **ceiling**, then take a
      fraction of it. Above `f*`, more size lowers median wealth (`03-stochastic.md` §20).
- [ ] Worst realised period × current size = currency loss. Compare with the risk budget explicitly.
- [ ] Tail clustering: assume simultaneous stops across correlated legs (`01-probability.md` §5).
- [ ] Feedback loop check: does size depend on equity/free margin? Then the chain
      `size → margin → free margin → next size` exists — decide if you want it (`03-stochastic.md` §18).
- [ ] Stress the size: what if the worst period is 1.5× the historical worst?

## H. Output

- [ ] Full metric set for baseline and each variant: total, per-period mean ± SE, `t`, win%,
      max DD, time-under-water, worst period, skew, cost ratio, `n`.
- [ ] The five required charts (`charts.md`).
- [ ] Explicit verdict: decidable-positive / decidable-negative / **not decidable** + path to decide.
- [ ] Default action stated: **change nothing** unless decidable-positive.
- [ ] Limitations section listing every unchecked box above.
