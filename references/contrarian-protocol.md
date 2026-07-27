# Contrarian Protocol — mandate: destroy the finding

Not a review. The brief is **refutation**. A finding that survives a genuine attempt to kill it is
worth something; a finding that was never attacked is a story.

**Spawn an independent adversarial agent** (fresh context, told to demolish, given the data and the
method) whenever the conclusion: moves capital, changes a live strategy, goes to a third party, or
is irreversible. Give it the raw numbers and let it recompute — an adversary that must trust your
summary can't attack your summary.

---

## Attack 1 — Measurability (kill it before believing it)
- Is the rule's characteristic distance below the data's resolution? (`SKILL.md` §2)
- Does the metric **converge** under refinement, or keep moving? Measure the refinement factor;
  compare to `√n` (random-walk null, `03-stochastic.md` §21).
- Was a *coarsening* test mistaken for a granularity test? (It proves nothing.)
- If the rule must be widened until measurable, does it collapse into the existing rule?
  **If yes, the idea is vacuous** — a legitimate and decisive kill.

## Attack 2 — Algebra, not edge
Ask: is the improvement **mechanically guaranteed** rather than empirical?
- **Contraction**: does the variant cap both tails? Then it beats the baseline on every window
  where the baseline is negative, by construction. Any claim resting on such a window is void.
- **Denominator games**: did the metric improve because risk was removed from the denominator
  (lower size, lower vol) rather than return added? Compare at **equal risk**.
- **Leverage identity**: scaling size scales mean and σ together — Sharpe unchanged. A "better"
  result from size alone is not a better strategy.
- **Truncation**: does the variant simply take fewer trades? Compare per-trade, not per-period.
- **Survivor arithmetic**: does the improvement come from *not trading* the bad days, identified
  with hindsight?

## Attack 3 — Selection & post-hoc
- Was the winner picked from `N` variants? Compute the search headroom `σ_period·√n·√(2 ln N)`
  and check whether the edge clears it.
- Was the deciding window chosen **after** seeing that the baseline did badly there?
- Was the same data used to form the hypothesis and to test it?
- Are the reported metrics the ones pre-registered, or the ones that looked best?
- Is the chosen parameter an isolated spike on the surface (overfit) or mid-plateau (robust)?

## Attack 4 — Concentration & fragility
- Jackknife: does the sign survive dropping each period? Report the worst deletion.
- What share of the effect comes from the top 3 periods?
- Does it survive excluding the exceptional periods (the ones that carry the baseline)?
- Bootstrap: what fraction of resamples show the opposite sign?
- Is `n` really `n`? (overlapping windows, same-day correlated legs → effective `n` lower)

## Attack 5 — Regime & missing states
- How many regimes are in the sample? (One contiguous month = **one** observation, not `n` days.)
- Which states are absent: sustained trend, crash, gap/halt, high-rate, low-liquidity, holiday?
- What is the *mechanism* by which the finding fails in each absent state? Name it concretely
  (e.g. "on a gap the trail never arms → full stop, three times").
- Does the strategy's own thesis (its SDE, §22) still hold in the missing regime?

## Attack 6 — Executability
- Poll rate / latency: does the tested parameter equal the executable parameter?
- Trigger basis: mark/oracle vs last. Slippage measured for **this** exit type?
- Min size, rounding, margin at arm time, order-type support (does the venue even offer it natively?).
- Can the rule be raced/gamed by other participants once it's visible in the book?
- What happens on a partial fill, a rejected order, a disconnect mid-position?

## Attack 7 — Costs & carry
- Measured from real fills or assumed? Split by exit type?
- Cost as a **fraction of target** — does the variant change the economics
  (e.g. 4-7% of target → 18-36%)?
- Adverse selection on momentum exits?
- Funding/borrow for longer holds — modelled? If not, the longer-holding variant is optimistic by
  an unquantified amount.
- Does an earlier claim ("costs are negligible") still hold at the new target size? Justifications
  don't transfer across parameter changes.

## Attack 8 — Presentation integrity (catches real mistakes)
- Is the headline metric of the **recommended** variant actually reported, or was a different
  variant's number quoted?
- Are bands reported, or points?
- Is anything stated as verified that was only assumed? (Every unverified claim must carry
  `⚠ non verificato`.)
- Do the numbers reconcile to an external source (exchange balance, on-chain state) to the cent?
- Is the default action stated (change nothing) or is adoption assumed?

---

## Verdict format (required)

Per attack: **holds / falls / falls as formulated**, with the quantified reason.
Then one of:

- **"The finding survives"** + what it is *precisely* limited to.
- **"Correct it as follows"** + the corrected statement.
- **"NOT DECIDABLE with this data"** + shortest path to decidability (data, resolution, duration,
  pre-registered criterion).

And separately: **what remains decidable right now** — usually a *negative* result (a claim
definitively killed) is the most valuable output, because it stops a change to a working system.
