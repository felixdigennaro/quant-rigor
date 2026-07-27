# Formulas 38-44 — Sizing, capacity, and execution reality

A real edge still loses money if it is sized or executed badly, and a validated backtest says nothing
about either. These seven decide **how much of a proven edge reaches the account** — and they are the
reason most strategies that pass `05-validation.md` still fail in production.

---

## 38. Kelly Criterion and fractional Kelly

```
continuous:  f* = (μ − r)/σ²                  μ, σ² per unit time, same units
discrete:    f* = p/a − q/b                   win fraction b w.p. p, lose fraction a w.p. q = 1 − p
growth:      g(f) = f·μ − f²·σ²/2             g(f*) = μ²/(2σ²) = SR²/2,   g(2f*) = 0
```

**What it does:** the fraction of capital that maximises the long-run **geometric** growth rate. Not
the fraction that maximises expected wealth (that is always "all of it") and not the fraction any
sane operator runs.

**Use in practice:**
- The continuous form is the one that matters: `f* = μ/σ²` in excess-return terms. With `μ = 15%` and
  `σ = 20%`, `f* = 0.15/0.04 = 3.75` — 375% of capital. The formula routinely recommends leverage no
  one can survive, which is the first clue about how to use it.
- `f* = μ/σ² = SR/σ`, so Kelly leverage is the Sharpe divided by the vol. Doubling assumed Sharpe
  doubles recommended size; halving assumed vol doubles it.
- The growth curve is **flat near the optimum and steep past it**. Fraction of maximum growth
  retained: `c = 1 → 100%`, `c = 0.5 → 75%`, `c = 0.25 → 43.75%`, where `f = c·f*`. Half Kelly buys
  75% of the growth for half the exposure. This asymmetry is the entire argument for fractional Kelly.

**How it lies:**
- **The drawdowns at full Kelly are not survivable.** For a GBM at `f = c·f*`, the probability of ever
  suffering a drawdown of `D` is `(1 − D)^(2/c − 1)`. At full Kelly (`c = 1`) that is `1 − D`: a **50%
  chance of ever halving your capital**, a 10% chance of ever losing 90%. At `c = 0.5` the chance of
  halving falls to 12.5%, at `c = 0.25` to 0.8%. The median worst drawdown is ~50% / ~21% / ~9% for
  `c = 1 / 0.5 / 0.25`.
- **`σ` and especially `μ` are estimated, and `f*` is hypersensitive to `μ`.** The relative error on
  `f*` equals the relative error on `μ`, which is `1/(SR·sqrt(T_years))`. With `SR = 1` and 3 years,
  that is ±58%. You can believe you are at `f*` while standing at `2f*`, where the growth rate is
  **exactly zero** and the drawdown distribution is catastrophic. Overbetting is not a smaller version
  of the right answer; past `2f*` it is a negative-growth strategy with a positive edge.
- `σ²` estimated in a calm window is too small, which inflates `f*` at precisely the wrong moment
  (see §40).
- The discrete form assumes a **known, stationary** payoff `(p, a, b)` and one bet at a time.
  Simultaneous correlated bets require the multivariate solution `f* = Σ⁻¹μ`, which inherits every
  pathology of the inverse covariance matrix (`02-linear-algebra.md` §10, §12) and is far less stable
  than the scalar version.
- Kelly optimises growth, and growth is not the objective if a drawdown triggers redemption, a margin
  call, or your own capitulation. Those are hard constraints, and Kelly does not know they exist.

**Do this:** compute `f*` from the **lower bound** of the `μ` confidence band, never the point
estimate; run at quarter-to-half Kelly; then take the **minimum** of that and the
drawdown-constrained size from §39 — and size on the minimum, not the average.

---

## 39. Drawdown-constrained sizing

```
P(max DD ≥ D  ever)  ≈  exp( −2μ·d / (f·σ²) )        with d = −ln(1 − D)
exact (GBM, log barrier):  P = (1 − D)^( 2μ/(f·σ²) − 1 )
inverted for f:      f = 2μ / [ σ²·( 1 + ln(p) / ln(1 − D) ) ]      p = tolerated probability
rule of thumb:       typical max DD  ≈  f·σ²/(2μ)  =  c/2   for f = c·f*
```

**What it does:** turns a drawdown tolerance into a position size. This is the constraint that should
usually bind, because a drawdown limit is a hard operational fact while a growth rate is a preference.

**Use in practice:**
- Same example as §38: `μ = 0.15`, `σ = 0.20`, `f* = 3.75`. Tolerate `D = 25%` with probability
  `p = 10%`: `ln(0.10)/ln(0.75) = 8.004`, so
  `f = 0.30 / [0.04·(1 + 8.004)] = 0.30/0.360 = 0.83`. **The drawdown constraint gives 0.83 versus
  Kelly's 3.75 — 22% of Kelly.** This is the normal outcome: the drawdown constraint is several times
  tighter than fractional Kelly, so Kelly is rarely the binding rule.
- The rule of thumb `D_typ ≈ c/2` is a useful sanity check and slightly pessimistic: it predicts 25%
  at half Kelly (exact median 21%) and 12.5% at quarter Kelly (exact 9.4%). Use it for mental
  arithmetic, the exact form for decisions.
- Halving `f` roughly halves the typical drawdown but cuts growth by only ~25% at the top of the
  curve — the same asymmetry as §38, read from the risk side.

**How it lies:**
- **The GBM assumption is the whole formula, and markets are not GBM.** Constant `μ`, constant `σ`,
  no jumps, no autocorrelation. Real drawdowns are worse than this predicts because vol clusters
  (losses arrive together) and prices gap (you do not get the intermediate prices). Treat the output
  as an optimistic bound and haircut it.
- **`μ` must be net of everything** — costs, fees, funding, financing, taxes. `f` scales roughly
  linearly in `μ`, so using gross drift overstates the safe size by the full cost ratio (§43).
- **The barrier probability is over an infinite horizon.** Over one year it is lower, which makes the
  infinite-horizon number conservative for annual planning — but the reverse error is fatal: quoting
  a one-year drawdown probability and then running the strategy for a decade.
- **A backtest's max drawdown is one draw from a distribution, not the distribution.** The realised
  worst is usually *milder* than the median of the true worst-case distribution for that path length.
  Block-bootstrap the return series and report the 50th and 90th percentile of max drawdown, not the
  single historical number.
- The formula is silent about **drawdown duration**, which is what actually ends strategies. A 20%
  drawdown lasting three months and one lasting three years are the same number here.
- It assumes you keep sizing at `f` through the drawdown. If you delever on the way down (§40) the
  drawdown is smaller and the recovery is slower — a different process than the one modelled.

**Do this:** state the tolerated `(D, p)` pair explicitly before sizing, compute `f` from net `μ`,
haircut it for jumps and vol clustering, and cross-check against a block-bootstrapped max-drawdown
distribution rather than the single backtested worst case.

---

## 40. Volatility targeting

```
f_t = σ_target / σ̂_t                        capped:  f_t = min( f_max, σ_target / max(σ̂_t, σ_floor) )
EWMA:  σ²_t = λ·σ²_{t−1} + (1 − λ)·r²_t      λ ≈ 0.94 daily
half-life = ln(0.5)/ln(λ) ≈ 11.2 days        centre of mass = λ/(1 − λ) ≈ 15.7 days
```

**What it does:** holds realised portfolio volatility roughly constant by scaling exposure inversely
to an estimate of current vol. It usually improves Sharpe and almost always makes the return
distribution *less* normal at the tails, not more.

**Use in practice:**
- Choose `λ` by the half-life you actually want, not by convention. `λ = 0.94` → 11 days;
  `λ = 0.97` → 23 days; `λ = 0.90` → 6.6 days. Shorter reacts faster and trades far more (§43).
- Always publish the realised-vol distribution *after* targeting. If you target 10% and realise a
  distribution centred at 13% with a fat right tail, the mechanism is not working and the rest of the
  sizing logic is built on a number that is not true.
- Report the turnover the scaler itself generates: roughly `E|Δf_t|` per period. A 10% daily move in
  `σ̂` trades 10% of the book for reasons unrelated to the signal.

**How it lies:**
- **Low measured vol is often low *information*, not low risk.** Close-to-close vol over a holiday
  week, a half-session, a market closure, or an illiquid overnight window is small because nothing
  traded. The scaler responds by **sizing up into the thinnest liquidity available** — the exact
  configuration in which impact (§41) and gap risk are worst.
- **Estimator latency guarantees you take the first leg of every crash at full size.** With an 11-day
  half-life, realised vol can triple before `f_t` halves. The mechanical sequence is: full size into
  the drop, then delever near the bottom, then re-lever after the recovery has happened. Vol targeting
  reliably crystallises the worst part of a fast selloff.
- **It systematically adds leverage in quiet regimes**, which is when carry, short-vol, and
  mean-reversion exposures build their largest positions — and quiet regimes are what precede jumps.
  The Sharpe improves and the left tail gets worse. Judge it on tail statistics, not Sharpe.
- **Vol is not risk for every exposure.** Anything with a discrete jump-to-zero risk (a venue, a
  counterparty, a peg) shows near-zero measured vol right up to the event, so the scaler grants it
  maximum size. Never let a vol scaler set the size of a jump-risk exposure; use a hard notional cap.
- Squared-return EWMA is a single-observation-per-period estimator and very noisy. Range- or
  realised-variance estimators from intraday data cut the noise substantially for the same half-life.
- The improvement is measured **gross**. Netting the scaler's own turnover cost removes a meaningful
  part of it and occasionally all of it.

**Do this:** cap `f_max`, floor `σ̂` at something like the 20th percentile of the trailing two years,
exclude non-trading and low-liquidity windows from the estimator, make the response **asymmetric**
(size down fast on rising vol, size up slowly on falling vol), and re-report the improvement net of
the turnover the scaler creates.

---

## 41. Market impact and non-linear costs

```
square-root law:   impact ≈ Y · σ · sqrt(Q/V)         Y ≈ 0.5-1 (calibrate; ~0.75 typical)
total cost:        C(Q) = Q·(spread/2 + fees) + Q·impact(Q)   →  the impact term grows as Q^1.5
```
`σ` = daily volatility in return units, `Q` = order size, `V` = daily volume in the same units.

**What it does:** describes how much the price moves against you as a function of how much of the
day's liquidity you consume. It is one of the most robust empirical regularities in market
microstructure and the single largest correction between a backtest and a real P&L.

**Use in practice:**
- With `σ = 1.5%` and `Y = 0.75`: participation `Q/V = 1%` costs `0.75·0.015·0.10 = 11 bps`;
  `Q/V = 10%` costs `0.75·0.015·0.316 = 36 bps`. **Ten times the size costs 3.2× per unit and 32× in
  total.** That non-linearity is the whole reason a strategy that works at small size stops working.
- Decompose the cost curve. At small `Q` the fixed half-spread dominates and cost per unit is roughly
  flat; past the crossover, impact dominates and cost per unit grows as `sqrt(Q)`. Find your crossover
  size — below it, size is nearly free; above it, every increment is paid for twice.
- Calibrate `Y` from your own fills (regress realised slippage on `sqrt(Q/V)`), do not inherit it.
  `Y` differs by venue, instrument, and time of day by a factor of two or more.

**How it lies:**
- **The law is calibrated on liquid large-cap metaorders executed over a full day.** It understates
  cost in thin books, in the first and last minutes of a session, in stressed markets, and for any
  order that is large relative to *visible depth* rather than to daily volume. In thin instruments
  impact can be closer to linear or worse.
- **`V` is endogenous to your signal.** If the signal fires when everyone is selling, the printed
  volume is high but the depth on your side is gone — so the formula reads "cheap" exactly when it is
  most expensive. Calibrate `V` conditional on signal days, not on all days.
- **It does not distinguish temporary from permanent impact.** Backtests that assume impact fully
  reverts book money that never existed. Assume the permanent component (typically a third to two
  thirds) is a pure loss.
- **A backtest filled at mid, at the close, or at the touch has already assumed zero impact.** The
  "alpha" it measured *contains* the impact you did not pay, so subtracting impact afterwards
  double-counts nothing — it corrects a number that was never gross of it.
- Splitting the order across time reduces impact but introduces timing risk `≈ σ·sqrt(execution
  horizon)`. There is a genuine trade-off; ignoring the timing-risk half makes slow execution look
  free.

**Do this:** put an explicit impact model in the backtest from the first run, not as a later
sensitivity; calibrate `Y` on your own fills and on signal-conditional volume; and report the P&L at
1×, 2×, and 5× the assumed cost so the reader can see where the edge dies.

---

## 42. Capacity

```
break-even size:   solve  α_per_unit = impact(Q)     →   Q_max = V · ( α / (Y·σ) )²
retention:         to keep fraction (1 − k) of the alpha,  Q = k² · Q_max
```

**What it does:** the size at which the expected alpha per unit is exactly consumed by the cost of
acquiring it. Beyond `Q_max` the strategy has a positive signal and a negative P&L.

**Use in practice:**
- `α = 20 bps` per round trip, `σ = 1.5%`, `Y = 0.75`, `V = $100M/day`:
  `α/(Yσ) = 0.0020/0.01125 = 0.178`, so `Q_max = 100M · 0.0316 = $3.16M` per trade.
- **The retention scaling is the number nobody computes.** `Q = k²·Q_max`, so keeping 50% of the alpha
  means trading at `0.25·Q_max = $790k`; keeping 75% means `0.0625·Q_max = $198k`; keeping 90% means
  `0.01·Q_max = $32k`. **The size at which the edge is mostly intact is one to two orders of magnitude
  below the break-even size.** Quoting `Q_max` as "capacity" describes the size at which you earn
  nothing.
- Set the operating limit at the retention level you require (typically 75-90%), then state the
  break-even size separately as the hard ceiling.

**How it lies:**
- **Capacity computed on average volume ignores that volume evaporates when the signal fires.**
  Recompute with the 10th-percentile daily volume of *signal days*; capacity typically drops by
  60-80%. During halts, limit-ups, or venue outages the tradeable capacity is zero at any price.
- **Capacity is not additive across strategies.** Two sleeves trading the same instrument in the same
  direction share one liquidity pool: joint capacity is roughly the single-strategy capacity, not the
  sum. Aggregate at the instrument level before summing.
- **It assumes `α` is independent of `Q` and of time.** Crowding decays `α` — and `Q_max` scales as
  `α²`, so a 30% decay in alpha halves capacity. Capacity estimated from a period when the trade was
  uncrowded is not the capacity you will get.
- **A backtest with mid-price fills has infinite capacity by construction**, which is why capacity is
  usually absent from backtest reports. Its absence is not evidence of abundance.
- Other constraints frequently bind before impact: borrow availability, margin and funding limits,
  position limits, settlement, and the sheer number of orders per second a venue accepts. Check which
  binds first rather than assuming it is impact.
- Capacity in *notional* terms silently changes when price changes: the same contract count is a
  different fraction of volume after a 40% price move.

**Do this:** report capacity as a band — mean volume and 10th-percentile signal-day volume — state
both the break-even size and the size retaining 75% of alpha, aggregate across sleeves at instrument
level, and re-measure from live fills once trading, never from the model alone.

---

## 43. Turnover and cost drag

```
net = gross − turnover × cost_per_turn
turnover (annual, one-way) = Σ_t Σ_i |Δw_{i,t}| / 2        state your convention
cost_per_turn = spread/2 + fees + impact(Q) + slippage + funding/borrow
```

**What it does:** the arithmetic that converts a signal's gross edge into what lands in the account.
It looks trivial and it is where most strategies actually die.

**Use in practice:**
- 250 round trips a year at 4 bps all-in is a **100 bps/yr** drag. Against a gross alpha of 300 bps,
  net is 200 bps and cost is **33% of gross**. Express cost as a fraction of gross, always — it is
  the only scale-free way to compare a daily strategy with a monthly one.
- The survival bar: gross must be **2-3× cost**, i.e. cost ≤ 33-50% of gross.

| gross / cost | interpretation |
|---|---|
| < 1 | negative net. Not a strategy at any size |
| 1-2 | positive only while costs behave; dies in the first stress regime |
| 2-3 | marginal. Requires live cost measurement before any sizing |
| 3-5 | survivable. Re-measure costs quarterly |
| > 5 | robust to cost error; check capacity (§42) instead — the edge is probably small in absolute terms |

- The 2-3× bar is not conservatism. It exists because (a) costs in stress are 2-4× benign-regime
  costs, exactly when you trade most, (b) gross alpha decays over time and costs do not, and (c) the
  gross estimate has a wide confidence band while the cost estimate is biased *low*.

**How it lies:**
- **`cost_per_turn` is almost always quoted as the exchange fee alone.** The half-spread, impact,
  funding, borrow, financing, and the cost of the rebalance trades are omitted. A "0.5 bps venue fee"
  becomes 4-8 bps all-in, which is an 8-16× error in the drag.
- **Backtests fill at mid or at the signal bar's close, so the half-spread is silently zero.** For a
  strategy whose per-trade edge is a few ticks, the half-spread *is* the edge.
- **Turnover computed on target weights, not executed weights**, misses the trades you actually did:
  partial fills, re-entries, forced closes, and the vol-scaler's own trading (§40).
- `cost_per_turn` is a **function of size** (§41), so a constant per-turn number is only valid at one
  size. Reporting one constant hides the fact that the net edge shrinks as you grow.
- **High-turnover strategies show the highest backtest Sharpes** precisely because per-trade edge is
  tiny and `n` is huge. That makes the entire result an artefact of the cost assumption: the same
  statistics with a realistic cost per turn are often negative. High Sharpe plus high turnover should
  raise suspicion, not confidence.
- Netting once, at the end, understates compounding drag: costs are paid from capital that would
  otherwise compound.

**Do this:** put full costs inside the backtest loop from the first run; report the gross/cost ratio
and the breakeven cost per turn at which the edge reaches zero; and reconcile modelled costs against
real fills and the account balance to the cent before believing either.

---

## 44. Latency and fill probability

```
E[P&L per signal] = P(fill)·E[edge | fill] − (1 − P(fill))·(opportunity cost)
adverse selection:  E[edge | fill]  <  E[edge]              always, for a resting limit order
latency slippage ≈ σ_per_second · sqrt(latency)  + queue effects
effective trigger of a polled path rule ≈ D + (range over the polling interval)
```

**What it does:** accounts for the fact that intending to trade is not trading. A backtest assigns
fills by rule; a market assigns them by queue and by who is in a hurry.

**Use in practice:**
- **A resting limit order fills when the queue ahead of you clears — which requires pressure through
  your price.** So conditional on being filled, the short-horizon move is against you. This is not an
  edge case; it is the structural reason market making is hard, and it means an unmodelled fill
  assumption imports a negative expectancy.
- **Backtest fill conventions are not a detail.** "Touch" (price traded *at* your level → assume
  fill) versus "trade-through" (price traded strictly *through* your level → fill). Touch is wildly
  optimistic for a limit at a level with a deep queue; trade-through is the conservative bound. Run
  both. If the result flips, the strategy is a bet on queue position, not on the signal, and the
  backtest cannot evaluate it.
- **Latency cost scales as `sqrt(latency)`.** With `σ_daily = 1.5%` over a 6.5-hour session,
  per-second `σ ≈ 0.98 bps`. A 200 ms round trip costs ~0.44 bps; 5 seconds costs ~2.2 bps. So 25×
  the latency is only 5× the cost — but on a 3 bps edge, 2.2 bps is most of it.
- **Polling widens every path-dependent rule.** A stop, trailing stop, or running-max rule evaluated
  on candle closes or every `Δt` seconds is *not* the rule you specified: its effective tolerance is
  the price range over `Δt`. A 3-tick stop polled on 1-minute candles whose median range is 4 ticks
  triggers somewhere in a 4-tick band — the measurement error on the rule's expectancy is as large as
  the expectancy. The distance `D` must be `≥ 3-5×` the median range of your finest *achievable*
  evaluation interval, or you are measuring your own blindness
  (`03-stochastic.md` §21). Note this is the **achievable polling** interval, not the data interval —
  tick data does not help if the bot checks once a second.
- **Re-emission compounds it.** Each cancel-and-repost cycle of duration `τ` lets price travel a
  further `~σ·sqrt(τ)`, so the realised trigger is `D + k·σ·sqrt(τ)` after `k` cycles.

**How it lies:**
- **Fill ratio is assumed to be 1 and never reported.** A strategy with a 40% fill rate has 40% of
  the trade count and therefore 40% of the statistical evidence — recompute `MinTRL`
  (`05-validation.md` §36) on realised fills, not intended signals.
- **The unfilled cases are not free.** They are selectively the profitable ones (the price ran away
  because you were right), so `E[edge | fill]` is biased down twice: adverse selection on the fills
  and positive selection on the misses.
- **Market orders remove fill uncertainty and replace it with impact and spread** (§41). "We'll just
  use market orders" converts an unmeasured fill assumption into a measurable cost — usually a large
  one, which is why the backtest was run with limits.
- Stops modelled as filling *at* the stop price are optimistic by definition; a stop is a market order
  triggered at the worst moment, so model it as touch-plus-slippage, gap-aware.
- Backtests that both enter and exit inside the same bar have decided the intrabar path. Any rule
  whose entry and exit can occur within one bar is untestable at that resolution regardless of how
  good the code is.

**Do this:** run every path-dependent backtest under both fill conventions and at the polling
resolution you can actually achieve; measure realised fill ratio and realised adverse selection live
on small size before sizing up; and if the edge exists only under touch fills or only below your
achievable polling interval, the verdict is "not decidable", not "promising".

---

Next: back to `05-validation.md` for the inference discipline, or `backtest-protocol.md` for the
order in which to run all of this.
