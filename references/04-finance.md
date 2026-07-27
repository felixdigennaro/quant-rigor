# Formulas 23-30 — Discounting, Options, Risk-Adjusted Return, Tail Risk

The applied layer. Also the layer where numbers get quoted to other people, so the failure modes
here are the ones that damage credibility rather than just an analysis.

---

## 23. Present Value

`PV = CF/(1+r)^t`  (continuous: `PV = CF·e^{−rt}`; multi-period: `Σ CFₜ/(1+r)^t`)

**What it does:** converts a future cash flow into today's terms. The single most under-used tool
in crypto-yield thinking.

**Use in practice:**
- **Every yield claim should be a PV.** "8% APY" from a protocol that might not exist in 12 months
  is not comparable to 4% from T-bills. Discount at a rate that includes the **default/venue
  hazard**: `r = r_f + spread`. If a venue has a 5%/yr hazard of loss, its 8% is worth ~3%
  risk-adjusted — often *less* than the safe alternative. This is how to compare a conservative venue with a higher-yielding riskier one, and it usually kills the "juicier" one.
- Airdrop valuation: `PV = P(receive)·E[value|receive]·e^{−rt}` with `t` = expected months to
  distribution. Makes an unpaid-airdrop pipeline valuable but discounted, not free money.
- Fee structures: the PV of "10% of gains above HWM" is not 10% of expected gains — it's an option
  (see §25) and worth more than the naive number when vol is high.

**How it lies:** the discount rate does all the work and is usually unstated. A hazard rate of
2% vs 10% changes conclusions completely. Also: discounting a *nominal* flow at a *real* rate
(or vice versa) is a silent 2-3%/yr error.

**Do this:** state `r` and its components (`r_f` + venue hazard + illiquidity) every time.
Show the PV at two hazard rates as a band.

---

## 24. Future Value

`FV = PV(1+r)^t`  (continuous compounding: `FV = PV·e^{rt}`)

**What it does:** compounding forward.

**Use in practice:** projections like "€10,000 over N years". Non-negotiables:
- Use the **geometric** (median) rate, not the arithmetic mean — i.e. `μ − σ²/2` (§19/§20).
  Projecting with the arithmetic mean overstates the typical outcome, and the error grows with
  horizon and vol.
- Net of **all** costs and fees, including the performance fee's drag on compounding.
- Show a **band** (e.g. 25th/50th/75th percentile paths from a block bootstrap), never a single
  curve. A single smooth curve is a marketing artefact.
- Label the rate's provenance and whether it's live or backtested.

**How it lies:** compounding a rate estimated over a short favourable window is the most common
way to produce an indefensible number. If the rate came from a 2.8-year window in which gold ran
+118%, the projection inherits that luck — say so **next to** the chart, not in a footnote.

**Do this:** in any FV projection state (1) the rate and its source, (2) that it's geometric,
(3) the band, (4) which regime it has never seen. The rate must come from your single source of numbers, or a live read.

---

## 25. Black-Scholes

`C = S₀N(d₁) − Ke^{−rT}N(d₂)`,
`d₁ = [ln(S₀/K) + (r + σ²/2)T]/(σ√T)`, `d₂ = d₁ − σ√T`

**What it does:** the no-arbitrage price of a European call under GBM + continuous hedging.
Its real value is as a **common language for optionality**, not as a price oracle.

**Use in practice — find the options hiding in your strategy:**
- A **stop-loss is a short put** on your own position (you cap the downside, you pay in whipsaw
  premium). A **take-profit is a short call** (you cap the upside, you collect certainty).
  Any change to stops/targets is therefore a change in option position — this reframing tells you
  immediately that "add a trailing stop" = "buy more downside protection, sell more upside", and
  the question becomes *is the premium fair*, not *does the backtest look better*.
- A **performance fee above a HWM is a call option** the manager holds on the fund's gains;
  its BS value is why fee structures reward vol. A performance-only fee aligns interest but is not free — compute its option value once, honestly.
- **Convexity budget:** if a sleeve's payoff is convex (small cost, large upside — the RWA overlay
  funded only by coupon), value it as a long option position. That's the correct framing for
  "upside at no capital risk": it's a call financed by yield.

**How it lies:** constant vol, no jumps, continuous costless hedging, log-normal tails. Real
option prices show a **volatility smile** precisely because the market knows BS's tails are too
thin. Implied vol is not a forecast; it's a price. And BS says nothing about liquidity.

**Do this:** use BS to *classify and price optionality qualitatively* (long/short convexity, and
roughly how much). For actual tail P&L, revalue at stressed prices (§16).

---

## 26. Delta

`Δ = ∂V/∂S`

**What it does:** first-order price sensitivity; the hedge ratio. For a call, `Δ = N(d₁) ∈ (0,1)`.

**Use in practice:**
- **If delta-neutrality is a core claim, it must be *measured*, not asserted**: `Δ_net = Σ (position notional × sensitivity to the common factor)`. Report residual
  delta in currency and as % of NAV, and re-measure after every rebalance.
- Basis/funding trades are delta-neutral in *price* but exposed to basis, funding, and
  **liquidation asymmetry** — hedged does not mean flat.
- Sizing a hedge: `hedge = −Δ_portfolio/Δ_hedge_instrument`.

**How it lies:**
- **Delta is local and moves** (that's gamma, §27). "Delta-neutral" is true at an instant.
- Two legs on different venues are neutral only if both survive; delta ignores counterparty.
- Delta computed against the *wrong* factor (e.g. hedging an equity-index perp against BTC) reads
  as neutral while carrying full exposure.

**Do this:** state which factor you are neutral to, measure residual delta in $, and list what
delta cannot see (venue, funding, oracle, liquidation).

---

## 27. Gamma

`Γ = ∂²V/∂S² = ∂Δ/∂S`

**What it does:** how fast delta changes — the curvature of the payoff.

**Use in practice:**
- **Long gamma** = you profit from movement, pay in decay/whipsaw (a straddle; a wide-stop trend
  follower). **Short gamma** = you collect steadily and lose in jumps (selling premium; a
  tight-stop mean-reverter; a fade strategy). **Know which you are.**
  A fade of an extension with a fixed stop is **short gamma**: many small wins, occasional
  large loss when the move continues. Every "improvement" that increases win rate while capping
  upside makes you *more* short gamma — i.e. sells more convexity. That is the correct lens for
  judging trailing stops and tight targets, and it predicts the observed shape (win rate up,
  tail unchanged) without any backtest.
- Hedging cost scales with `Γ·σ²` — a high-gamma book in high vol bleeds through rebalancing.

**How it lies:** gamma is huge near the strike/stop and near expiry — exactly where the linear
approximation fails and where slippage is worst. Backtests that assume clean fills at the stop
systematically understate the cost of being short gamma.

**Do this:** label every strategy long or short gamma, and check that the P&L distribution matches
(short gamma must show negative skew — if a backtest shows short-gamma structure with positive
skew, the tail is missing from the sample, not absent from the world).

---

## 28. Sharpe Ratio

`SR = (R_p − R_f)/σ_p`

**What it does:** excess return per unit of volatility.

**Use in practice:** comparing strategies **of similar shape** on the same horizon. Annualise
correctly: `SR_annual = SR_period·√(periods/yr)` — and only if returns are iid.

**How it lies — more ways than any other metric here:**
- **Assumes symmetric risk.** A short-gamma strategy (frequent small wins, rare large loss) posts a
  flattering Sharpe until the tail arrives. Selling deep OTM options has a superb Sharpe until it
  doesn't. Always report **skew and max drawdown** beside Sharpe.
- **Standard error is large:** `SE(SR) ≈ √((1+SR²/2)/n)`. With `n = 15` periods, `SE ≈ 0.27` —
  so `SR = 0.6` vs `SR = 0.55` is **indistinguishable**. Never rank variants on Sharpe differences
  below ~2·SE.
- **Annualisation inflates** with autocorrelated or overlapping returns (√n is wrong); it also
  makes higher-frequency strategies look better for purely mechanical reasons.
- **Trivially gameable** by smoothing marks (illiquid or model-priced assets), by truncating
  reporting periods, or by cherry-picking the window.
- `R_f` is not optional: in a 4-5% rate world, subtracting it changes rankings materially.

**Do this:** report `SR ± SE(SR)`, with skew, max drawdown, and time-under-water. Prefer
**Sortino** (downside deviation) or **return/max-DD** for asymmetric strategies. Never make a
decision on a Sharpe gap smaller than its standard error.

---

## 29. CAPM

`E(Rᵢ) = R_f + βᵢ(E(R_m) − R_f)`,  `βᵢ = Cov(Rᵢ,R_m)/Var(R_m)`

**What it does:** splits expected return into compensation for systematic risk (`β`) and everything
else (`α`).

**Use in practice — the honesty test for any "alpha" claim:**
- Regress the strategy's returns on the relevant market factor(s). If `β` explains the return,
  you have **beta with extra steps**, not edge. This is the single most useful check on a new
  strategy, and it's cheap.
- For crypto-yield: the factors are usually market beta, **funding/carry**, and **liquidity**.
  A basis trade's return is largely a carry factor — legitimate, but it should be *called* carry,
  priced as carry, and expected to compress.
- Residual `α` must clear its own standard error before being called alpha.

**How it lies:**
- **Empirically weak**: β doesn't explain the cross-section of returns; low-vol and momentum
  anomalies contradict it. Treat CAPM as an accounting decomposition, not a prediction.
- **β is unstable** and regime-dependent; estimated on short samples it's mostly noise.
- **"The market" may be undefined** for the strategy (what is the market portfolio for a Nasdaq
  perp fade on a builder DEX?). Choose the factor deliberately, state it, and show `R²`.
- Non-linear exposures (options-like payoffs) are invisible to a linear β — a short-gamma sleeve
  can show `β ≈ 0` and still be long the market's tail.

**Do this:** for every claimed edge, run the factor regression, report `α ± SE`, `β`, `R²`, and
name the factor. If `R²` is high, rename the strategy after its factor.

---

## 30. Value at Risk (VaR)

`VaR_α = −μ + z_α·σ` (parametric) — the loss not exceeded with probability `α` over a horizon

**What it does:** a single-number loss threshold. Widely required, deeply flawed.

**Use in practice:**
- Use the **empirical quantile** of actual (or bootstrapped) P&L, not the normal formula, unless
  you've verified normality — which for trading returns you generally can't.
- Prefer **CVaR / Expected Shortfall** = `E[loss | loss > VaR]`: it describes the tail's *size*,
  is **coherent/subadditive**, and is what actually threatens survival.
- Scale honestly: `√t` only under iid; horizons longer than a few days need block bootstrap.

**How it lies — know these before quoting a VaR:**
1. **Says nothing about severity beyond the threshold.** Two books with identical VaR can have
   wildly different worst cases. A strategy engineered to a VaR limit will hide risk *past* the
   quantile — that is the standard way VaR gets gamed.
2. **Not subadditive**: `VaR(A+B)` can exceed `VaR(A)+VaR(B)`. Aggregating VaR across sleeves is
   not sound; ES is.
3. **Normal-parametric VaR understates tails badly** (fat tails, §7).
4. **Estimated from history it has never seen the event that matters.** For a lending-heavy fund
   the dominant risk is a venue going to zero — a discrete, non-historical, non-normal event.
   Its VaR contribution is ~0 right until it's 100%. **VaR is structurally the wrong tool for
   protocol/counterparty risk; use scenario P&L (§9) and per-venue caps.**
5. **Backtest the VaR itself**: count exceedances. If a 95% daily VaR is breached far more (or
   less) than ~5% of days, the model is wrong — say so.

**Do this:** report empirical VaR **and** CVaR **and** the worst realised period, plus a scenario
table for the non-statistical risks. State plainly which binding risks (counterparty insolvency,
oracle failure, regulatory action) are **not** in any VaR: those are controlled structurally by
position caps, per-venue limits and timelocks — never by a risk number.
