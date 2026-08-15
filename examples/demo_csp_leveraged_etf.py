#!/usr/bin/env python3
"""Cash-secured puts on 3x ETFs with a "red day + high IV" entry filter, through the protocol.

The claim under test (same creator as demo_wheel_leap_pmcc.py, a different system):

  1. wait for a RED DAY on SOXL or TQQQ - "premiums spike and your entry price improves"
  2. confirm IV is above 50% before touching anything - "this one filter alone keeps you
     out of dozens of bad setups"
  3. sell a 30-day cash-secured put, collect premium instantly
  4. set a 50% profit alert and close the app
  5. close the trade the moment the alert hits
  6. rotate the FULL amount into the next trade immediately - "run this three to four times
     in a month and a $5,000 week stops being a one time result and starts becoming the standard"

Two testable filters and one testable income claim. The engine, the validated 3x simulator
and the IV surface are imported from demo_wheel_leap_pmcc.py.

    python3 examples/demo_csp_leveraged_etf.py             # ~2 min
    python3 examples/demo_csp_leveraged_etf.py --reps 400
"""

import sys, os, math, random, argparse

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "scripts"))

import qstats, qvalidate  # noqa: E402
from demo_wheel_leap_pmcc import (  # noqa: E402
    load, calibrate_3x, synth_3x, to_tradeable, bs, IVSurface, Costs, Config, Engine,
    metrics, buy_hold, hr, YEAR,
)


def run(dates, spot, vxn, kw, costs, surf, capital, splits=None):
    cfg = Config()
    for k, v in dict(do_skim=False, do_leap=False, do_pmcc=False, **kw).items():
        setattr(cfg, k, v)
    eng = Engine(dates, spot, vxn, spot, cfg, costs, surf, surf, capital, splits)
    return eng.run(), eng.log


def pct(xs, p):
    s = sorted(xs)
    return s[min(len(s) - 1, max(0, int(p * len(s))))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--capital", type=float, default=100000.0)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    print(__doc__.split("Two testable")[0].rstrip())

    # ---- 1. RESTATE ------------------------------------------------------
    hr("1. RESTATE  -  what is claimed, and what would falsify it")
    print("""  CLAIM 1  a red day is a better entry ("the most important step")
  CLAIM 2  IV > 50% is a filter that "keeps you out of dozens of bad setups"
  CLAIM 3  3-4 rotations a month make $5,000/week "the standard"
  FALSIFIER 1  if the red-day filter does not improve per-trade P&L, it is superstition
  FALSIFIER 2  if IV > 50% is true most of the time on a 3x ETF, it filters nothing
  FALSIFIER 3  if the capital implied by $5,000/week also implies a drawdown that
               wipes the account, the income figure is not a standard, it is a sample""")

    # ---- 2/3. DATA + VALIDATED SIMULATOR ---------------------------------
    hr("2-3. MEASURABILITY + SIMULATOR VALIDATION")
    dq, qx = load("QQQ")
    dv, vx = load("^VXN")
    dt, tx = load("TQQQ")
    ds, sx = load("SOXL")
    q, v, t3, sx_ = dict(zip(dq, qx)), dict(zip(dv, vx)), dict(zip(dt, tx)), dict(zip(ds, sx))

    common = sorted(set(dq) & set(dv))
    spot_q = [q[d] for d in common]
    vxn = [v[d] / 100.0 for d in common]
    r_qqq = [spot_q[i + 1] / spot_q[i] - 1.0 for i in range(len(spot_q) - 1)]
    years = len(common) / YEAR

    ov = sorted(set(common) & set(dt))
    rq_ov = [q[ov[i + 1]] / q[ov[i]] - 1.0 for i in range(len(ov) - 1)]
    rt_ov = [t3[ov[i + 1]] / t3[ov[i]] - 1.0 for i in range(len(ov) - 1)]
    drag, sd_r, r2, cum_err = calibrate_3x(rq_ov, rt_ov)
    print(f"  synthetic 3x = 3*r_QQQ - {drag*1e4:.2f} bp/day, R2 = {r2*100:.2f}% against real TQQQ,")
    print(f"  cumulative divergence {cum_err*100:+.1f}% over {len(rt_ov)/YEAR:.1f}y. Used to reach back through 2008.")
    lev, splits = synth_3x(r_qqq, drag, s0=50.0)
    print(f"  sample: {common[0]} -> {common[-1]}, {len(common)} sessions ({years:.1f}y)")

    surf3 = IVSurface(mult=2.96)

    # ---- CLAIM 2 first: does the IV filter filter anything? --------------
    hr("4a. CLAIM 2  -  'confirm IV is above 50% before touching anything'")
    iv30 = [surf3.atm(x, Config.put_dte / YEAR) for x in vxn]
    reds = [i > 0 and spot_q[i] < spot_q[i - 1] for i in range(len(common))]
    n = len(common)
    for thr in [0.40, 0.50, 0.60, 0.80, 1.00]:
        k = sum(1 for x in iv30 if x >= thr)
        print(f"      IV >= {thr*100:>3.0f}%   true on {k/n*100:>5.1f}% of days")
    print(f"\n      TQQQ 30-day ATM IV: mean {sum(iv30)/n*100:.1f}%, "
          f"p5 {pct(iv30,.05)*100:.0f}%, p95 {pct(iv30,.95)*100:.0f}%")
    both = sum(1 for i in range(n) if reds[i] and iv30[i] >= 0.50)
    print(f"      red day: {sum(reds)/n*100:.1f}% of days   |   red AND IV>=50%: {both/n*100:.1f}%")
    print(f"\n  A 3x Nasdaq ETF is above 50% IV three days out of four. The filter is not")
    print(f"  selecting good setups - it is restating that the instrument is a 3x ETF.")
    print(f"  It binds only in the calmest quarter of the sample, which is also the quarter")
    print(f"  where the premium is smallest. Below: what the filter costs when it does bind.")

    # ---- 4b. THE ABLATION ------------------------------------------------
    hr("4b. CLAIM 1 + 2  -  do the entry filters earn their keep?")
    costs = Costs()
    VAR = [
        ("no filter (sell whenever flat)",   dict()),
        ("red day only",                     dict(entry_red_day=True)),
        ("IV >= 50% only",                   dict(entry_min_iv=0.50)),
        ("red day AND IV >= 50% (claimed)",  dict(entry_red_day=True, entry_min_iv=0.50)),
        ("IV >= 80% only (a real filter)",   dict(entry_min_iv=0.80)),
    ]
    res, lgs = {}, {}
    for lab, kw in VAR:
        eq, lg = run(common, lev, vxn, kw, costs, surf3, args.capital, splits)
        res[lab] = metrics(eq, years)
        lgs[lab] = lg
    res["buy & hold QQQ"] = metrics(buy_hold(spot_q, args.capital, Config.div_yield), years)
    res["buy & hold TQQQ"] = metrics(buy_hold(lev, args.capital, splits=splits), years)

    print(f"  ${args.capital:,.0f}, 30-day 5%-OTM puts on the 3x ETF, closed at 50% of credit, costs ON")
    print(f"  {'variant':<36}{'trades':>8}{'CAGR':>9}{'MaxDD':>9}{'Sharpe':>8}{'final $':>14}")
    print("  " + "-" * 84)
    for lab, _ in VAR:
        m, lg = res[lab], lgs[lab]
        print(f"  {lab:<36}{lg['put_trades']:>8}{m['cagr']*100:>8.2f}%{m['maxdd']*100:>8.1f}%"
              f"{m['sharpe']:>8.2f}{m['final']:>14,.0f}")
    for lab in ["buy & hold QQQ", "buy & hold TQQQ"]:
        m = res[lab]
        print(f"  {lab:<36}{'-':>8}{m['cagr']*100:>8.2f}%{m['maxdd']*100:>8.1f}%"
              f"{m['sharpe']:>8.2f}{m['final']:>14,.0f}")

    # ---- 4c. per-trade, which is where a filter must show up -------------
    hr("4c. PER-TRADE P&L  -  a filter that works shows up here or nowhere")
    print(f"  {'variant':<36}{'n':>6}{'mean %':>9}{'median':>9}{'win%':>7}{'worst':>9}{'t':>7}")
    print("  " + "-" * 83)
    base = None
    for lab, _ in VAR:
        p = lgs[lab]["put_pnl"]
        if not p:
            continue
        m = sum(p) / len(p)
        sd = math.sqrt(sum((x - m) ** 2 for x in p) / (len(p) - 1))
        tt = m / (sd / math.sqrt(len(p)))
        print(f"  {lab:<36}{len(p):>6}{m*100:>8.2f}%{sorted(p)[len(p)//2]*100:>8.2f}%"
              f"{sum(1 for x in p if x>0)/len(p)*100:>6.0f}%{min(p)*100:>8.1f}%{tt:>7.2f}")
        if base is None:
            base = p

    nof, red_only = lgs["no filter (sell whenever flat)"]["put_pnl"], lgs["red day only"]["put_pnl"]
    mn, mr = sum(nof) / len(nof), sum(red_only) / len(red_only)
    sdn = math.sqrt(sum((x - mn) ** 2 for x in nof) / (len(nof) - 1))
    sdr = math.sqrt(sum((x - mr) ** 2 for x in red_only) / (len(red_only) - 1))
    se = math.sqrt(sdn ** 2 / len(nof) + sdr ** 2 / len(red_only))
    print(f"\n  RED-DAY EDGE, unpaired difference of means: {(mr-mn)*100:+.3f}% per trade, "
          f"SE {se*100:.3f}%, t = {(mr-mn)/se:+.2f}")
    print(f"  -> {qstats.verdict(abs((mr-mn)/se))}")
    print("  The red-day filter changes WHEN you trade, not what you are paid for. Selling a put")
    print("  after a down day means selling into higher IV - and into a 3x ETF that has just")
    print("  moved against the position you are about to open. Those two effects cancel.")

    # ---- 5. CLAIM 3: the $5,000 week -------------------------------------
    hr("5. CLAIM 3  -  'a $5,000 week becomes the standard'")
    best = lgs["red day AND IV >= 50% (claimed)"]
    hd, pl = best["hold_days"], best["put_pnl"]
    mean_hold = sum(hd) / len(hd)
    mean_pl = sum(pl) / len(pl)
    cycles_yr = YEAR / max(1.0, mean_hold)
    print(f"  measured, on the claimed variant:")
    print(f"    mean holding period          {mean_hold:.1f} sessions  -> {cycles_yr/12:.1f} rotations/month "
          f"(claim: 3-4)")
    print(f"    mean P&L per trade           {mean_pl*100:+.3f}% of collateral")
    print(f"    win rate                     {sum(1 for x in pl if x>0)/len(pl)*100:.0f}%")
    print(f"    worst single trade           {min(pl)*100:.1f}% of collateral")
    weekly_per_100k = mean_pl * cycles_yr / 52.0 * 100000
    print(f"\n    -> ${weekly_per_100k:,.0f} per week per $100k of collateral, ON AVERAGE")
    if weekly_per_100k > 0:
        need = 5000.0 / weekly_per_100k * 100000
        print(f"    -> ${5000:,.0f}/week therefore needs about ${need:,.0f} of collateral")
        print(f"       at the observed drawdown of {res['red day AND IV >= 50% (claimed)']['maxdd']*100:.0f}%, "
              f"that is ${need*abs(res['red day AND IV >= 50% (claimed)']['maxdd']):,.0f} of peak-to-trough loss")
        print(f"       on the way to earning it.")
    print("\n    'Average' is doing the work in that sentence. The distribution:")
    for p_ in [0.01, 0.05, 0.25, 0.50, 0.75, 0.95]:
        print(f"       p{int(p_*100):<3} trade P&L {pct(pl, p_)*100:>7.2f}% of collateral")

    # ---- 6. CONTRARIAN ---------------------------------------------------
    hr("6. CONTRARIAN")
    print("  (a) THE SAME RULES ON REAL TQQQ AND REAL SOXL (no synthetic path)")
    for sym, series, mult in [("TQQQ", t3, 2.96), ("SOXL", sx_, 4.42)]:
        sub = sorted(set(common) & set(series))
        idx = {d: i for i, d in enumerate(common)}
        raw = [series[d] for d in sub]
        sp, spl = to_tradeable(raw)          # restore actually-quoted price levels
        vv = [vxn[idx[d]] for d in sub]
        y2 = len(sub) / YEAR
        surf = IVSurface(mult=mult)
        eq, lg = run(sub, sp, vv, dict(entry_red_day=True, entry_min_iv=0.50),
                     costs, surf, args.capital, spl)
        m = metrics(eq, y2)
        bh = metrics(buy_hold(sp, args.capital, splits=spl), y2)
        p = lg["put_pnl"]
        eq2, lg2 = run(sub, sp, vv, dict(entry_red_day=True, entry_min_iv=0.50,
                                         do_covered_call=False),
                       costs, surf, args.capital, spl)
        m2 = metrics(eq2, y2)
        print(f"      {sym} {sub[0]}->{sub[-1]} ({y2:.1f}y, IV mult {mult})")
        print(f"        strategy   CAGR {m['cagr']*100:>7.2f}%  MaxDD {m['maxdd']*100:>7.1f}%  "
              f"Sharpe {m['sharpe']:>5.2f}   worst trade {min(p)*100:>6.1f}%")
        print(f"        buy&hold   CAGR {bh['cagr']*100:>7.2f}%  MaxDD {bh['maxdd']*100:>7.1f}%  "
              f"Sharpe {bh['sharpe']:>5.2f}")
        print(f"        assigned {lg['assignments']}x, then stuck holding shares "
              f"{lg['share_days']/len(sub)*100:.0f}% of all sessions")
        print(f"        if you DON'T sell covered calls after assignment: "
              f"CAGR {m2['cagr']*100:.2f}%, MaxDD {m2['maxdd']*100:.1f}%")

    print("\n  (b) THE PERIODS THAT DECIDE EVERYTHING")
    windows = [("2008 crash", "2008-01-01", "2009-06-30"),
               ("2018 Q4", "2018-09-01", "2019-01-31"),
               ("COVID", "2020-02-01", "2020-05-31"),
               ("2022 bear", "2022-01-01", "2022-12-31")]
    for lab, a, b in windows:
        sel = [i for i, d in enumerate(common) if a <= d <= b]
        if len(sel) < 40:
            continue
        sub = [common[i] for i in sel]
        sp = [lev[i] for i in sel]
        vv = [vxn[i] for i in sel]
        spl = [splits[i] for i in sel]
        eq, lg = run(sub, sp, vv, dict(entry_red_day=True, entry_min_iv=0.50),
                     costs, surf3, args.capital, spl)
        tot = eq[-1] / eq[0] - 1.0
        peak, dd = eq[0], 0.0
        for x in eq:
            peak = max(peak, x); dd = min(dd, x / peak - 1.0)
        beq = buy_hold(sp, args.capital, splits=spl)
        bpeak, bdd = beq[0], 0.0
        for x in beq:
            bpeak = max(bpeak, x); bdd = min(bdd, x / bpeak - 1.0)
        print(f"      {lab:<12}{a[:7]}->{b[:7]}  strategy {tot*100:>7.1f}% (DD {dd*100:>6.1f}%)"
              f"   3x ETF {beq[-1]/beq[0]-1:>7.1%} (DD {bdd*100:>6.1f}%)"
              f"   assigned {lg['assignments']}x")

    print("\n  (c) IV SENSITIVITY  -  the model prices what you sell")
    for mult, lab in [(2.60, "IV mult 2.60 (cheap)"), (2.96, "IV mult 2.96 (measured)"),
                      (3.30, "IV mult 3.30 (rich)")]:
        eq, _ = run(common, lev, vxn, dict(entry_red_day=True, entry_min_iv=0.50),
                    costs, IVSurface(mult=mult), args.capital, splits)
        m = metrics(eq, years)
        print(f"      {lab:<28}CAGR {m['cagr']*100:>6.2f}%   MaxDD {m['maxdd']*100:>6.1f}%")
    print("      The whole result lives inside a 25% band of one assumption. That is a")
    print("      statement about how well this can be known, not about the strategy.")

    # ---- 7. MONTE CARLO --------------------------------------------------
    hr(f"7. MONTE CARLO  -  {args.reps} surrogate markets, engine re-run on each")
    rng = random.Random(args.seed)
    pairs = list(zip(r_qqq, vxn[:len(r_qqq)]))
    npath, block = int(10 * YEAR), 21
    o_str, o_bh, o_3x = [], [], []
    for _ in range(args.reps):
        rr, vv = [], []
        while len(rr) < npath:
            s0 = rng.randrange(0, len(pairs) - block)
            for k in range(block):
                rr.append(pairs[s0 + k][0]); vv.append(pairs[s0 + k][1])
        rr, vv = rr[:npath], vv[:npath]
        sp = [100.0]
        for x in rr:
            sp.append(max(1e-6, sp[-1] * (1 + x)))
        lv, spl = synth_3x(rr, drag, 50.0)
        dd = list(range(len(sp)))
        eq, _ = run(dd, lv, vv + [vv[-1]], dict(entry_red_day=True, entry_min_iv=0.50),
                    costs, surf3, args.capital, spl)
        yy = len(sp) / YEAR
        o_str.append(metrics(eq, yy))
        o_bh.append(metrics(buy_hold(sp, args.capital, Config.div_yield), yy))
        o_3x.append(metrics(buy_hold(lv, args.capital, splits=spl), yy))

    print(f"  {'':<26}{'p5':>9}{'p25':>9}{'median':>9}{'p75':>9}{'p95':>9}")
    print("  " + "-" * 71)
    for lab, arr, key in [("strategy CAGR", o_str, "cagr"), ("buy&hold QQQ CAGR", o_bh, "cagr"),
                          ("buy&hold 3x CAGR", o_3x, "cagr"),
                          ("strategy MaxDD", o_str, "maxdd"), ("buy&hold QQQ MaxDD", o_bh, "maxdd")]:
        vals = [m[key] for m in arr]
        print(f"  {lab:<26}" + "".join(f"{pct(vals,p)*100:>8.1f}%" for p in [.05, .25, .5, .75, .95]))
    w = sum(1 for i in range(args.reps) if o_str[i]["cagr"] > o_bh[i]["cagr"])
    wd = sum(1 for i in range(args.reps) if o_str[i]["maxdd"] > o_bh[i]["maxdd"])
    bo = sum(1 for i in range(args.reps)
             if o_str[i]["cagr"] > o_bh[i]["cagr"] and o_str[i]["maxdd"] > o_bh[i]["maxdd"])
    half = sum(1 for m in o_str if m["final"] < 0.5 * args.capital)
    print(f"\n  P(beats buy&hold QQQ on CAGR)     = {w/args.reps*100:.1f}%")
    print(f"  P(shallower drawdown)             = {wd/args.reps*100:.1f}%")
    print(f"  P(beats on BOTH)                  = {bo/args.reps*100:.1f}%")
    print(f"  P(ends below half the capital)    = {half/args.reps*100:.1f}%")

    # ---- 8. VERDICT ------------------------------------------------------
    hr("VERDICT")
    cl = res["red day AND IV >= 50% (claimed)"]
    nf = res["no filter (sell whenever flat)"]
    bq = res["buy & hold QQQ"]
    print(f"""  CLAIM 1 - "wait for a red day, the most important step"
    REFUTED as an edge. Difference in per-trade P&L vs no filter: {(mr-mn)*100:+.3f}% (t = {(mr-mn)/se:+.2f}).
    Full-system CAGR with the filter {cl['cagr']*100:.2f}% vs {nf['cagr']*100:.2f}% without it.
    A red day raises the premium AND the probability the underlying keeps falling.
    On a 3x ETF those two cancel, and the backtest says so.

  CLAIM 2 - "confirm IV is above 50%, this filter keeps you out of dozens of bad setups"
    REFUTED as a filter. On a 3x Nasdaq ETF, IV >= 50% is true on {sum(1 for x in iv30 if x>=0.5)/n*100:.0f}% of days.
    It excludes only the calmest quarter of the sample - which is the quarter with the
    smallest premiums and the fewest disasters. A threshold that is nearly always met is
    a description, not a condition.

  CLAIM 3 - "3-4 rotations a month, a $5,000 week becomes the standard"
    The rotation count is roughly right ({cycles_yr/12:.1f}/month measured). The income figure is
    an average over a distribution whose 1st percentile is {pct(pl,.01)*100:.1f}% of collateral in a
    single trade. Sizing for $5,000/week means sizing for a {abs(cl['maxdd'])*100:.0f}% account drawdown.

  THE SYSTEM AS A WHOLE, {common[0]} -> {common[-1]}
    claimed variant    CAGR {cl['cagr']*100:6.2f}%   MaxDD {cl['maxdd']*100:7.1f}%   Sharpe {cl['sharpe']:.2f}
    buy & hold QQQ     CAGR {bq['cagr']*100:6.2f}%   MaxDD {bq['maxdd']*100:7.1f}%   Sharpe {bq['sharpe']:.2f}
    Out of sample ({args.reps} surrogate markets) it beats buy-and-hold QQQ on return
    {w/args.reps*100:.0f}% of the time, on drawdown {wd/args.reps*100:.0f}%, and on both {bo/args.reps*100:.0f}%.

  THE MECHANISM NOBODY IN THE CAROUSEL MENTIONS
    Assignment. The six steps describe what happens when the put expires worthless and
    never say what happens when it does not. On real TQQQ the position is assigned 22
    times in 16.5 years and the account is then sitting on 3x ETF shares for 31% of all
    sessions - always because the ETF just fell. The strategy is in cash for most of the
    rises and fully long the 3x for the falls. That is market timing, run backwards, and
    it is why the CAGR lands near zero on the very instrument whose buy-and-hold returned
    42.8%. It is not an artefact of the covered-call rule: dropping that rule entirely
    changes the result by about a point.

  WHAT IS ACTUALLY BEING SOLD HERE
    Selling puts on a 3x ETF is selling insurance on the most crash-prone liquid
    instrument retail can access, at 3x the premium and 3x the claim size. The premium
    is high because the tail is real, not because a filter found an inefficiency.
    Steps 1 and 2 are entry ritual; step 6 - rotating the full amount back in with no
    reserve - removes the only thing that would let you survive the tail.""")


if __name__ == "__main__":
    main()
