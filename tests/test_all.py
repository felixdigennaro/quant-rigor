#!/usr/bin/env python3
"""quant-rigor test suite. Run: python3 tests/test_all.py

STDLIB ONLY, except one section that needs numpy: it runs when numpy is installed and is
listed as SKIP when it is not. Exits 1 if any check fails.

Every numeric assertion is against an INDEPENDENT reference (published statistical tables,
closed-form results, or analytically known values) — never against the library's own output.
That discipline is what caught a 1% p-value error in the t-tail at df=1.
"""
import math
import random
import os
import sys
import xml.dom.minidom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import qstats as S          # noqa: E402
import qvalidate as V       # noqa: E402
import qsize as Z           # noqa: E402
import qcharts as C         # noqa: E402
import qmicro as MI      # noqa: E402
import qportfolio as PF  # noqa: E402
import qclaim as QC      # noqa: E402

PASS, FAIL, SKIP = [], [], []


def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(f"{name}{(' — ' + detail) if detail else ''}")


def close(name, got, want, tol, detail=""):
    d = abs(got - want)
    ok(name, d <= tol, detail or f"got {got:.8g}, want {want:.8g}, |Δ|={d:.2e} (tol {tol:.1e})")


# ══════════════════════════════════ 0. the suite can fail
# A copy of this file printed SUITE RED and still exited 0 (found 2026-10-25), so a CI check
# on it could never go red. A check written below the exit never runs. Both break one rule:
# the last statement of this file is the exit.
import ast  # noqa: E402
_last = ast.unparse(ast.parse(open(os.path.abspath(__file__), encoding="utf-8").read()).body[-1])
if _last != "sys.exit(1 if FAIL else 0)":
    sys.exit(f"tests/test_all.py must end with `sys.exit(1 if FAIL else 0)`, "
             f"not `{_last.splitlines()[0]}`")


# ══════════════════════════════════ 1. t distribution vs published tables
# Source: standard Student-t critical value tables (two-sided).
T_TABLE = {(1, .05): 12.706, (1, .10): 6.314, (2, .05): 4.303, (3, .05): 3.182,
           (3, .01): 5.841, (5, .05): 2.571, (5, .01): 4.032, (10, .05): 2.228,
           (14, .05): 2.145, (14, .01): 2.977, (20, .05): 2.086, (30, .05): 2.042,
           (30, .01): 2.750, (60, .05): 2.000, (120, .05): 1.980, (100, .001): 3.390}
for (df, alpha), tc in T_TABLE.items():
    close(f"t_sf df={df} α={alpha}", 2 * S.t_sf(tc, df), alpha, 1e-3)

ok("t_sf t=0 gives p=1", all(abs(2 * S.t_sf(0, df) - 1.0) < 1e-6 for df in (1, 2, 5, 30)))
ok("t_sf monotone in t", S.t_sf(1, 10) > S.t_sf(2, 10) > S.t_sf(3, 10))
close("t_sf df=200 ≈ normal", 2 * S.t_sf(1.959964, 200), 0.0512, 2e-3)
# large df must not overflow (math.gamma dies above ~340 → any multi-year daily series)
for df in (300, 755, 1500, 5000, 100000):
    p = 2 * S.t_sf(2.0, df)
    ok(f"t_sf df={df} finite (no gamma overflow)", math.isfinite(p) and 0 < p < 1, f"p={p}")
close("t_sf df=100000 → normal limit", 2 * S.t_sf(2.0, 100000), 0.0455003, 1e-5)
ok("t_sf monotone in df at fixed t", 2*S.t_sf(2.0, 10) > 2*S.t_sf(2.0, 100) > 2*S.t_sf(2.0, 5000))
# a long daily-style series must produce a usable p-value end to end
long_series = [0.0003 + 0.01 * math.sin(i * 0.7) for i in range(756)]
dl = S.describe(long_series)
ok("describe on 756-point series", math.isfinite(dl["p"]) and math.isfinite(dl["t"]),
   f"t={dl['t']:.2f} p={dl['p']:.4f}")

# ══════════════════════════════════ 2. normal helpers vs known quantiles
for p, z in [(0.975, 1.959964), (0.995, 2.575829), (0.95, 1.644854),
             (0.999, 3.090232), (0.5, 0.0), (0.025, -1.959964)]:
    close(f"norm_ppf({p})", V.norm_ppf(p), z, 1e-6)
close("norm_cdf(0)", V.norm_cdf(0), 0.5, 1e-12)
close("norm_cdf(1.959964)", V.norm_cdf(1.959964), 0.975, 1e-6)
ok("norm_ppf/cdf round-trip", all(abs(V.norm_cdf(V.norm_ppf(p)) - p) < 1e-9
                                  for p in (0.01, 0.3, 0.5, 0.87, 0.999)))

# ══════════════════════════════════ 3. describe() vs hand-computed
x = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]        # mean 5, population sd 2, sample sd 2.13809
d = S.describe(x)
close("describe mean", d["mean"], 5.0, 1e-12)
close("describe sd (sample, n-1)", d["sd"], math.sqrt(32/7), 1e-12)
close("describe SE", d["se"], math.sqrt(32/7) / math.sqrt(8), 1e-12)
close("describe t", d["t"], 5.0 / (math.sqrt(32/7) / math.sqrt(8)), 1e-12)
ok("describe win_rate", abs(d["win_rate"] - 100.0) < 1e-9)
ok("describe worst/best", d["worst"] == 2.0 and d["best"] == 9.0)
# symmetric data => skew 0
close("skew of symmetric data", S.describe([-2, -1, 0, 1, 2])["skew"], 0.0, 1e-12)
ok("skew sign (right tail)", S.describe([1, 1, 1, 1, 10])["skew"] > 0)

# ══════════════════════════════════ 4. paired test vs hand-computed
a = [1.0, 2.0, 3.0, 4.0]
b = [2.0, 4.0, 4.0, 6.0]                              # diffs 1,2,1,2 → mean 1.5, sd 0.57735
pr = S.paired(a, b)
close("paired mean", pr["mean"], 1.5, 1e-12)
close("paired sd", pr["sd"], math.sqrt(1/3), 1e-12)
close("paired t", pr["t"], 1.5 / (math.sqrt(1/3) / 2), 1e-12)
ok("paired counts", pr["improved"] == 4 and pr["worsened"] == 0)
try:
    S.paired([1, 2], [1, 2, 3]); ok("paired length guard", False)
except ValueError:
    ok("paired length guard", True)

# ══════════════════════════════════ 5. bootstrap / jackknife determinism & sanity
b1 = S.bootstrap(x, reps=500, seed=42)
b2 = S.bootstrap(x, reps=500, seed=42)
ok("bootstrap deterministic with seed", b1 == b2)
ok("bootstrap band brackets point estimate", b1["p05"] <= sum(x) <= b1["p95"],
   f"{b1['p05']:.1f} <= {sum(x)} <= {b1['p95']:.1f}")
ok("bootstrap block mode runs", S.bootstrap(x, reps=100, seed=1, block=3)["p50"] is not None)
jk = S.jackknife(x)
close("jackknife full == sum", jk["full"], sum(x), 1e-12)
ok("jackknife sign survives all", jk["sign_survives"] == f"{len(x)}/{len(x)}")
# a series whose sign flips when one point is dropped
jk2 = S.jackknife([10.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0])
ok("jackknife detects fragility", jk2["sign_survives"] != "10/10", jk2["sign_survives"])

# ══════════════════════════════════ 6. concentration / drawdown / VaR analytics
cc = S.concentration([10.0, 1.0, 1.0, 1.0], top=1)
close("concentration top1 share", cc["share_pct"], 100 * 10 / 13, 1e-9)
eq = S.equity_from_returns([1.0, -3.0, 1.0])           # equity 1, -2, -1 → peak 1 → max dd -3
dd = S.drawdown(eq)
close("max drawdown", dd["max_dd"], -3.0, 1e-12)
ok("periods under water", dd["periods_under_water"] == 2, str(dd))
# VaR/CVaR on a known ordered set: 20 values 1..20, alpha .95 → k=1 → VaR=2nd smallest? verify def
vc = S.var_cvar([float(i) for i in range(1, 21)], alpha=0.95)
ok("VaR is a low quantile", vc["VaR_empirical"] <= 2.0, str(vc["VaR_empirical"]))
ok("CVaR <= VaR", vc["CVaR_expected_shortfall"] <= vc["VaR_empirical"])
ok("worst is the min", vc["worst"] == 1.0)
close("kelly f* = mu/sigma^2", S.kelly_fraction(0.02, 0.10), 0.02 / 0.01, 1e-12)

# ══════════════════════════════════ 7. measurability gate (analytic)
m = S.measurability([10.0] * 9 + [1000.0], 40.0)       # median 10 → R = 4
close("measurability R", m["ratio_R"], 4.0, 1e-12)
ok("measurability verdict MEASURABLE", m["verdict"].startswith("MEASURABLE"), m["verdict"])
m2 = S.measurability([67.0, 67.0, 67.0], 20.0)         # R = 0.298
ok("measurability verdict FICTION", m2["verdict"].startswith("FICTION"), m2["verdict"])
m3 = S.measurability([40.0] * 5, 60.0)                 # R = 1.5
ok("measurability verdict NOT MEASURABLE", m3["verdict"].startswith("NOT"), m3["verdict"])
close("needed distance for R=3", m2["needed_distance_for_R3"], 201.0, 1e-9)

# ══════════════════════════════════ 8. search headroom & required n (closed form)
close("headroom = sd*sqrt(n)*sqrt(2 ln N)", S.headroom(10.0, 25, 100),
      10.0 * 5 * math.sqrt(2 * math.log(100)), 1e-9)
ok("headroom 0 for single variant", S.headroom(10.0, 25, 1) == 0.0)
close("required_n = (t*sd/eff)^2", S.required_n(287.0, 52.0, 2.0), (2 * 287 / 52) ** 2, 1e-9)
ok("required_n infinite for zero effect", S.required_n(10, 0) == float('inf'))


# ══════════════════ 8b. la soglia ESATTA della ricerca, e la deriva di sqrt(2 ln N)
#
# Nato da un errore reale: sqrt(2 ln N) era usato come se fosse "la soglia di
# significativita' dopo aver cercato N volte". Non lo e'. Non e' nemmeno la MEDIA del
# massimo, ed e' un quantile che DERIVA con N — 22,5% a N=2, 15,4% a N=8, 8,5% a N=1e4.
# Chi lo usa come test al 5% sta usando uno screening al 10-20%.
close("search_threshold(8, 5%) = 2.490", S.search_threshold(8, 0.05), 2.490, 0.003)
close("search_threshold(462, 5%) = 3.693", S.search_threshold(462, 0.05), 3.693, 0.003)
close("search_threshold(462, 1%) = 4.088", S.search_threshold(462, 0.01), 4.088, 0.003)
close("search_threshold(462, 10%) = 3.505", S.search_threshold(462, 0.10), 3.505, 0.003)
close("search_quantile(8) = 15.4%", S.search_quantile(8), 0.154, 0.002)
close("search_quantile(462) = 10.1%", S.search_quantile(462), 0.101, 0.002)
close("search_quantile(10000) = 8.5%", S.search_quantile(10000), 0.085, 0.002)
ok("il quantile di sqrt(2 ln N) DERIVA con N (non e' fisso)",
   S.search_quantile(2) > S.search_quantile(8) > S.search_quantile(462) > S.search_quantile(10000))
ok("sqrt(2 ln N) sta SOPRA la media del massimo e SOTTO la soglia al 5%",
   all(S.search_threshold(n, 0.05) > math.sqrt(2*math.log(n)) for n in (8, 100, 462, 10000)))
ok("soglia piu' severa = t piu' alto", S.search_threshold(462, 0.01) > S.search_threshold(462, 0.05)
   > S.search_threshold(462, 0.10))
ok("stabile nella coda profonda (niente cancellazione a t alto)",
   0.0 < S._max_exceed_prob(6.0, 462) < 1e-6)

# ══════════════════════════════════ 9. HAC / Newey-West
iid = [1.0, -1.0] * 40                                  # strongly negatively autocorrelated
h = V.hac_t(iid)
ok("HAC runs on alternating series", math.isfinite(h["t_hac"]))
ok("HAC deflates SE for negative autocorr", h["se_inflation"] < 1.0, f"{h['se_inflation']:.3f}")
persistent = []
v = 0.0
for i in range(200):                                    # AR(1) with block innovations -> genuinely positive autocorrelation
    v = 0.85 * v + (0.3 if (i // 7) % 2 == 0 else -0.3)
    persistent.append(v + 0.05)
hp = V.hac_t(persistent)
ok("HAC inflates SE for positive autocorr", hp["se_inflation"] > 1.0, f"{hp['se_inflation']:.3f}")
ok("effective n <= n", V.effective_sample_size(persistent) <= len(persistent))
close("effective n of iid ≈ n", V.effective_sample_size([1.0, 2.0, 1.5, 1.8, 1.2] * 20), 100, 35)

# ══════════════════════════════════ 10. Deflated Sharpe / MinTRL
sr0_1 = V.expected_max_sharpe(10, 0.5)
sr0_2 = V.expected_max_sharpe(1000, 0.5)
ok("SR0 grows with number of trials", sr0_2 > sr0_1, f"{sr0_1:.3f} -> {sr0_2:.3f}")
close("SR0 scales linearly in sd", V.expected_max_sharpe(100, 1.0),
      2 * V.expected_max_sharpe(100, 0.5), 1e-9)
d_hi = V.deflated_sharpe(1.2, 500, 10, 0.2)
d_lo = V.deflated_sharpe(0.15, 20, 5000, 0.5)
ok("DSR high for strong SR / few trials", d_hi["dsr"] > 0.95, f"{d_hi['dsr']:.3f}")
ok("DSR low for weak SR / many trials", d_lo["dsr"] < 0.5, f"{d_lo['dsr']:.3f}")
ok("DSR is a probability", 0.0 <= d_hi["dsr"] <= 1.0)
mt = V.min_track_record_length(0.5, 0.0, 3.0)
close("MinTRL closed form", mt,
      1 + (1 + (3-1)/4*0.5**2) * (V.norm_ppf(0.95) / 0.5) ** 2, 1e-9)
ok("MinTRL infinite when SR <= benchmark", V.min_track_record_length(0.0) == float('inf'))
ok("MinTRL falls as SR rises", V.min_track_record_length(1.0) < V.min_track_record_length(0.3))

# ══════════════════════════════════ 11. splits: walk-forward & purged k-fold
wf = V.walk_forward_splits(100, n_folds=4, anchored=True)
ok("walk-forward count", len(wf) == 4, str(len(wf)))
ok("walk-forward train precedes test", all(tr[1] == te[0] for tr, te in wf))
ok("walk-forward anchored starts at 0", all(tr[0] == 0 for tr, _ in wf))
wfr = V.walk_forward_splits(100, n_folds=4, anchored=False)
ok("rolling train is bounded", all(tr[0] > 0 for tr, _ in wfr[1:]))
ok("no test beyond n", all(te[1] <= 100 for _, te in wf))

pk = V.purged_kfold_splits(50, n_folds=5, label_span=1, embargo=0)
ok("purged kfold count", len(pk) == 5)
ok("purged kfold covers all test indices once",
   sorted(i for _, te in pk for i in te) == list(range(50)))
ok("no train/test overlap", all(not (set(tr) & set(te)) for tr, te in pk))
pk2 = V.purged_kfold_splits(50, n_folds=5, label_span=5, embargo=3)
ok("purging+embargo shrinks train", len(pk2[2][0]) < len(pk[2][0]),
   f"{len(pk2[2][0])} < {len(pk[2][0])}")
ok("embargo removes the bars right after test",
   all(i not in pk2[0][0] for i in range(10, 13)))

# ══════════════════════════════════ 12. PBO
# Case A: one genuinely better config → PBO should be low.
good = [[1.0 + 0.1 * ((i * 7) % 5) for i in range(40)],          # best in every period
        [0.2 + 0.1 * ((i * 3) % 5) for i in range(40)],
        [0.1 + 0.1 * ((i * 5) % 5) for i in range(40)]]
pa = V.pbo(good, n_splits=6, max_combos=40)
ok("PBO low when one config truly dominates", pa["pbo"] < 0.25, f"{pa['pbo']:.2f}")
# Case B: pure noise, deterministic pseudo-random, no config is better → PBO should be high.
noise = [[math.sin(c * 2.7 + t * 1.3) for t in range(40)] for c in range(8)]
pb = V.pbo(noise, n_splits=6, max_combos=40)
ok("PBO higher on noise than on a real edge", pb["pbo"] > pa["pbo"], f"{pb['pbo']:.2f} > {pa['pbo']:.2f}")
ok("PBO in [0,1]", 0.0 <= pb["pbo"] <= 1.0)
ok("PBO guards tiny input", math.isnan(V.pbo([[1, 2, 3]])["pbo"]))

# ══════════════════════════════════ 13. FDR / Bonferroni
pv = [0.001, 0.008, 0.02, 0.04, 0.2, 0.5, 0.7, 0.9]
bh = V.benjamini_hochberg(pv, q=0.10)
# BH by hand: thresholds i/m*q = .0125,.025,.0375,.05,... → largest i with p_(i)<=thr is i=2
ok("BH selects the right set", bh["significant"] == [0, 1, 2, 3], str(bh["significant"]))
close("BH cutoff", bh["cutoff"], 4 / 8 * 0.10, 1e-12)
bf = V.bonferroni(pv, alpha=0.05)
ok("Bonferroni stricter than BH", len(bf["significant"]) <= len(bh["significant"]))
close("Bonferroni cutoff", bf["cutoff"], 0.05 / 8, 1e-12)
ok("BH empty input safe", V.benjamini_hochberg([])["significant"] == [])

# ══════════════════════════════════ 14. sizing
close("kelly_continuous", Z.kelly_continuous(0.01, 0.1), 1.0, 1e-12)
# discrete Kelly: p=.6, win=1, loss=1 → f* = .6/1 - .4/1 = .2
close("kelly_discrete", Z.kelly_discrete(0.6, 1.0, 1.0), 0.2, 1e-12)
close("fractional kelly", Z.fractional_kelly(1.0, 0.25), 0.25, 1e-12)
# growth is zero at twice Kelly
f_star = Z.kelly_continuous(0.01, 0.1)
close("growth(2f*) = 0", Z.growth_rate(2 * f_star, 0.01, 0.1), 0.0, 1e-12)
ok("growth(f*) > growth(2f*)", Z.growth_rate(f_star, 0.01, 0.1) > Z.growth_rate(2 * f_star, 0.01, 0.1))
p_dd = Z.prob_hit_drawdown(0.2, 1.0, 0.01, 0.1)
close("prob_hit_drawdown closed form", p_dd, math.exp(-2 * 0.01 * 0.2 / 0.01), 1e-12)
ok("prob_hit_drawdown falls with bigger dd",
   Z.prob_hit_drawdown(0.4, 1.0, 0.01, 0.1) < Z.prob_hit_drawdown(0.1, 1.0, 0.01, 0.1))
ok("prob_hit = 1 when mean <= 0", Z.prob_hit_drawdown(0.2, 1.0, -0.01, 0.1) == 1.0)
fdd = Z.dd_constrained_fraction(0.2, 0.01, 0.1, prob=0.10)
close("dd_constrained inverts the barrier",
      Z.prob_hit_drawdown(0.2, fdd, 0.01, 0.1), 0.10, 1e-9)
ss = Z.sizing_summary(0.01, 0.1, dd_tolerated=0.2)
ok("sizing recommends the min of the two", ss["recommended_f"] ==
   min(ss["kelly_fractional"], ss["dd_constrained"]))
ok("binding constraint labelled", ss["binding_constraint"] in ("drawdown", "Kelly fraction"))

# ══════════════════════════════════ 15. vol targeting, impact, capacity, costs
close("ewma_vol of constant-magnitude series", Z.ewma_vol([0.01, -0.01] * 50, lam=0.94),
      0.01, 2e-3)
close("vol_target_size", Z.vol_target_size(0.10, 0.20), 0.5, 1e-12)
ok("vol_target cap applies", Z.vol_target_size(0.10, 0.001, cap=3.0) == 3.0)
ok("vol_target zero on zero vol", Z.vol_target_size(0.1, 0.0) == 0.0)
# square-root law: quadrupling size doubles impact
i1 = Z.impact_sqrt(1_000, 1_000_000, 0.02)
i4 = Z.impact_sqrt(4_000, 1_000_000, 0.02)
close("impact scales as sqrt(size)", i4 / i1, 2.0, 1e-9)
# capacity is the inverse of impact
cap = Z.capacity(alpha_per_unit=i1, daily_volume=1_000_000, sigma_daily=0.02)
close("capacity inverts impact", cap, 1000.0, 1e-6)
td = Z.turnover_drag(0.10, 12, 0.002)                  # cost .024, net .076, surv 4.17
close("turnover net", td["net"], 0.076, 1e-12)
close("turnover survivability", td["survivability"], 0.10 / 0.024, 1e-9)
ok("turnover verdict acceptable", td["verdict"] == "acceptable")
td2 = Z.turnover_drag(0.03, 20, 0.002)                 # cost .04 > gross → fragile
ok("turnover flags fragile", td2["verdict"].startswith("fragile"), td2["verdict"])
cr = Z.cost_ratio(7.25, 40)
close("cost ratio pct", cr["ratio_pct"], 100 * 7.25 / 40, 1e-9)
ok("cost ratio verdict material/dominant", "material" in cr["verdict"] or "dominant" in cr["verdict"],
   cr["verdict"])
er = Z.effective_rule_distance(40, 15, 0.8)
close("effective rule distance", er["effective"], 40 + 12, 1e-12)
af = Z.adverse_fill_adjustment(0.9, 0.5)
close("adverse fill adjustment", af["adverse_adjusted"], 0.45, 1e-12)

# ══════════════════════════════════ 16. charts: valid SVG, inside viewBox, degenerate inputs
def svg_checks(name, svg):
    try:
        xml.dom.minidom.parseString(svg)
        ok(f"{name}: valid XML", True)
    except Exception as e:
        ok(f"{name}: valid XML", False, str(e)[:60]); return
    import re
    bad = []
    for m in re.finditer(r'(?:x|cx)="(-?\d+\.?\d*)"', svg):
        if not -2 <= float(m.group(1)) <= C.W + 2:
            bad.append(m.group(1))
    for m in re.finditer(r'(?:y|cy)="(-?\d+\.?\d*)"', svg):
        if not -2 <= float(m.group(1)) <= C.H + 2:
            bad.append(m.group(1))
    for m in re.finditer(r'points="([^"]+)"', svg):
        for pt in m.group(1).split():
            if "," in pt:
                px, py = pt.split(",")
                if not (-2 <= float(px) <= C.W + 2 and -2 <= float(py) <= C.H + 2):
                    bad.append(pt)
    ok(f"{name}: inside viewBox", not bad, str(bad[:3]))
    ok(f"{name}: no NaN/inf in output", "nan" not in svg.lower() and "inf" not in svg.lower())


series = [-161.8, -321.8, 410.2, 96.8, -45.8, 233.5]
svg_checks("equity", C.equity_chart({"base": series, "variant": [v * 1.2 for v in series]}))
svg_checks("drawdown", C.drawdown_chart({"base": series}))
svg_checks("contribution", C.contribution_chart(series))
svg_checks("distribution", C.dist_chart(series))
svg_checks("surface", C.surface_chart([0, 25, 50, 75], [10, 20, 15, 12], headroom=8))
# degenerate: constant series, single point, all zeros, all negative
svg_checks("equity(constant)", C.equity_chart({"flat": [0.0] * 5}))
svg_checks("contribution(single)", C.contribution_chart([1.0]))
svg_checks("distribution(constant)", C.dist_chart([2.0] * 6))
svg_checks("drawdown(all negative)", C.drawdown_chart({"neg": [-1.0] * 5}))
svg_checks("surface(flat)", C.surface_chart([1, 2, 3], [5, 5, 5]))
ok("sparkline length sane", 0 < len(C.sparkline(series)) <= 48)
ok("sparkline empty safe", C.sparkline([]) == "")

# ══════════════════════════════════ 16b. report generator
import qreport as Q                                   # noqa: E402
rep = (Q.Report("Test report", subtitle="unit test")
       .assumptions(data="synthetic", periods=6, resolution="1 bar",
                    costs="fee 1.0", not_modelled=["funding", "gaps"],
                    extra={"Zone": "40 points"})
       .measurability([10.0] * 5, {"stop 40": 40.0, "trail 5": 5.0})
       .add_series("base", series)
       .add_series("variant", [v * 1.1 for v in series])
       .results()
       .compare("base", "variant", n_variants=9)
       .contraction_check("base", "variant")
       .overfitting(series_name="base", perf_matrix=[series, [v * 1.1 for v in series],
                                                     [v * 0.9 for v in series]], n_variants=9)
       .sizing(0.01, 0.1, dd_tolerated=0.15)
       .charts(param_grid=([1, 2, 3], [10, 20, 15]), headroom=5)
       .verdict("not decidable", "test reason", path_to_decide="collect more",
                does_not_support=["anything at all"]))
h = rep.to_html()
ok("report: valid single-root HTML", h.startswith("<!doctype html>") and h.rstrip().endswith("</html>"))
ok("report: all 8 sections present", all(f"{i} · " in h for i in range(1, 9)),
   [i for i in range(1, 9) if f"{i} · " not in h])
ok("report: SVG figures inlined", h.count("<figure>") == 5, str(h.count("<figure>")))
ok("report: zero external dependencies",
   "src=" not in h and "@import" not in h and "<script" not in h and 'href="http' not in h)
ok("report: theme-aware + print styles", "prefers-color-scheme" in h and "@media print" in h)
ok("report: disclaimer present", "not investment advice" in h.lower())
ok("report: verdict rendered", "not decidable" in h and "Change nothing" in h)
ok("report: limits section rendered", "does not support" in h)
ok("report: no NaN leaked as text", "nan" not in h.lower().replace("financial", ""))
ok("report: no promotional wording", not any(w in h.lower() for w in
   ("amazing", "guaranteed", "revolutionary", "outperform", "superior", "unlock", "maximize your")))
ok("report: HTML-escapes user text", "&lt;b&gt;" in Q.Report("<b>x</b>").to_html())
md = rep.to_markdown()
ok("report: markdown sibling", md.startswith("# Test report") and "Verdict" in md)
import tempfile as _tf
_p = os.path.join(_tf.mkdtemp(), "r.html")
ok("report: writes to disk", rep.write(_p) == _p and os.path.getsize(_p) > 4000)
# a report with only a verdict must still be valid
ok("report: minimal report valid",
   Q.Report("m").verdict("decidable-negative", "r").to_html().count("</html>") == 1)
xml_ok = True
try:                                                  # the body must be well-formed XML too
    import re as _re
    frag = _re.search(r"<body>(.*)</body>", h, _re.S).group(1)
    xml.dom.minidom.parseString("<root>" + frag.replace("&nbsp;", " ") + "</root>")
except Exception as _e:
    xml_ok = False
ok("report: body is well-formed markup", xml_ok)

# ══════════════════════════════════ 16c. portfolio (formulas 9-15) vs analytic results
import qportfolio as P                                # noqa: E402
# two assets, known correlation: build returns with rho = 0 exactly by construction
xa = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
xb = [1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0]     # orthogonal to xa
nm, R = P.corr_matrix({"a": xa, "b": xb})
close("corr of orthogonal series = 0", R[0][1], 0.0, 1e-12)
close("corr diagonal = 1", R[0][0], 1.0, 1e-12)
nm, Cv = P.cov_matrix({"a": xa, "b": xb})
close("cov of unit-magnitude alternating series", Cv[0][0], 8 / 7, 1e-12)
# perfectly correlated pair -> corr 1, PC1 share 1, infinite condition number
nm, R2 = P.corr_matrix({"a": xa, "b": [2 * v for v in xa]})
close("corr of proportional series = 1", R2[0][1], 1.0, 1e-9)
_, C2 = P.cov_matrix({"a": xa, "b": [2 * v for v in xa]})
close("PC1 share = 1 for collinear assets", P.pc1_share(C2), 1.0, 1e-6)
ok("condition number huge for collinear assets", P.condition_number(C2) > 1e6,
   f"{P.condition_number(C2):.2e}")
# portfolio vol: independent equal-variance assets, equal weights -> sigma/sqrt(N)
I2 = [[1.0, 0.0], [0.0, 1.0]]
close("portfolio vol of 2 independent unit-var assets", P.portfolio_vol([0.5, 0.5], I2),
      1 / math.sqrt(2), 1e-12)
I4 = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
close("portfolio vol of 4 independent, equal weight", P.portfolio_vol([0.25] * 4, I4), 0.5, 1e-12)
# risk contributions must sum to portfolio vol (Euler decomposition)
cov3 = [[0.04, 0.01, 0.00], [0.01, 0.09, 0.02], [0.00, 0.02, 0.16]]
rc, sp = P.marginal_risk_contribution([0.4, 0.4, 0.2], cov3)
close("risk contributions sum to sigma_p", sum(rc), sp, 1e-12)
ok("risk contributions all positive here", all(x > 0 for x in rc), str([round(x, 5) for x in rc]))
# risk parity: contributions equal, and the highest-vol asset gets the lowest weight
rp = P.risk_parity_weights(cov3)
rcp, spp = P.marginal_risk_contribution(rp, cov3)
ok("risk parity equalises contributions", max(rcp) - min(rcp) < 1e-6 * max(1.0, spp),
   str([round(x, 6) for x in rcp]))
close("risk parity weights sum to 1", sum(rp), 1.0, 1e-9)
ok("risk parity underweights the riskiest asset", rp[2] < rp[0], str([round(x, 3) for x in rp]))
# eigen: diagonal matrix has its diagonal as eigenvalues, descending
vals, vecs = P.eigen_symmetric([[3.0, 0.0], [0.0, 1.0]])
close("eigenvalue 1 of diagonal matrix", vals[0], 3.0, 1e-10)
close("eigenvalue 2 of diagonal matrix", vals[1], 1.0, 1e-10)
# eigen: known 2x2 [[2,1],[1,2]] -> eigenvalues 3 and 1
vals2, _ = P.eigen_symmetric([[2.0, 1.0], [1.0, 2.0]])
close("eigenvalues of [[2,1],[1,2]]", vals2[0], 3.0, 1e-9)
close("second eigenvalue of [[2,1],[1,2]]", vals2[1], 1.0, 1e-9)
# trace is preserved by eigendecomposition
close("trace preserved", sum(P.eigen_symmetric(cov3)[0]), sum(cov3[i][i] for i in range(3)), 1e-9)
# Marchenko-Pastur edge grows with N/T
ok("MP edge grows as T shrinks", P.marchenko_pastur_edge(10, 20) > P.marchenko_pastur_edge(10, 1000))
close("MP edge at q=1", P.marchenko_pastur_edge(10, 10, 1.0), 4.0, 1e-12)
# shrinkage moves towards the identity and lowers the condition number
sh, delta = P.shrink_cov(C2, n_obs=8)
ok("shrinkage reduces condition number", P.condition_number(sh) < P.condition_number(C2),
   f"{P.condition_number(sh):.2e} < {P.condition_number(C2):.2e}")
ok("shrinkage delta in (0,1]", 0 < delta <= 1)
close("shrinkage preserves the trace", sum(sh[i][i] for i in range(2)),
      sum(C2[i][i] for i in range(2)), 1e-9)
# min-variance with a cap: cap must bind and be respected
mv, mvvol = P.solve_min_variance(cov3, cap=0.5)
ok("min-variance weights respect the cap", all(x <= 0.5 + 1e-6 for x in mv), str([round(x,3) for x in mv]))
close("min-variance weights sum to 1", sum(mv), 1.0, 1e-6)
ok("min-variance vol <= equal-weight vol", mvvol <= P.portfolio_vol([1/3]*3, cov3) + 1e-9,
   f"{mvvol:.5f} vs {P.portfolio_vol([1/3]*3, cov3):.5f}")
ok("binding caps detected", isinstance(P.binding_caps(mv, 0.5), list))
# tail-conditional correlation runs and returns a proper matrix
_, Rt = P.tail_conditional_corr({"a": xa, "b": xb, "c": [v * 0.5 for v in xa]})
ok("tail correlation is a square matrix with unit diagonal",
   len(Rt) == 3 and all(abs(Rt[i][i] - 1.0) < 1e-9 for i in range(3)))
ok("portfolio report runs", P.report({"a": xa, "b": xb}, cap=0.6).startswith("==="))
# tail dependence: analytic benchmarks — independent -> lambda ~ quantile; identical -> 1
import random as _rnd
_r = _rnd.Random(11)
_ind = {f"s{i}": [_r.gauss(0, 1) for _ in range(600)] for i in range(3)}
_, _, lam_i, q_i = P.tail_dependence(_ind, 0.10)
close("tail dependence of independent series ~ quantile", lam_i, q_i, 0.06)
_a = [_r.gauss(0, 1) for _ in range(400)]
_, _, lam_s, _ = P.tail_dependence({"a": _a, "b": [3 * v for v in _a]}, 0.10)
close("tail dependence of proportional series = 1", lam_s, 1.0, 1e-9)
_, Lm, _, _ = P.tail_dependence({"a": _a, "b": [3 * v for v in _a]}, 0.10)
ok("tail dependence diagonal is 1", all(abs(Lm[i][i] - 1.0) < 1e-12 for i in range(2)))
ok("tail dependence in [0,1]", all(0 <= Lm[i][j] <= 1 for i in range(2) for j in range(2)))
# lambda must rise when a common stress factor is added
_base = {f"s{i}": [_r.gauss(0, 1) for _ in range(600)] for i in range(3)}
_shock = [(-4.0 if _r.random() < 0.08 else 0.0) for _ in range(600)]
_stressed = {k: [v + s for v, s in zip(vals, _shock)] for k, vals in _base.items()}
_, _, lam_b, _ = P.tail_dependence(_base, 0.10)
_, _, lam_st, _ = P.tail_dependence(_stressed, 0.10)
ok("tail dependence detects a common shock", lam_st > lam_b + 0.2, f"{lam_b:.3f} -> {lam_st:.3f}")

# ══════════════════════════════════ 16d. monte carlo (45-49) vs known nulls
import qmontecarlo as MC                              # noqa: E402
_rr = _rnd.Random(5)
_noise = [_rr.gauss(0, 0.01) for _ in range(200)]
_edge = [_rr.gauss(0.004, 0.01) for _ in range(200)]
sf_n = MC.sign_flip_test(_noise, "mean", 1500, 1)
sf_e = MC.sign_flip_test(_edge, "mean", 1500, 1)
ok("sign-flip: noise not significant", sf_n["p_value"] > 0.10, f"p={sf_n['p_value']:.3f}")
ok("sign-flip: real edge significant", sf_e["p_value"] < 0.01, f"p={sf_e['p_value']:.3f}")
ok("sign-flip p is never exactly 0", sf_e["p_value"] > 0)
ok("permutation refuses mean/total", math.isnan(MC.permutation_test(_edge, "mean")["p_value"]))
pm_iid = MC.permutation_test(_noise, "sharpe", 1000, 1)
ok("permutation on iid series: order does not matter", pm_iid["p_value"] > 0.5,
   f"p={pm_iid['p_value']:.3f}")
_timing = [0.01 if i < 100 else -0.01 for i in range(200)]
pm_t = MC.permutation_test(_timing, "max_dd", 1000, 1)
ok("permutation detects timing structure", pm_t["p_value"] < 0.01, f"p={pm_t['p_value']:.3f}")
for _m in ("iid_bootstrap", "block_bootstrap", "gbm", "ou"):
    sg = MC.surrogate_paths(_edge, _m, "sharpe", 800, 1)
    ok(f"surrogate({_m}) rejects a real edge", sg["p_value"] < 0.05, f"p={sg['p_value']:.3f}")
    sgn = MC.surrogate_paths(_noise, _m, "sharpe", 800, 1)
    ok(f"surrogate({_m}) accepts noise", sgn["p_value"] > 0.10, f"p={sgn['p_value']:.3f}")
ok("surrogate reports its calibration",
   "observed_sd" in MC.surrogate_paths(_edge, "gbm", "sharpe", 200, 1)["calibration"])
_rc_noise = [[_rr.gauss(0, 0.01) for _ in range(150)] for _ in range(8)]
_rc_good = [[_rr.gauss(0.006, 0.01) for _ in range(150)]] + \
           [[_rr.gauss(0, 0.01) for _ in range(150)] for _ in range(7)]
rc1, rc2 = MC.reality_check(_rc_noise, 500, 1), MC.reality_check(_rc_good, 500, 1)
ok("reality check rejects a screen of pure noise", rc1["p_value"] > 0.10, f"p={rc1['p_value']:.3f}")
ok("reality check finds a genuinely superior strategy", rc2["p_value"] < 0.05, f"p={rc2['p_value']:.3f}")
ok("reality check identifies WHICH strategy", rc2["best_strategy_index"] == 0)
sp1, sp2 = MC.spa_test(_rc_noise, 500, 1), MC.spa_test(_rc_good, 500, 1)
ok("SPA rejects noise", sp1["p_value"] > 0.10, f"p={sp1['p_value']:.3f}")
ok("SPA finds the real one", sp2["p_value"] < 0.05, f"p={sp2['p_value']:.3f}")
ok("SPA drops hopeless candidates", sp2["n_retained"] <= sp2["n_strategies"])
ok("reality check guards empty input", math.isnan(MC.reality_check([])["p_value"]))
sm = MC.stress_moments(_edge, cost_per_period=0.001)
ok("stress has all scenarios", len(sm["rows"]) == 7, str(len(sm["rows"])))
ok("stress: strong edge survives", sm["n_scenarios_surviving_t1.6"] >= 5,
   f"{sm['n_scenarios_surviving_t1.6']}/6")
sm2 = MC.stress_moments(_noise, cost_per_period=0.001)
ok("stress: noise survives nothing", sm2["n_scenarios_surviving_t1.6"] == 0,
   f"{sm2['n_scenarios_surviving_t1.6']}/6")
ok("monte carlo report runs", MC.report(_edge, 0.001, reps=200).startswith("==="))

# ══════════════════════════════════ 16e. data loader (offline paths only)
import qdata as DL                                   # noqa: E402
_csv = os.path.join(_tf.mkdtemp(), "bars.csv")
with open(_csv, "w") as _f:
    _f.write("time,open,high,low,close,volume\n")
    for _i in range(10):
        _f.write(f"{1700000000 + _i*300},{100+_i},{101+_i},{99+_i},{100.5+_i},{1000+_i}\n")
_rows = DL.from_csv(_csv, v_col=5)
ok("from_csv reads all rows", len(_rows) == 10, str(len(_rows)))
ok("from_csv parses OHLC", _rows[0]["o"] == 100.0 and _rows[9]["c"] == 109.5)
ok("from_csv sorts by time", all(_rows[i]["t"] <= _rows[i+1]["t"] for i in range(9)))
_ret = DL.to_returns(_rows)
ok("to_returns length is n-1", len(_ret) == 9)
close("to_returns first value", _ret[0], 101.5/100.5 - 1, 1e-12)
_st = os.path.join(_tf.mkdtemp(), "store.csv")
_a = DL.append_csv(_st, _rows[:5])
_b = DL.append_csv(_st, _rows)                        # overlapping -> must dedupe
ok("append_csv is idempotent/dedupes", _a["total_bars"] == 5 and _b["total_bars"] == 10,
   f"{_a['total_bars']} then {_b['total_bars']}")
ok("append_csv reports the window", _b["first"] is not None and _b["last"] is not None)
ok("RETENTION documents every interval used",
   all(k in DL.RETENTION for k in ("1m", "5m", "1h", "1d")))
ok("qdata is the only networked module",
   all("urllib" not in open(os.path.join(HERE, "..", "scripts", m)).read()
       for m in ("qstats.py", "qvalidate.py", "qsize.py", "qportfolio.py", "qmontecarlo.py")))

# ══════════════════════════════════ 17. degenerate inputs must not crash
ok("describe n=1 no crash", S.describe([1.0])["n"] == 1)
ok("describe zero variance", S.describe([3.0, 3.0, 3.0])["sd"] == 0.0)
ok("var_cvar single value", S.var_cvar([5.0])["worst"] == 5.0)
ok("drawdown single value", S.drawdown([1.0])["max_dd"] == 0.0)
ok("effective_n tiny input", S.effective_n([1.0, 2.0]) == 2)
ok("newey_west tiny input", math.isnan(V.newey_west_se([1.0])[0]))
ok("kelly zero sd", Z.kelly_continuous(0.1, 0.0) == float('inf'))
ok("impact zero volume", Z.impact_sqrt(100, 0, 0.02) == float('inf'))
ok("report runs", S.report([1.0, -2.0, 3.0], "x", 5).startswith("==="))
ok("validate report runs", V.report([0.01, -0.02, 0.03], 5, 0.2).startswith("==="))
ok("size report runs", Z.report(0.01, 0.1).startswith("==="))

# ══════════════════════════════════ 18. qmicro + RAS (v1.4)
import random as _rnd
import qmicro as M          # noqa: E402

# --- Corwin-Schultz e Roll contro uno spread NOTO ---------------------------------
_r = _rnd.Random(11)
_TRUE = 0.02
_mid, _hi, _lo, _cl = 100.0, [], [], []
for _ in range(3000):
    _mid *= math.exp(_r.gauss(0, 0.0008))
    _rng = abs(_r.gauss(0, 0.0008)) * _mid
    _hi.append(_mid + _rng + _TRUE / 2)
    _lo.append(_mid - _rng - _TRUE / 2)
    _cl.append(_mid + (_TRUE / 2 if _r.random() < 0.5 else -_TRUE / 2))
_cs = M.corwin_schultz_median(_hi, _lo, window=1)
ok("corwin_schultz produce una stima", _cs is not None)
ok("corwin_schultz in range plausibile", _cs is not None and 0.0 <= _cs <= 0.01, f"got {_cs}")
_rs = M.roll_spread(_cl)
ok("roll_spread positivo con bounce", _rs is not None and _rs > 0)
ok("roll_spread nell ordine dello spread vero",
   _rs is not None and abs(_rs - _TRUE) < _TRUE, f"got {_rs:.5f} vs {_TRUE}")
ok("roll_spread ritorna None su trend lineare (cov solo rumore numerico)",
   M.roll_spread([100.0 + 0.1 * i for i in range(200)]) is None)
_ar, _p, _dp = [], 100.0, 0.0
for _ in range(500):                      # AR(1) con phi>0: incrementi persistenti
    _dp = 0.7 * _dp + _r.gauss(0, 0.05)
    _p += _dp; _ar.append(_p)
ok("roll_spread ritorna None su incrementi persistenti (phi>0)", M.roll_spread(_ar) is None)
ok("roll_spread None su input minuscolo", M.roll_spread([1.0, 2.0]) is None)
ok("corwin_schultz lista vuota su input minuscolo", M.corwin_schultz_spread([1.0], [1.0]) == [])

# --- cost cliff: identita' verificabili -------------------------------------------
_cc = M.cost_cliff(2.0, 1000, 1.0, 10.0)
ok("cost_cliff netto < lordo", _cc["net_sharpe"] < _cc["gross_sharpe"])
close("cost_cliff costo zero => netto == lordo", M.cost_cliff(2.0, 1000, 0.0, 10.0)["net_sharpe"], 2.0, 1e-9)
_f = M.max_trade_frequency(2.0, 1.0, 1.0, 10.0)
close("max_trade_frequency inverte cost_cliff", M.cost_cliff(2.0, _f, 1.0, 10.0)["net_sharpe"], 1.0, 1e-6)
close("required_edge_per_trade al 20%", M.required_edge_per_trade(1.0, 0.20), 5.0, 1e-9)
ok("cost_cliff segnala il non-sopravvissuto", not M.cost_cliff(0.5, 5000, 2.0, 10.0)["survives"])

# --- decomposizione slippage: il residuo DEVE chiudere a zero ---------------------
_d = M.decompose_slippage(100.0, 100.5, 110.0, 109.0, side=+1,
                          commission_per_round_trip=2.0, multiplier=1.0)
close("slippage pnl teorico", _d["theoretical_pnl"], 10.0, 1e-9)
close("slippage pnl reale", _d["realised_pnl"], 6.5, 1e-9)
close("slippage entrata", _d["entry_slippage"], -0.5, 1e-9)
close("slippage uscita", _d["exit_slippage"], -1.0, 1e-9)
close("slippage residuo chiude a zero", _d["unexplained_residual"], 0.0, 1e-9)
_ds = M.decompose_slippage(110.0, 109.5, 100.0, 101.0, side=-1, commission_per_round_trip=2.0)
close("slippage short: residuo chiude a zero", _ds["unexplained_residual"], 0.0, 1e-9)

# --- Rademacher / RAS -------------------------------------------------------------
_r2 = _rnd.Random(3)
_T = 600
_base = [_r2.gauss(0, 0.01) for _ in range(_T)]
_corr = [[b + _r2.gauss(0, 0.001) for b in _base] for _ in range(30)]
_indep = [[_r2.gauss(0, 0.01) for _ in range(_T)] for _ in range(30)]
_rc = V.rademacher_complexity(_corr, reps=400, seed=2)
_ri = V.rademacher_complexity(_indep, reps=400, seed=2)
ok("Rademacher: correlate < indipendenti", _rc < _ri, f"{_rc:.5f} < {_ri:.5f}")
_mb = V.massart_bound(30, _T)
ok("Rademacher empirico <= Massart (correlate)", _rc <= _mb * 1.05)
ok("Rademacher empirico <= Massart (indipendenti)", _ri <= _mb * 1.05)
ok("Massart cresce con N", V.massart_bound(100, _T) > V.massart_bound(10, _T))
ok("trial efficaci: correlate < indipendenti",
   V.effective_trials(30, _T, _rc) < V.effective_trials(30, _T, _ri))
_hi_ras = V.ras_sharpe_adjustment(5.0, _corr, reps=400, seed=2)
_lo_ras = V.ras_sharpe_adjustment(0.05, _corr, reps=400, seed=2)
ok("RAS: Sharpe alto sopravvive", _hi_ras["sharpe_floor_ann"] > 0)
ok("RAS: Sharpe marginale non sopravvive", _lo_ras["sharpe_floor_ann"] < 0)
ok("RAS: Massart >= empirico", _hi_ras["massart_ann"] >= _hi_ras["rademacher_ann"])
ok("RAS gestisce input minuscolo", "error" in V.ras_sharpe_adjustment(1.0, [[0.1, 0.2]]))
ok("rademacher su lista vuota non crasha", math.isnan(V.rademacher_complexity([])))
ok("qmicro non tocca la rete", "urllib" not in open(os.path.join(HERE, "..", "scripts", "qmicro.py")).read())

# --- guardia di frequenza: lo stimatore va applicato alla risoluzione giusta -------
# Discriminante corretto = TICK impliciti, non quota del range di barra (il primo
# tentativo usava il range e non discriminava: ~50% su daily, ~29% a 1 minuto).
_pl_ok = M.plausibility_flag(1.57e-4, mid_price=20000.0, tick_size=0.25)      # ~1.3 tick
_pl_no = M.plausibility_flag(62.5e-4, mid_price=4200.0, tick_size=0.25)        # ~105 tick
ok("plausibility: stima a 1 tick e' credibile", _pl_ok["plausible"] is True,
   f"{_pl_ok['implied_ticks']:.1f} tick")
ok("plausibility: stima a 100+ tick NON e' credibile", _pl_no["plausible"] is False,
   f"{_pl_no['implied_ticks']:.0f} tick")
ok("plausibility: input mancante -> None, non un falso ok",
   M.plausibility_flag(None, 100.0, 0.25)["plausible"] is None)
ok("bounce_check espone la plausibilita' di entrambi gli stimatori",
   all(k in M.bounce_check(_cl, _hi, _lo, 0.02)
       for k in ("roll_plausibility", "corwin_schultz_plausibility")))

# ══════════════════════════════════ 19. qpairs (v1.4)
import qpairs as PR        # noqa: E402

_rp = _rnd.Random(5)
# (a) coppia COINTEGRATA sintetica: y = 2x + OU(phi=0.92)
_x, _sp, _X, _Y = 100.0, 0.0, [], []
for _ in range(1500):
    _x += _rp.gauss(0, 0.5)
    _sp = 0.92 * _sp + _rp.gauss(0, 0.4)
    _X.append(_x); _Y.append(2.0 * _x + _sp)
_eg = PR.engle_granger(_Y, _X)
close("engle_granger recupera l hedge ratio vero", _eg["hedge_ratio"], 2.0, 0.15)
ok("engle_granger: coppia cointegrata riconosciuta", _eg["cointegrated_at_5pct"])
_HL_TRUE = -math.log(2) / math.log(0.92)
close("half-life vicina alla teorica", _eg["half_life"], _HL_TRUE, 0.5 * _HL_TRUE)

# (b) CONTROLLO NEGATIVO — e' un test al 5%: una SINGOLA estrazione fallisce il 5% delle
#     volte per costruzione. Si verifica quindi la CALIBRAZIONE su molte repliche, che e'
#     l'unica proprieta' davvero interessante di una procedura statistica.
def _mc_rate(make, predicate, reps=120, seed0=4000):
    hits = 0
    for _s in range(reps):
        hits += 1 if predicate(make(_rnd.Random(seed0 + _s))) else 0
    return hits / reps


def _two_walks(rr, n=800):
    a = b = 100.0; A = []; B = []
    for _ in range(n):
        a += rr.gauss(0, 0.5); b += rr.gauss(0, 0.5); A.append(a); B.append(b)
    return A, B


_A, _B = _two_walks(_rnd.Random(77), n=1500)   # coppia di riferimento per i confronti
_fp_eg = _mc_rate(_two_walks, lambda ab: PR.engle_granger(ab[0], ab[1])["cointegrated_at_5pct"])
ok("Engle-Granger: falsi positivi su random walk vicini al 5%", 0.0 <= _fp_eg <= 0.15,
   f"tasso {100*_fp_eg:.1f}% su 120 repliche")
ok("Engle-Granger usa soglie piu severe dell ADF semplice",
   PR.engle_granger(_A, _B)["critical_values"]["5%"] < PR.adf_stat(_A)["critical_values"]["5%"])


def _one_walk(rr, n=800):
    w = [100.0]
    for _ in range(n):
        w.append(w[-1] + rr.gauss(0, 0.5))
    return w


_fp_adf = _mc_rate(_one_walk, lambda w: PR.adf_stat(w)["stationary_at_5pct"], seed0=6000)
ok("ADF: falsi positivi su random walk vicini al 5% (calibrazione)", 0.0 <= _fp_adf <= 0.15,
   f"tasso {100*_fp_adf:.1f}% su 120 repliche")

# (c) Hurst: mediana su repliche, non singola estrazione
def _median(v):
    v = sorted(v); return v[len(v) // 2]


_h_rw, _h_ou, _h_tr = [], [], []
for _s in range(40):
    _r2 = _rnd.Random(8000 + _s)
    _h_rw.append(PR.hurst(_one_walk(_r2))["hurst"])
    _z = [0.0]
    for _ in range(1200):
        _z.append(0.9 * _z[-1] + _r2.gauss(0, 1))
    _h_ou.append(PR.hurst(_z)["hurst"])
    _h_tr.append(PR.hurst([float(_i) * 0.5 + _r2.gauss(0, 3.0) for _i in range(800)])["hurst"])
close("Hurst mediano su random walk ~ 0.5", _median(_h_rw), 0.5, 0.08)
ok("Hurst mediano su OU sotto 0.45", _median(_h_ou) < 0.45, f"{_median(_h_ou):.3f}")
ok("Hurst mediano su trend sopra 0.55", _median(_h_tr) > 0.55, f"{_median(_h_tr):.3f}")

# (d) half-life: la MEDIANA su repliche deve centrare la formula chiusa
_HL_AR = -math.log(2) / math.log(0.9)
_hls = []
for _s in range(60):
    _r3 = _rnd.Random(9000 + _s); _z = [0.0]
    for _ in range(2000):
        _z.append(0.9 * _z[-1] + _r3.gauss(0, 1))
    _hls.append(PR.ou_half_life(_z)["half_life"])
close("half-life: mediana su 60 repliche vs formula AR(1)", _median(_hls), _HL_AR, 0.20 * _HL_AR)
ok("half-life: nessuna replica assurda", all(1.0 < v < 30.0 for v in _hls))
ok("half-life inf/lunga se non c'e reversione",
   PR.ou_half_life(_A)["half_life"] > 50 or math.isinf(PR.ou_half_life(_A)["half_life"]))
ok("ou_half_life gestisce campione minuscolo",
   math.isnan(PR.ou_half_life([1.0, 2.0, 3.0])["half_life"]))
_ou = [0.0]                                   # OU di riferimento (phi=0.9)
for _ in range(2000):
    _ou.append(0.9 * _ou[-1] + _rnd.Random(4242).gauss(0, 1) if False else 0.9 * _ou[-1])
_ourng = _rnd.Random(4242)
_ou = [0.0]
for _ in range(2000):
    _ou.append(0.9 * _ou[-1] + _ourng.gauss(0, 1))
ok("ADF non inventa un p-value", "p_value" not in PR.adf_stat(_ou))
ok("ADF: OU stazionario al 5%", PR.adf_stat(_ou)["stationary_at_5pct"])

# (f) Kalman: su beta costante deve convergere al beta vero
_kx = [float(i) + _rp.gauss(0, 1) for i in range(400)]
_ky = [3.0 * v + _rp.gauss(0, 0.5) for v in _kx]
_k = PR.kalman_hedge_ratio(_ky, _kx)
close("kalman converge al beta vero", _k["beta_final"], 3.0, 0.4)
ok("kalman restituisce il percorso di beta", len(_k["beta"]) == len(_kx))
ok("qpairs non tocca la rete",
   "urllib" not in open(os.path.join(HERE, "..", "scripts", "qpairs.py")).read())

# ============================================================ v1.5 — le guardie del silenzio
# Ogni assert qui riproduce un errore REALMENTE commesso il 2026-07-30. Se uno fallisce, la
# guardia corrispondente non prende piu' il bug che esiste per prendere.
print("\n-- v1.5: guardie contro i fallimenti silenziosi --")

# soglia irraggiungibile: z su finestra 4 con sd di popolazione ha |z|max = sqrt(3) = 1.732
_r = S.reachable_threshold(4, 2.0, ddof=0)
ok("reachable: W=4 ddof=0 soglia 2 e' irraggiungibile", _r["reachable"] is False)
close("reachable: |z|max = sqrt(W-1)", _r["max_abs"], math.sqrt(3), 1e-12)
ok("reachable: propone una finestra minima", _r["min_window_needed"] >= 5)
close("reachable: ddof=1 da (W-1)/sqrt(W)", S.reachable_threshold(4, 2.0, ddof=1)["max_abs"],
      1.5, 1e-12)
ok("reachable: W=240 soglia 2 va bene", S.reachable_threshold(240, 2.0, 0)["reachable"] is True)
ok("reachable: soglia bassa non segnalata", S.reachable_threshold(100, 1.0, 0)["reachable"] is True)

# drawdown in scala log spacciato per percentuale: -0.652 log = -47.9% reale, non -65.2%
close("log_to_simple: -0.652 log = -47.9% reale", 100 * S.log_to_simple(-0.652), -47.9, 0.1)
close("log_to_simple: 0 -> 0", S.log_to_simple(0.0), 0.0, 1e-15)

# tasso di evento: ORB che rientrava 7.5x per un flag mancante
_e = S.event_rate_check(54920, 10.0, 730)
ok("event_rate: prende 7.5x troppi eventi", _e["ok"] is False and _e["ratio"] > 7)
ok("event_rate: zero eventi non e' un negativo", "ZERO EVENTI" in S.event_rate_check(0, 10.0, 730)["note"])
ok("event_rate: tasso corretto passa", S.event_rate_check(7300, 10.0, 730)["ok"] is True)

# due guardie di plausibilita' su casi DISGIUNTI
_a = MI.plausibility_flag(62e-4, 2128.0, 0.25)                       # ES daily
close("guardia A: ES daily implica 53 tick", _a["implied_ticks"], 52.77, 0.1)
ok("guardia A: applicabile su ES (tick grossolano)", _a["guard_a_applicable"] is True)
ok("guardia A: respinge ES daily", _a["plausible"] is False)
ok("guardia A sola -> coverage 'A only'", _a["coverage"] == "A only")
_ba = MI.plausibility_flag(17.87e-4, 58000.0, 0.01)                  # BTC, tick inutile
ok("BTC: il conteggio in tick e' enorme ma privo di senso", _ba["implied_ticks"] > 4000)
ok("BTC: la guardia A si dichiara NON APPLICABILE", _ba["guard_a_applicable"] is False)
ok("BTC: senza la B non si emette giudizio", _ba["plausible"] is None)
ok("BTC: e lo dice", "NON VERIFICABILE" in _ba["note"])
_b = MI.plausibility_flag(17.87e-4, 58000.0, 0.01, median_bar_range_rel=6.15e-4)
ok("guardia B: su BTC copre da sola (A non applicabile)", _b["coverage"] == "B only")
ok("guardia B: respinge BTC sulla fisica", _b["plausible"] is False and _b["guard_b_ok"] is False)
close("guardia B: 17.87bp = 2.9x il range di 6.15bp", _b["range_share"], 17.87 / 6.15, 1e-9)
_g = MI.plausibility_flag(1.0e-4, 58000.0, 0.01, median_bar_range_rel=6.15e-4)
ok("guardia: 1bp su BTC (lo spread VERO) viene ACCETTATA", _g["plausible"] is True)
ok("guardia: e passa per la B, non per la A", _g["coverage"] == "B only")

# il numero che decide il verdetto: soglia di ribaltamento
_f = lambda pen: 0.872 * 6.0 - (1 - 0.872) * pen
close("sign_flip: il penalty del box si ribalta a 40.9c",
      MI.sign_flip_parameter(_f, 0.0, 100.0), 40.875, 0.05)
ok("sign_flip: None se non cambia segno nell'intervallo",
   MI.sign_flip_parameter(_f, 0.0, 10.0) is None)
_bp = MI.borrowed_parameter("penalty", 9.0, "README terzo", _f, 0.0, 50.0)
ok("in prestito non misurato -> NON consegnare il verdetto", "NON consegnare" in _bp["deliverable"])
_bm = MI.borrowed_parameter("penalty", 9.0, "README terzo", _f, 0.0, 50.0, measured=44.3)
ok("in prestito misurato -> verdetto negativo", _bm["verdict_at_used"] == "negative")
ok("in prestito: segnala la differenza di categoria", "CATEGORIA" in _bm.get("note", ""))

# breadth: rifiutare la formula con rho <= 0
_rng = random.Random(11)
_base = [_rng.gauss(0, 1) for _ in range(600)]
_pos = [[0.7 * _base[i] + 0.7 * _rng.gauss(0, 1) for i in range(600)] for _ in range(6)]
_rp = PF.effective_breadth(_pos)
ok("breadth: valida con rho > 0", _rp["valid_as_breadth"] is True and _rp["breadth"] is not None)
ok("breadth: sta fra 1 e N", 1.0 < _rp["breadth"] < 6.0)
_x = [_rng.gauss(0, 1) for _ in range(600)]
_ra = PF.effective_breadth([_x, [-v for v in _x]])
ok("breadth: due serie che si annullano hanno rho ~ -1", _ra["mean_corr"] < -0.9)
ok("breadth: RIFIUTA di chiamarla breadth con rho <= 0",
   _ra["valid_as_breadth"] is False and _ra["breadth"] is None)
ok("breadth: lo dichiara nella nota", "NON va letta come breadth" in _ra["note"])

# diversificare con mu < 0 peggiora lo Sharpe
ok("segno: mu > 0 e' sicuro", PF.diversification_sign_check(0.001)["safe_to_diversify"] is True)
_bd = PF.diversification_sign_check(-0.0004)
ok("segno: mu < 0 avverte", _bd["safe_to_diversify"] is False and "PEGGIORA" in _bd["note"])


# ============================================================ v1.6 — le affermazioni generate
# Ogni assert riproduce una frase REALMENTE consegnata e sbagliata (lug-ago 2026). Se uno
# fallisce, la funzione non impedisce piu' l'errore che esiste per impedire.
print("\n-- v1.6: la prosa generata dai numeri --")

# "11.71 resta sotto 8.37" — scritto due righe sotto una tabella che mostrava il contrario
_c = QC.cmp_(11.71, 8.37, "breadth 7 famiglie", "breadth 34 strumenti")
ok("cmp_: dice SOPRA quando a > b", "SOPRA" in _c)
ok("cmp_: non dice SOTTO", "SOTTO" not in _c)
ok("cmp_: simmetrico", "SOTTO" in QC.cmp_(8.37, 11.71, "x", "y"))
ok("cmp_: uguaglianza riconosciuta", "UGUALE" in QC.cmp_(2.0, 2.0, "x", "y"))
ok("cmp_: riporta entrambi i valori", "11.71" in _c and "8.37" in _c)

# "3.50 e' la media del rumore" quando la media era 3.013
ok("rel: soglia falsa e' marcata", "NON regge" in QC.rel(3.503, "<=", 3.013, "sqrt(2 ln N)"))
ok("rel: soglia vera passa", "VERO" in QC.rel(2.02, ">=", 2.0, "t"))
ok("rel: nessun avviso se vera", "NON regge" not in QC.rel(2.02, ">=", 2.0, "t"))
try:
    QC.rel(1.0, ">", 2.0, "x", raise_on_fail=True); _raised = False
except QC.ClaimError:
    _raised = True
ok("rel: solleva ClaimError se richiesto", _raised)

# "aggiornata la casella del gate" quando grep -c diceva 0
ok("rel: postcondizione a zero e' falsa", "NON regge" in QC.rel(0, ">", 0, "occorrenze grep"))

# "salva il 60-70%" quando la misura diceva 25%
_a = QC.audit("il coordinatore salva il 60-70% dei costi", {"risparmio": 25.0})
ok("audit: prende i numeri inventati", _a["ok"] is False)
ok("audit: li elenca", 60.0 in _a["untraceable"])
_b = QC.audit("il t vale 2.02 contro una soglia di 1.6", {"t": 2.02, "soglia": 1.6})
ok("audit: passa se ogni numero e' tracciabile", _b["ok"] is True)
ok("audit: tollera percentuali (0.25 vs 25)", QC.audit("il 25%", {"q": 0.25})["ok"] is True)

# variazione col segno giusto
ok("pct_change: un calo si chiama calo", "cala" in QC.pct_change(25, 68, "risparmio"))
ok("pct_change: una crescita si chiama crescita", "aumenta" in QC.pct_change(68, 25, "x"))

# il verdetto si deriva dai criteri dichiarati, non si sceglie
_v = QC.verdict({"t": 0.68, "floor": -1.02},
                [("t sotto 1.6", lambda m: m["t"] < 1.6, "DECIDIBILE-NEGATIVO"),
                 ("floor negativo", lambda m: m["floor"] < 0, "NON ESTABLISHED")])
ok("verdict: applica il primo criterio che scatta", _v["verdict"] == "DECIDIBILE-NEGATIVO")
ok("verdict: dichiara quale", _v["matched"].startswith("t sotto"))
_v2 = QC.verdict({"t": 3.0}, [("t sotto 1.6", lambda m: m["t"] < 1.6, "NEG")])
ok("verdict: default se nessun criterio scatta", _v2["verdict"] == "NON DECIDIBILE")
ok("verdict: un criterio rotto non passa in silenzio",
   QC.verdict({}, [("chiave assente", lambda m: m["manca"] < 1, "X")])["verdict"] == "NON DECIDIBILE")


# ══════════════════════════════════ ndarray attraversa il confine lista (R13)
# Errore REALE del 2026-08-10: bootstrap(block=5) su un ndarray. Nel ramo di avvolgimento
# `x[i:] + x[:k]` il `+` deve CONCATENARE; su un ndarray e' una SOMMA elementwise. Quando le
# due fette hanno lunghezza diversa esplode — ed e' il caso fortunato. Quando hanno la stessa
# lunghezza (o una delle due e' lunga 1, e allora numpy trasmette) somma i numeri e restituisce
# una banda plausibile e SBAGLIATA, senza alcun errore. Il tipo dichiarato non cambia mai.
# This section and the qconfluence and qregime ones need numpy; everything above needs only
# the stdlib. Without numpy they are SKIPs, not a crash that hides every result above them.
try:
    import numpy as _np
except ImportError as e:
    _np, _no_np = None, str(e)
if _np is None:
    SKIP.append(f"ndarray attraversa il confine lista (R13): {_no_np}")
else:
    _arr = _np.array([1.0, -2.0, 3.0, -1.0, 2.0, 0.5, -0.5, 1.5, -1.5, 2.5])
    _lst = [float(v) for v in _arr]
    _ba = S.bootstrap(_arr, stat=_np.mean, reps=200, seed=3, block=4)
    _bl = S.bootstrap(_lst, stat=_np.mean, reps=200, seed=3, block=4)
    ok("bootstrap: ndarray con blocchi non esplode", isinstance(_ba["p05"], float) or _ba["p05"] == _ba["p05"])
    ok("bootstrap: ndarray e lista danno lo STESSO risultato", abs(_ba["p50"] - _bl["p50"]) < 1e-12)
    ok("bootstrap: la banda resta dentro l'intervallo dei dati",
       min(_lst) <= _ba["p05"] <= max(_lst) and min(_lst) <= _ba["p95"] <= max(_lst))
    _ja = S.jackknife(_arr, stat=_np.mean)
    _jl = S.jackknife(_lst, stat=_np.mean)
    ok("jackknife: ndarray non esplode e coincide con la lista",
       abs(_ja["full"] - _jl["full"]) < 1e-12 and _ja["sign_survives"] == _jl["sign_survives"])
    import qmontecarlo as QMC
    _ma = QMC.surrogate_paths(_arr, model="block_bootstrap", block=4, reps=100, seed=3, drift_zero=False)
    ok("surrogate_paths: block_bootstrap con drift_zero=False su ndarray non esplode",
       isinstance(_ma, dict))


# ══════════════════════════════════ qconfluence: the five gates (§5b)
# References are closed forms computed here or effects planted in the data, never the
# module's own output.
if _np is None:
    SKIP.append(f"qconfluence: {_no_np}")
else:
    import statistics as _st
    import qconfluence as QF

    close("sizing_cost: 1x/3x on a coin flip = 2/sqrt(5)",
          QF.sizing_cost([1.0, 3.0] * 50), 2 / math.sqrt(5), 1e-12)
    close("sizing_cost: 3x a quarter of the time = 1.5/sqrt(3)",
          QF.sizing_cost([3.0, 1.0, 1.0, 1.0] * 25), 1.5 / math.sqrt(3), 1e-12)
    close("sizing_cost: equal weights cost nothing", QF.sizing_cost([2.0] * 10), 1.0, 1e-12)
    # The formula against what it claims to be: Sharpe of weighted over unweighted returns,
    # with weights orthogonal to the returns by construction (period 4).
    _w = _np.tile([1.0, 3.0, 1.0, 3.0], 250)
    _r = 0.01 + _np.tile([1.0, 1.0, -1.0, -1.0], 250)
    close("sizing_cost: equals Sharpe(w·r)/Sharpe(r) when the weights are noise",
          QF.sizing_cost(_w), ((_w * _r).mean() / (_w * _r).std()) / (_r.mean() / _r.std()), 1e-4)

    # z of [0, 1, 2] with the sample SD is [-1, 0, 1], so the weights are 1 + slope·z.
    _sw = QF.sizing_from_slope([0.0, 1.0, 2.0], 0.5, 1.0)
    ok("sizing_from_slope: hand case [0.5, 1, 1.5]", _np.allclose(_sw, [0.5, 1.0, 1.5], atol=1e-12), str(_sw))
    _sw = QF.sizing_from_slope([0.0, 1.0, 2.0], 10.0, 1.0)
    ok("sizing_from_slope: capped at [0.5, 2]", _np.allclose(_sw, [0.5, 1.0, 2.0], atol=1e-12), str(_sw))
    ok("sizing_from_slope: zero base mean gives flat weights",
       _np.allclose(QF.sizing_from_slope([0.0, 1.0, 2.0], 0.5, 0.0), 1.0))

    _sc = [float(v) for v in range(200)]
    _dr = QF.dose_response(_sc, [2.0 * v for v in _sc])
    close("dose_response: slope per SD of a noiseless line = 2·stdev(score)",
          _dr.get("slope_per_sd", float("nan")), 2.0 * _st.stdev(_sc), 1e-9)
    ok("dose_response: every observation lands in exactly one bucket", sum(_dr.get("n", [])) == 200, str(_dr.get("n", [])))
    ok("dose_response: a noiseless rising line is monotone up",
       _dr.get("monotone_up") and not _dr.get("monotone_down"))
    close("dose_response: Spearman of a monotone line is 1", _dr.get("spearman", float("nan")), 1.0, 1e-12)
    # The trap the module documents: a discrete 0-3 count with a REAL planted effect. Skewed
    # the way confluence counts are (most periods align few conditions), quintile edges come
    # out 0,0,0,1,1 and merge 1-3 into one bucket; one bucket per value keeps all four.
    _g = _np.random.default_rng(11)
    _sc = _g.choice([0.0, 1.0, 2.0, 3.0], size=4000, p=[0.5, 0.3, 0.15, 0.05])
    _pn = 0.25 * _sc + _g.normal(0, 1, 4000)
    _dr = QF.dose_response(_sc, _pn)
    ok("dose_response: a discrete 0-3 count gets one bucket per value",
       _dr.get("ok") and len(_dr.get("buckets", [])) == 4, str(_dr.get("n")))
    ok("dose_response: bucket means within 4 SE of the planted 0.25·k",
       _dr.get("ok") and all(abs(m - 0.25 * k) <= 4 * s for k, (m, s) in enumerate(zip(_dr.get("buckets", []), _dr.get("se", [])))),
       str(_dr.get("buckets")))
    ok("dose_response: the planted effect reads as monotone up", _dr.get("ok") and _dr.get("monotone_up"))
    ok("dose_response: a value with too few periods is dropped, not kept as nan",
       QF.dose_response([0.0] * 100 + [1.0] * 100 + [2.0] * 100 + [3.0] * 5,
                        [float(i % 7) for i in range(305)]).get("n") == [100, 100, 100])

    def _avg_ranks(v):
        order = sorted(range(len(v)), key=v.__getitem__)
        ranks, i = [0.0] * len(v), 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in order[i:j + 1]:
                ranks[k] = (i + j) / 2
            i = j + 1
        return ranks

    def _pearson(a, b):
        ma, mb = sum(a) / len(a), sum(b) / len(b)
        return (sum((x - ma) * (y - mb) for x, y in zip(a, b))
                / math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)))

    # Swapping each pair of 0..99 makes every |d| = 1, so rho = 1 - 6n/(n(n²-1)).
    _swap = [float(v + 1 if v % 2 == 0 else v - 1) for v in range(100)]
    close("dose_response: Spearman without ties = 1 - 6/(n²-1)",
          QF.dose_response([float(v) for v in range(100)], _swap).get("spearman", float("nan")), 1 - 6 / (100 ** 2 - 1), 1e-12)
    _sl, _pl = [float(v) for v in _sc[:400]], [float(v) for v in _pn[:400]]
    close("dose_response: Spearman with tied scores uses average ranks",
          QF.dose_response(_sc[:400], _pn[:400]).get("spearman", float("nan")), _pearson(_avg_ranks(_sl), _avg_ranks(_pl)), 1e-9)

    # 20 disjoint rules on centred noise: the best t is the max of 20 t(49), whose quantiles
    # are the t quantiles of 0.5^(1/20) and 0.95^(1/20) (t from z: Abramowitz-Stegun 26.7.5).
    _g = _np.random.default_rng(0)
    _pn0 = _g.normal(0, 1, 1000)
    _pn0 -= _pn0.mean()
    _mk = _np.zeros((20, 1000), bool)
    for _i in range(20):
        _mk[_i, 50 * _i:50 * (_i + 1)] = True
    _sn = QF.search_null(_mk, _pn0, n_rep=400, seed=0)

    def _t_quantile(p, df):
        z = _st.NormalDist().inv_cdf(p)
        return z + (z ** 3 + z) / (4 * df)

    close("search_null: median best-of-20 t = max of 20 t(49)", _sn.get("median", float("nan")), _t_quantile(0.5 ** (1 / 20), 49), 0.1)
    close("search_null: 95th percentile of the best-of-20 t", _sn.get("q95", float("nan")), _t_quantile(0.95 ** (1 / 20), 49), 0.3)
    ok("search_null: no rule with enough periods is not a verdict",
       QF.search_null(_mk, _pn0, min_n=51).get("ok") is False)

    _jk = QF.jackknife_sign([float(v) for v in range(1, 21)])
    close("jackknife_sign: worst leave-one-out mean of 1..20 = 190/19", _jk.get("worst_mean", float("nan")), 10.0, 1e-12)
    close("jackknife_sign: best leave-one-out mean of 1..20 = 209/19", _jk.get("best_mean", float("nan")), 11.0, 1e-12)
    close("jackknife_sign: top-3 share of 1..20 = 57/210", _jk.get("top3_share", float("nan")), 57 / 210, 1e-12)
    ok("jackknife_sign: the sign of 1..20 survives every deletion", _jk.get("sign_survives_all"))
    _jk = QF.jackknife_sign([-1.0] * 19 + [100.0])
    ok("jackknife_sign: one outlier carrying the mean is caught",
       not _jk.get("sign_survives_all") and _jk.get("worst_mean", float("nan")) == -1.0)

    # The whole gate, both directions, on base strategies whose own mean is not zero.
    def _book(seed, effect, base, n=2000):
        g = _np.random.default_rng(seed)
        sc = g.integers(0, 4, n).astype(float)
        pn = base + effect * (sc - 1.5) + g.normal(0, 1, n)
        masks = _np.array([sc >= k for k in range(4)] + [g.random(n) < 0.5 for _ in range(6)])
        rs = g.integers(0, 4, n).astype(float)
        return sc, pn, masks, (rs, base + effect * (rs - 1.5) + g.normal(0, 1, n))

    _sc, _pn, _mk, _rep = _book(1, 0.25, 1.0)
    _cg = QF.confluence_gate(_sc, _pn, masks=_mk, replication=_rep, n_rep=200)
    ok("confluence_gate: a planted effect on a profitable base passes all five gates",
       all(_cg["gates"].values()) and _cg["verdict"].startswith("UTILIZZABILE"), str(_cg["gates"]))
    _cg = QF.confluence_gate(_sc, _pn, masks=_mk, n_rep=200)
    ok("confluence_gate: a gate not run gives NON DECIDIBILE, never a pass",
       _cg["verdict"].startswith("NON DECIDIBILE") and "4_replicates" in _cg["verdict"], _cg["verdict"])
    _sc, _pn, _mk, _rep = _book(2, 0.0, -0.3)
    _cg = QF.confluence_gate(_sc, _pn, masks=_mk, replication=_rep, n_rep=200)
    ok("confluence_gate: noise on a losing base does not beat its own null",
       _cg["gates"]["2_beats_own_null"] is False and _cg["verdict"].startswith("NON USARE"), str(_cg["gates"]))


# ══════════════════════════════════ qregime: no edge works forever (§5c)
if _np is None:
    SKIP.append(f"qregime: {_no_np}")
else:
    import statistics as _st
    import qregime as QR

    _bs = QR.block_shuffle(_np.arange(103, dtype=float), 10, _np.random.default_rng(0))
    ok("block_shuffle: same values, none lost or repeated (103 = 10 blocks + 3)",
       sorted(_bs) == list(range(103)))
    _bs = QR.block_shuffle(_np.arange(100, dtype=float), 10, _np.random.default_rng(0))
    ok("block_shuffle: every block of 10 stays whole and in order",
       all(_bs[i] % 10 == 0 and list(_bs[i:i + 10]) == list(range(int(_bs[i]), int(_bs[i]) + 10))
           for i in range(0, 100, 10)))
    ok("block_shuffle: the order of the blocks changes", list(_bs) != list(range(100)))

    _g = _np.random.default_rng(7)
    _xs = _g.normal(0, 1, 700)
    _xs[300:400] += 0.6
    _ws = QR.window_scan(_xs, 100, n_null=50)
    _ts = [_st.mean(_xs[i:i + 100]) / (_st.stdev(_xs[i:i + 100]) / 10) for i in range(0, 601, 25)]
    close("window_scan: best window t matches a stdlib scan", _ws.get("best", float("nan")), max(_ts), 1e-9)
    ok("window_scan: the same windows, one every w/4", _ws.get("n_win", float("nan")) == len(_ts), f"{_ws.get('n_win')} vs {len(_ts)}")
    close("window_scan: share of windows above t = 2", _ws.get("frac_good", float("nan")), sum(t > 2 for t in _ts) / len(_ts), 1e-12)
    ok("window_scan: too short is not a verdict", QR.window_scan(_xs[:250], 100).get("ok") is False)

    # P3 by hand: +1 for 500 days, then -1 for 500, trailing 250-day mean. The switch stays on
    # until day 625, where the window holds as many losses as gains: 250 gains and 125 losses
    # selected, mean 1/3, against -1/3 always on. A window that includes today changes both.
    _x = _np.array([1.0] * 500 + [-1.0] * 500)
    _sw = QR.regime_switch(_x, look=250)
    ok("regime_switch: on for 375 of 750 days", _sw.get("duty", float("nan")) == 0.5, str(_sw.get("duty", float("nan"))))
    close("regime_switch: the trailing rule earns 1/3 a day", _sw.get("mean_rule", float("nan")), 1 / 3, 1e-12)
    close("regime_switch: always on earns -1/3 a day", _sw.get("mean_always", float("nan")), -1 / 3, 1e-12)
    ok("regime_switch: beats always-on and the random switch", _sw.get("beats_always") and _sw.get("beats_random"))
    _sw = QR.regime_switch(_x, look=250, detector=[1.0 if i < 400 else -1.0 for i in range(1000)])
    ok("regime_switch: a given detector is used as is (on days 250-399 only)",
       _sw.get("duty", float("nan")) == 0.2 and _sw.get("mean_rule", float("nan")) == 1.0, f"duty {_sw.get('duty')}, mean {_sw.get('mean_rule')}")
    ok("regime_switch: too short is not a verdict", QR.regime_switch(_x[:300], look=250).get("ok") is False)

    # The module's own power check: an edge confined to a quarter of the sample must come back
    # real and identifiable in advance; pure noise must not.
    _vp = QR.validate_power()
    ok("validate_power: a regime-limited edge beats the block null (P2)", _vp["regime_edge"]["P2_more_than_chance"])
    ok("validate_power: and is identifiable in advance (P3)", _vp["regime_edge"]["P3_identifiable_ex_ante"])
    ok("validate_power: pure noise does not beat the block null", not _vp["no_edge"]["P2_more_than_chance"])
    ok("regime_report: the verdict follows P2 and P3",
       _vp["regime_edge"]["verdict"].startswith("ROTAZIONE SOSTENUTA")
       and not _vp["no_edge"]["verdict"].startswith(("ROTAZIONE", "REGIMI REALI")), _vp["no_edge"]["verdict"])

# ══════════════════════════════════ summary
print(f"\n{'='*66}\nquant-rigor test suite\n{'='*66}")
for f in FAIL:
    print(f"  FAIL  {f}")
for s in SKIP:
    print(f"  SKIP  {s}")

print(f"\n  passed {len(PASS)} · failed {len(FAIL)} · sections skipped {len(SKIP)}")
print(f"  {'ALL GREEN' if not FAIL else 'SUITE RED — fix before shipping'}")
sys.exit(1 if FAIL else 0)
