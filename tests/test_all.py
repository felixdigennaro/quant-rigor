#!/usr/bin/env python3
"""quant-rigor test suite. Run: python3 tests/test_all.py

STDLIB ONLY, except three sections that need numpy (ndarray inputs, qconfluence, qregime): they
run when numpy is installed and are listed as SKIP when it is not. Exits 1 if any check fails.

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
with open(os.path.abspath(__file__), encoding="utf-8") as _fh:
    _last = ast.unparse(ast.parse(_fh.read()).body[-1])
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

    # The whole gate on base strategies whose own mean is not zero. The score is drawn before
    # the P&L, so it is pre-registered: gate 2 keeps it and shuffles only the P&L.
    def _book(seed, effect, base, n=2000):
        g = _np.random.default_rng(seed)
        sc = g.integers(0, 4, n).astype(float)
        pn = base + effect * (sc - 1.5) + g.normal(0, 1, n)
        rs = g.integers(0, 4, n).astype(float)
        return sc, pn, (rs, base + effect * (rs - 1.5) + g.normal(0, 1, n))

    def _gate(*a, **kw):
        try:
            return QF.confluence_gate(*a, **kw)
        except Exception as e:  # a raise is a FAIL on the check that reads it, not a crash
            return {"error": repr(e)}

    def _gate2(res):
        return res.get("gates", {}).get("2_beats_own_null", res.get("error", "missing"))

    def _null(res, key):
        return res.get("search_null", {}).get(key, float("nan"))

    _sc, _pn, _rep = _book(1, 0.25, 1.0)
    _cg = _gate(_sc, _pn, replication=_rep, preregistered=True)
    ok("confluence_gate: a planted effect on a profitable base passes all five gates",
       all(_cg.get("gates", {}).values()) and _cg.get("verdict", "").startswith("UTILIZZABILE"),
       str(_cg.get("gates")))
    _cg = _gate(_sc, _pn, preregistered=True)
    ok("confluence_gate: a gate not run gives NON DECIDIBILE, never a pass",
       _cg.get("verdict", "").startswith("NON DECIDIBILE") and "4_replicates" in _cg.get("verdict", ""),
       _cg.get("verdict"))

    # The real statistic against the HC3 slope t in plain Python (DiCiccio & Romano 2017):
    # z the centred score, b = Σz(y-ȳ)/Σz², e the residuals, h = 1/n + z²/Σz², and
    # var(b) = Σ z²e²/(1-h)² / (Σz²)². The plain OLS t assumed the P&L's variance ignores the
    # score, and noise whose variance rose with the score passed gate 2 at 17% (see below).
    def _hc3_t(x, y):
        n = len(x)
        mx, my = sum(x) / n, sum(y) / n
        z = [a - mx for a in x]
        szz = sum(a * a for a in z)
        b = sum(a * (v - my) for a, v in zip(z, y)) / szz
        e = [v - my - b * a for a, v in zip(z, y)]
        var = sum(a * a * r * r / (1 - 1 / n - a * a / szz) ** 2 for a, r in zip(z, e)) / szz ** 2
        return b / math.sqrt(var)

    def _p_hand(res):
        """the permutation p-value recomputed from the exposed draws"""
        sn = res.get("search_null", {})
        a = [float(v) for v in sn.get("null_abs_t", [])]
        return (1 + sum(v >= sn.get("real_abs_t", float("inf")) for v in a)) / (1 + len(a))

    def _rule(res):
        """gate 2 is (p <= 0.05); p is the permutation p-value of the exposed draws; q95 and
        median are those draws' percentiles (numpy's linear interpolation). Pre-registered,
        the draws are the n_rep shuffles. Replayed, they are the shuffles on which the search
        found something: n_rep of them, or at least 19 when 10·n_rep shuffles ran out first,
        and every shuffle tried is either a draw or a constant rebuild"""
        sn = res.get("search_null", {})
        a = sorted(float(v) for v in sn.get("null_abs_t", []))
        n_rep = sn.get("n_rep", 0)
        if not a:
            return False
        if sn.get("method") == "replay":
            drawn, const = sn.get("n_drawn", -1), sn.get("n_rebuilt_constant", -1)
            if const < 0 or drawn != len(a) + const or not (
                    (len(a) == n_rep and drawn <= 10 * n_rep)
                    or (19 <= len(a) < n_rep and drawn == 10 * n_rep)):
                return False
        elif len(a) != n_rep:
            return False
        pos = 0.95 * (len(a) - 1)
        q95 = a[int(pos)] + (pos - int(pos)) * (a[int(pos) + 1] - a[int(pos)])
        return (abs(_p_hand(res) - sn.get("p", float("nan"))) < 1e-12
                and _gate2(res) is (sn.get("p", 1.0) <= 0.05)
                and abs(q95 - sn.get("q95", float("nan"))) < 1e-12
                and abs(_st.median(a) - sn.get("median", float("nan"))) < 1e-12)

    def _not_run(res, word):
        """gate 2 not run, with a note that names the reason"""
        return _gate2(res) is None and word in res.get("search_null", {}).get("note", "")

    _sc, _pn, _rep = _book(3, 0.1, 0.5)
    _pre = _gate(_sc, _pn, preregistered=True)
    close("confluence_gate: real_abs_t is the HC3 |t| of the P&L on the score",
          _null(_pre, "real_abs_t"), abs(_hc3_t(_sc.tolist(), _pn.tolist())), 1e-9)
    ok("confluence_gate: gate 2 is p <= 0.05 with p, q95 and median recomputed from null_abs_t (pre-registered)",
       _rule(_pre), f"p {_null(_pre, 'p')}, gate {_gate2(_pre)}")
    # Every way gate 2 cannot run is None with a note that names the reason, never a verdict.
    # One NaN in the P&L used to come back as a degenerate null blamed on variance or block.
    _cg = _gate(_np.full(2000, 2.0), _pn, preregistered=True)
    ok("confluence_gate: a constant score leaves gate 2 not run, with a note on variance",
       _not_run(_cg, "senza varianza"), str(_cg.get("search_null")))
    _cg = _gate(_sc, _np.full(2000, 1.0), preregistered=True)
    ok("confluence_gate: a constant P&L leaves gate 2 not run, with a note on variance",
       _not_run(_cg, "senza varianza"), str(_cg.get("search_null")))
    _cg = _gate(_sc, _pn, preregistered=True, block=5000)
    ok("confluence_gate: a block longer than the series makes every shuffle the series itself, so the null ties the real and gate 2 is not run",
       _not_run(_cg, "pareggia"), str(_cg.get("search_null")))
    _bad = _pn.copy()
    _bad[10] = float("nan")
    _cg = _gate(_sc, _bad, preregistered=True)
    ok("confluence_gate: a NaN in the P&L leaves gate 2 not run, with a note naming non-finite values",
       _not_run(_cg, "non finiti") and _cg.get("verdict", "").startswith("NON DECIDIBILE"),
       f"{_cg.get('search_null')} / {_cg.get('verdict')}")
    _bad[10] = float("inf")
    _cg = _gate(_sc, _bad, preregistered=True)
    ok("confluence_gate: an inf in the P&L leaves gate 2 not run, with a note naming non-finite values",
       _not_run(_cg, "non finiti"), str(_cg.get("search_null")))
    _bad = _sc.copy()
    _bad[5] = float("inf")
    with _np.errstate(invalid="ignore"):   # the sizing's z-score of an inf score warns; not under test
        _cg = _gate(_bad, _pn, preregistered=True)
    ok("confluence_gate: an inf in the score leaves gate 2 not run, with a note naming non-finite values",
       _not_run(_cg, "non finiti"), str(_cg.get("search_null")))
    _cg = _gate(_pn, _pn, preregistered=True)
    ok("confluence_gate: a P&L that is an exact linear function of the score has no slope t, so gate 2 is not run",
       _not_run(_cg, "non è definito"), str(_cg.get("search_null")))
    # 10 periods in 5 blocks of 2, the P&L the score with its first two blocks swapped: the
    # draw that swaps them back (index 138 of 200 at seed 0) fits the score exactly, in
    # integer arithmetic, and has no t. Pre-registered, that draw is not a 0 either.
    _cg = _gate(_np.arange(10, dtype=float),
                _np.array([2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
                preregistered=True, block=2)
    ok("confluence_gate: a pre-registered draw whose P&L fits the score exactly has no slope t: gate 2 is not run, with a note naming the draw",
       _not_run(_cg, "non è definito") and _not_run(_cg, "estrazione"), str(_cg.get("search_null")))
    _cg = _gate(_sc, _pn, preregistered=True, n_rep=18)
    ok("confluence_gate: n_rep=18 cannot reach p = 0.05, so gate 2 is not run with a note naming n_rep",
       _not_run(_cg, "n_rep"), str(_cg.get("search_null")))

    # The pre-registered null on iid noise, shuffled one period at a time, is |t(1998)|:
    # median 0.6745, 95th percentile 1.961. The deviation has two parts: the Monte Carlo error
    # of 1000 draws (0.025 and 0.06) and the dataset's own permutation distribution, which
    # moves with its 2000 values. Measured over seeds 2-11 (2026-10-25): |Δmedian| up to 0.067,
    # |Δq95| up to 0.123. The tolerances sit above both; the rule check pins the record's
    # median and q95 to the exposed draws exactly.
    _g = _np.random.default_rng(2)
    _cg = _gate(_g.integers(0, 4, 2000).astype(float), _g.normal(0, 1, 2000),
                preregistered=True, n_rep=1000, block=1)
    close("confluence_gate: the pre-registered null's median |t| on iid noise = 0.6745",
          _null(_cg, "median"), _t_quantile(0.75, 1998), 0.1)
    close("confluence_gate: the pre-registered null's q95 of |t| on iid noise = 1.961",
          _null(_cg, "q95"), _t_quantile(0.975, 1998), 0.2)

    # Pre-registered noise on a losing base, 20 series: gate 2 runs on every one and passes
    # at most 4 (nominal 10/201; measured 2 of these 20 and 13 of 200 on seeds 100-299,
    # 2026-10-25; a 5% rate exceeds 4 of 20 with probability 0.26%). The calibration checks
    # above pin the null's scale; this one pins that the verdict follows it on a base whose
    # mean is not zero.
    _v = [_gate2(_gate(*_book(s, 0.0, -0.3)[:2], preregistered=True))
          for s in range(100, 120)]
    ok("confluence_gate: pre-registered noise on a losing base passes gate 2 at most 4 times in 20",
       all(isinstance(v, bool) for v in _v) and sum(v is True for v in _v) <= 4, str(_v))

    # A score the SEARCH built: the 3 of 100 random 50% masks with the highest one-sample t
    # of the P&L on the first 60% of periods, score = how many of the 3 are on. Replaying the
    # search on the shuffled P&L puts the selection bias into the null (q95 about 4, not 2),
    # so noise passes at most 5 of 20 (nominal 5/101, exceeded with probability 0.03%;
    # measured 2 of these 20, 18 of 200 on seeds 200-399 and 24 of 400 on seeds 1000-1399,
    # and at n_rep=400 the p-value is uniform: 7, 14, 28 and 54 of 100 below 0.05, 0.10,
    # 0.25 and 0.50; 2026-10-25).
    def _search3(masks, pnl):
        sub = masks[:, :1200].astype(float)
        p = pnl[:1200]
        cnt = sub.sum(1)
        mean = (sub @ p) / cnt
        var = ((sub @ (p ** 2)) / cnt - mean ** 2) * cnt / (cnt - 1)
        return masks[_np.argsort(mean / _np.sqrt(var / cnt))[-3:]].sum(0).astype(float)

    def _searched(seed, effect, n=2000):
        g = _np.random.default_rng(seed)
        masks = g.random((100, n)) < 0.5
        pn = g.normal(0, 1, n) + effect * masks[:3].sum(0)
        return masks, pn, _search3(masks, pn)

    _v = []
    for _s in range(200, 220):
        _mk, _pn, _sc = _searched(_s, 0.0)
        _v.append(_gate2(_gate(_sc, _pn, masks=_mk, build_score=_search3, n_rep=100)))
    ok("confluence_gate: a searched noise score with its search replayed passes gate 2 at most 5 times in 20",
       all(isinstance(v, bool) for v in _v) and sum(v is True for v in _v) <= 5, str(_v))

    _mk, _pn, _sc = _searched(300, 0.3)
    _rep_ = _gate(_sc, _pn, masks=_mk, build_score=_search3, n_rep=100)
    ok("confluence_gate: a planted effect on the 3 masks the search finds beats the replayed null",
       _gate2(_rep_) is True, f"|t| {_null(_rep_, 'real_abs_t'):.2f} vs q95 {_null(_rep_, 'q95'):.2f}")
    ok("confluence_gate: gate 2 is p <= 0.05 with p, q95 and median recomputed from null_abs_t (replay)",
       _rule(_rep_), f"p {_null(_rep_, 'p')}, gate {_gate2(_rep_)}")
    ok("confluence_gate: the search_null record names its method",
       _null(_rep_, "method") == "replay" and _null(_pre, "method") == "preregistered"
       and _null(_rep_, "two_sided") is True and _null(_rep_, "block") == 25,
       f"{_null(_rep_, 'method')} / {_null(_pre, 'method')}, two_sided {_null(_rep_, 'two_sided')}")
    _cg = _gate(_sc, _pn, masks=_mk, n_rep=100)
    ok("confluence_gate: a searched score with neither build_score nor preregistered leaves gate 2 not run",
       _gate2(_cg) is None and _cg.get("verdict", "").startswith("NON DECIDIBILE")
       and "2_beats_own_null" in _cg.get("verdict", ""), _cg.get("verdict"))
    _cg = _gate(_sc, _pn, masks=_mk, n_rep=100, build_score=lambda m, p: _search3(m, p)[::-1])
    ok("confluence_gate: a build_score that does not rebuild the score leaves gate 2 not run",
       _not_run(_cg, "non ricostruisce"), str(_cg.get("search_null")))
    _cg = _gate(_sc, _pn, masks=_mk, n_rep=100, build_score=lambda m, p: _search3(m, p)[:5])
    ok("confluence_gate: a build_score of the wrong shape on the real P&L leaves gate 2 not run",
       _not_run(_cg, "non ricostruisce"), str(_cg.get("search_null")))
    _cg = _gate(_sc, _pn, masks=_mk, n_rep=100,
                build_score=lambda m, p: _sc if _np.array_equal(p, _pn) else _np.ones(5))
    ok("confluence_gate: a rebuilt score of the wrong shape on a shuffle leaves gate 2 not run, with a note naming the draw",
       _not_run(_cg, "forma") and _not_run(_cg, "estrazione"), str(_cg.get("search_null")))
    _cg = _gate(_sc, _pn, masks=_mk, n_rep=100,
                build_score=lambda m, p: _sc if _np.array_equal(p, _pn) else p)
    ok("confluence_gate: a shuffle whose P&L fits the rebuilt score exactly has no slope t: gate 2 is not run, that draw is not a 0",
       _not_run(_cg, "non è definito") and _not_run(_cg, "estrazione"), str(_cg.get("search_null")))
    # A builder that ignores the P&L keeps the searched score fixed against a null that only
    # shuffles the P&L, which is the pre-registered null with the selection bias left in: on
    # a searched noise score it passed 185 of 200 as "replay".
    _mk, _pn, _sc = _searched(300, 0.0)
    _cg = _gate(_sc, _pn, masks=_mk, n_rep=100, build_score=lambda m, p: _sc)
    ok("confluence_gate: a build_score that ignores the P&L does not replay the search: gate 2 is not run and the note points at preregistered=True",
       _not_run(_cg, "preregistered=True"), str(_cg.get("search_null")))

    # build_score gets a copy of the P&L: a builder that centred in place moved the caller's
    # mean to 0 and shifted every later gate.
    def _centring(masks, pnl):
        pnl -= pnl.mean()
        return _search3(masks, pnl)

    _mk, _pn, _ = _searched(310, 0.3)
    _sc = _centring(_mk, _pn.copy())
    _was = _pn.copy()
    _cg = _gate(_sc, _pn, masks=_mk, build_score=_centring, n_rep=100)
    ok("confluence_gate: a builder that centres in place leaves the caller's P&L as it was",
       _np.array_equal(_pn, _was) and _gate2(_cg) is True,
       f"mean before {_was.mean():.3f}, after {_pn.mean():.3f}, gate {_gate2(_cg)}")
    # Centring leaves the slope t unchanged, so it cannot tell whether each SHUFFLE was
    # copied too: a builder that clips its argument in place must get the same null as the
    # pure builder that clips a copy (without copies the draws differ by up to 1.6).
    def _clipping(masks, pnl):
        pnl[:] = _np.minimum(pnl, 0.0)
        return _search3(masks, pnl)

    def _clipped(masks, pnl):
        return _search3(masks, _np.minimum(pnl, 0.0))

    _mk, _pn, _ = _searched(330, 0.3)
    _sc = _clipped(_mk, _pn)
    _a = _gate(_sc, _pn, masks=_mk, build_score=_clipping, n_rep=100)
    _b = _gate(_sc, _pn, masks=_mk, build_score=_clipped, n_rep=100)
    ok("confluence_gate: a builder that clips its P&L argument in place gets the same null as the pure builder, so every shuffle's t is computed on the untouched shuffle",
       _gate2(_a) is not None and _np.array_equal(_a.get("search_null", {}).get("null_abs_t", [0.0]),
                                                  _b.get("search_null", {}).get("null_abs_t", [1.0])),
       f"p {_null(_a, 'p')} vs {_null(_b, 'p')}")

    # build_score gets a copy of masks too: a builder that marked its picks as used in place
    # emptied all 100 of the caller's rules, and every rebuild after the first was constant.
    # A list of rows is deep-copied, since a shallow copy shares the rows.
    def _marking(masks, pnl):
        sc = _search3(masks, pnl)
        masks[:] = False
        return sc

    def _marking_rows(rows, pnl):
        sc = _search3(_np.stack(rows), pnl)
        rows[0][:] = False
        return sc

    _mk, _pn, _sc = _searched(300, 0.3)
    _was = _mk.copy()
    _cg = _gate(_sc, _pn, masks=_mk, build_score=_marking, n_rep=100)
    ok("confluence_gate: a builder that empties its masks argument in place leaves the caller's masks as they were, and the gate runs on them",
       _np.array_equal(_mk, _was) and _gate2(_cg) is True,
       f"rules on in {int(_mk.sum())} of {int(_was.sum())} cells, gate {_gate2(_cg)}")
    _rows = [_row.copy() for _row in _mk]
    _cg = _gate(_sc, _pn, masks=_rows, build_score=_marking_rows, n_rep=100)
    ok("confluence_gate: masks given as a list of rows are deep-copied: a builder that empties the first row in place leaves the caller's row",
       _np.array_equal(_np.stack(_rows), _was) and _gate2(_cg) is True,
       f"row 0 on in {int(_rows[0].sum())} of {int(_was[0].sum())} cells, gate {_gate2(_cg)}")

    # A threshold search: the score counts the rules whose on-vs-off t on the train part is
    # above a threshold. On shuffled noise it often finds nothing and the rebuilt score is
    # constant: that shuffle is not a null draw. Counted as |t| = 0 it let noise through on
    # the datasets where the real search had found something, the only ones the gate is run
    # on: 89 of 928 runs at threshold 2.5, 194 of 503 at 3.0, 191 of 191 at 3.5 (independent
    # re-review, 2026-10-25). The gate draws until n_rep shuffles have a score with variance,
    # within 10·n_rep, and the record counts both. At threshold 3 the search finds something
    # on 12-14% of shuffled noise P&Ls (1000 shuffles on each of seeds 400-404, 2026-10-25);
    # on seed 320 the 100 draws took 873 shuffles, 773 of them constant.
    def _thresh(masks, pnl, thr=3.0):
        on = masks[:, :1200].astype(float)
        off = 1.0 - on
        p = pnl[:1200]
        n1, n0 = on.sum(1), off.sum(1)
        m1, m0 = (on @ p) / n1, (off @ p) / n0
        v1 = ((on @ (p ** 2)) / n1 - m1 ** 2) * n1 / (n1 - 1)
        v0 = ((off @ (p ** 2)) / n0 - m0 ** 2) * n0 / (n0 - 1)
        return masks[(m1 - m0) / _np.sqrt(v1 / n1 + v0 / n0) > thr].sum(0).astype(float)

    _mk, _pn, _ = _searched(320, 0.3)
    _sc = _thresh(_mk, _pn)
    _cg = _gate(_sc, _pn, masks=_mk, build_score=_thresh, n_rep=100)
    _sn = _cg.get("search_null", {})
    ok("confluence_gate: a threshold search with a planted effect passes gate 2 on 100 found draws, the shuffles that find nothing skipped and counted",
       _sc.max() == 3 and _gate2(_cg) is True and len(_sn.get("null_abs_t", [])) == 100
       and _sn.get("n_rebuilt_constant", 0) > 0 and _rule(_cg),
       f"found {int(_sc.max())} rules, gate {_gate2(_cg)}, |t| {_null(_cg, 'real_abs_t'):.2f}, "
       f"p {_null(_cg, 'p')}, draws {len(_sn.get('null_abs_t', []))} of {_sn.get('n_drawn')} shuffles, "
       f"constant {_sn.get('n_rebuilt_constant')}")

    # At threshold 3.5 the search finds something on about 2.3% of shuffled noise P&Ls (19-30
    # of 1000 on seeds 400-404) and on 3 of 200 noise datasets (seeds 400-599, 2026-10-25).
    # Seed 438 is one of the 3: the real search found one rule. With n_rep=50 its 500
    # shuffles yield 11 found draws (8-13 over gate seeds 0-9), under the 19 the p-value
    # needs: the gate is not run. With n_rep=100 the 1000 shuffles yield 16-26: not run on 3
    # gate seeds of 10, p 0.48-0.64 on the other 7, never a pass. With the constant rebuilds
    # counted as |t| = 0 it passed.
    _mk, _pn, _ = _searched(438, 0.0)
    _sc = _thresh(_mk, _pn, 3.5)
    _cg = _gate(_sc, _pn, masks=_mk, build_score=lambda m, p: _thresh(m, p, 3.5), n_rep=50)
    ok("confluence_gate: a strict threshold search on noise that found a rule, with too few shuffles on which it finds anything, leaves gate 2 not run with a note on the count",
       _sc.max() == 1 and _not_run(_cg, "troppo pochi"),
       f"found {int(_sc.max())} rules, gate {_gate2(_cg)}, {_cg.get('search_null')}")

    # Conditional calibration, the measurement that exposed the defect: threshold 2.5 on
    # noise, only the datasets where the real search found something. The first 20 such
    # datasets from seed 500 (seeds 501-542) pass 0; 200 of them (seeds 501-895) pass 8,
    # at the nominal 5/101 (2026-10-25). A calibrated gate exceeds 5 of 20 with probability
    # 0.03%. With the constant rebuilds counted as |t| = 0 the re-review measured 9.6%.
    _v, _s = [], 500
    while len(_v) < 20:
        _mk, _pn, _ = _searched(_s, 0.0)
        _sc = _thresh(_mk, _pn, 2.5)
        if _sc.max() > 0:
            _v.append(_gate2(_gate(_sc, _pn, masks=_mk, n_rep=100,
                                   build_score=lambda m, p: _thresh(m, p, 2.5))))
        _s += 1
    ok("confluence_gate: a threshold search on noise, on the datasets where it found something, passes gate 2 at most 5 times in 20",
       all(isinstance(v, bool) for v in _v) and sum(v is True for v in _v) <= 5,
       f"{_v}, last seed {_s - 1}")

    # Noise whose variance rises with a skewed 0-3 count (sd 1, 1.2, 1.6, 2.5 at p .5, .3,
    # .15, .05), pre-registered, mean 0. The plain OLS slope t is inflated on it: 68 of 400
    # passed gate 2 (17%); with the HC3 t 22 of 400 (5.5%; seeds 10000-10399, 2026-10-25).
    # On these 60 seeds OLS passes 10 and HC3 5; a calibrated 5% gate exceeds 7 of 60 with
    # probability 0.98%.
    _sd = _np.array([1.0, 1.2, 1.6, 2.5])
    _v = []
    for _s in range(60):
        _g = _np.random.default_rng(10_000 + _s)
        _sc = _g.choice([0.0, 1.0, 2.0, 3.0], size=2000, p=[0.5, 0.3, 0.15, 0.05])
        _v.append(_gate2(_gate(_sc, _g.normal(0, 1, 2000) * _sd[_sc.astype(int)], preregistered=True)))
    ok("confluence_gate: noise whose variance rises with the score passes gate 2 at most 7 times in 60",
       all(isinstance(v, bool) for v in _v) and sum(v is True for v in _v) <= 7,
       f"{sum(v is True for v in _v)} passes, {sum(v is None for v in _v)} not run")
    try:
        QF.confluence_gate(_sc, _pn, masks=_mk, build_score=_search3, preregistered=True)
        _err = None
    except Exception as e:  # a wrong exception type is a FAIL, not a crash
        _err = e
    ok("confluence_gate: build_score with preregistered=True is refused",
       isinstance(_err, ValueError), repr(_err))

    _sc, _pn, _rep = _book(4, -0.25, 1.0)
    _cg = _gate(_sc, _pn, preregistered=True)
    ok("confluence_gate: a falling dose-response beats the two-sided null",
       _gate2(_cg) is True, f"|t| {_null(_cg, 'real_abs_t'):.2f} vs q95 {_null(_cg, 'q95'):.2f}")

    # The rule at its edges. n_rep=19 is the smallest null whose p-value reaches 0.05: a
    # planted effect above every draw has p = 1/20 and passes. On noise seed 406 one draw of
    # 29 reaches the real: p = 2/30 fails, where a bar at the 93rd percentile of the draws
    # would pass. With 2 blocks the series itself is half the shuffles and ties the real, so
    # even an evident effect cannot pass and the gate is not run; on seed 5 the batch kernel
    # rounds those draws one ulp below the single-row real, so only the tolerance sees the
    # ties. With 5 blocks it is 1 shuffle in 120; on seed 11 three of 200 are the series
    # itself, rounded the same way (2026-10-25): with exact equality p would be 1/201, not
    # 4/201.
    _sc, _pn = _book(1, 0.25, 1.0)[:2]
    _cg = _gate(_sc, _pn, preregistered=True, n_rep=19)
    ok("confluence_gate: n_rep=19 with the real above every draw is p = 0.05 exactly, which passes",
       _null(_cg, "p") == 0.05 and _gate2(_cg) is True and _rule(_cg),
       f"p {_null(_cg, 'p')}, gate {_gate2(_cg)}")
    _cg = _gate(*_book(406, 0.0, -0.3)[:2], preregistered=True, n_rep=29)
    ok("confluence_gate: one draw of 29 at or above the real is p = 2/30, which fails",
       _null(_cg, "p") == 2 / 30 and _gate2(_cg) is False and _rule(_cg),
       f"p {_null(_cg, 'p')}, gate {_gate2(_cg)}")
    _cg = _gate(*_book(5, 0.25, 1.0)[:2], preregistered=True, block=1000)
    ok("confluence_gate: with 2 blocks the null ties an evident effect on half its draws, so gate 2 is not run, never failed",
       _not_run(_cg, "pareggia"), str(_cg.get("search_null")))
    _cg = _gate(*_book(11, 0.25, 1.0)[:2], preregistered=True, block=400, seed=11)
    _real = _null(_cg, "real_abs_t")
    _ties = sum(abs(float(v) - _real) <= 1e-9 * _real
                for v in _cg.get("search_null", {}).get("null_abs_t", []))
    ok("confluence_gate: with 5 blocks the shuffles that are the series itself tie the real and count as >=",
       _ties >= 1 and _null(_cg, "p") == (1 + _ties) / 201 and _gate2(_cg) is True,
       f"ties {_ties}, p {_null(_cg, 'p')}, gate {_gate2(_cg)}")

    # A score with period 5 on daily data under block 25, five trading weeks: every block
    # shuffle keeps each P&L value on its weekday, the slope t of every draw equals the real
    # one, and an effect of |t| 14 came back p 1.0, gate 2 False, NON USARE (the re-review's
    # |t| 19.96 the same way, 2026-10-25). A null that ties the real on more than 5% of its
    # draws cannot pass whatever the data: that is the caller's block. Under block 24 the
    # same series passes at p 1/201. Replayed, a search that ranks the weekdays by their
    # train-part mean rebuilds the same ranks on 191 of 200 block-25 shuffles (the blocks
    # that cross the train cut move a few means): the same note, not the one that blames a
    # builder for ignoring the P&L.
    _wd = (_np.arange(2000) % 5).astype(float)
    _g = _np.random.default_rng(77)
    _wpn = 0.25 * (_wd - 2.0) + _g.normal(0, 1, 2000)
    _cg = _gate(_wd, _wpn, preregistered=True, block=25)
    ok("confluence_gate: a period-5 score under block 25 ties the real t on every draw: gate 2 is not run, with a note naming block",
       _not_run(_cg, "pareggia") and _not_run(_cg, "block")
       and _cg.get("verdict", "").startswith("NON DECIDIBILE"), str(_cg.get("search_null")))
    _cg = _gate(_wd, _wpn, preregistered=True, block=24)
    ok("confluence_gate: the same effect under block 24 passes gate 2",
       _gate2(_cg) is True and _rule(_cg), f"p {_null(_cg, 'p')}, gate {_gate2(_cg)}")

    def _weekday_rank(masks, pnl):
        wd = _np.arange(len(pnl)) % 5
        means = _np.array([pnl[:1200][wd[:1200] == w].mean() for w in range(5)])
        rank = _np.empty(5)
        rank[_np.argsort(means)] = _np.arange(5)
        return rank[wd]

    _cg = _gate(_wd, _wpn, build_score=_weekday_rank, block=25)
    ok("confluence_gate: a replayed weekday search under block 25 gets the tie note, not the one on a builder that ignores the P&L",
       _not_run(_cg, "pareggia") and not _not_run(_cg, "preregistered=True"),
       str(_cg.get("search_null")))
    _cg = _gate(_wd, _wpn, build_score=_weekday_rank, block=24)
    ok("confluence_gate: the same replayed search under block 24 passes gate 2",
       _gate2(_cg) is True and _rule(_cg), f"p {_null(_cg, 'p')}, gate {_gate2(_cg)}")

    # Gate 5 examines the bucket the sizing would weight UP: the top-score bucket on a rising
    # relation, the bottom-score bucket on a falling one. It examined the top bucket always,
    # and on a falling relation that bucket's negative total made top3_share negative and the
    # gate trivially true: a sized-up setup whose gains lived in 3 periods passed 49 of 50
    # (lead's measurement, 2026-10-25). It also passes only if that bucket makes money. The
    # fixtures put a falling 0.1-per-level relation on noise of sd 0.1 and rewrite the
    # score-0 bucket: 3 gains of 10 in 3 different jackknife blocks, or 20 gains of 0.2
    # inside one block (periods 135-154 of its 523, block 5 of 20), or gains spread evenly.
    def _falling(seed, bucket0):
        g = _np.random.default_rng(seed)
        sc = g.integers(0, 4, 2000).astype(float)
        pn = -0.1 * (sc - 1.5) + g.normal(0, 0.1, 2000)
        low = sc == 0.0
        if bucket0 == "three gains":
            pn[low] = 0.02 + g.normal(0, 0.1, low.sum())
            pn[_np.flatnonzero(low)[[5, 200, 400]]] = 10.0
        elif bucket0 == "one block":
            pn[low] = -0.002 + g.normal(0, 0.02, low.sum())
            pn[_np.flatnonzero(low)[135:155]] += 0.2
        elif bucket0 == "spread gains":
            pn[low] = 0.2 + g.normal(0, 0.1, low.sum())
        elif bucket0 == "all buckets lose":
            pn -= 0.8
        return sc, pn

    def _gate5(res):
        return res.get("gates", {}).get("5_not_concentrated", res.get("error", "missing"))

    _cg = _gate(*_falling(5, "three gains"), preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: a falling relation whose sized-up bucket owes its gains to 3 periods fails gate 5 on that bucket",
       _gate5(_cg) is False and _c.get("bucket") == "bottom" and _c.get("bucket_mean", 0.0) > 0
       and _c.get("sign_survives_all") is True and _c.get("top3_share", 0.0) >= 0.5,
       f"gate {_gate5(_cg)}, bucket {_c.get('bucket')}, mean {_c.get('bucket_mean')}, top3 {_c.get('top3_share')}")
    _cg = _gate(*_falling(5, "one block"), preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: a sized-up bucket whose mean one jackknife block carries fails gate 5 with top3_share below 0.5",
       _gate5(_cg) is False and _c.get("bucket") == "bottom" and _c.get("bucket_mean", 0.0) > 0
       and _c.get("sign_survives_all") is False and _c.get("top3_share", 1.0) < 0.5,
       f"gate {_gate5(_cg)}, bucket {_c.get('bucket')}, mean {_c.get('bucket_mean')}, worst {_c.get('worst_mean')}, top3 {_c.get('top3_share')}")
    _cg = _gate(*_falling(5, "spread gains"), preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: the same relation with the gains spread passes gate 5 on the bottom bucket",
       _gate5(_cg) is True and _c.get("bucket") == "bottom",
       f"gate {_gate5(_cg)}, bucket {_c.get('bucket')}, mean {_c.get('bucket_mean')}, top3 {_c.get('top3_share')}")
    _cg = _gate(*_falling(5, "all buckets lose"), preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: a base that loses in every bucket fails gate 5: an edge you size up must make money",
       _gate5(_cg) is False and _c.get("bucket_mean", 0.0) < 0 and _c.get("sign_survives_all") is True,
       f"gate {_gate5(_cg)}, bucket {_c.get('bucket')}, mean {_c.get('bucket_mean')}")
    _c = _gate(_sc, _pn, preregistered=True).get("concentration", {})
    ok("confluence_gate: a rising relation is examined on its top bucket, and the record says so",
       _c.get("bucket") == "top" and _c.get("bucket_n", 0) > 0 and _c.get("bucket_mean", 0.0) > 1.0,
       f"bucket {_c.get('bucket')}, n {_c.get('bucket_n')}, mean {_c.get('bucket_mean')}")

    # Persistent data, no effect: P&L AR(0.6), score = quartile of an AR(0.97) series. The
    # slope t is inflated on such data, so the null must be too. Blocks of 25 keep the
    # persistence; an iid shuffle destroys it (q95 ratio 1.46 to 2.30 over 40 series).
    def _ar(g, phi, n):
        x = _np.empty(n)
        x[0] = g.normal()
        for i in range(1, n):
            x[i] = phi * x[i - 1] + g.normal()
        return x

    _g = _np.random.default_rng(7)
    _pn = _ar(_g, 0.6, 2000)
    _x = _ar(_g, 0.97, 2000)
    _sc = _np.searchsorted(_np.quantile(_x, [0.25, 0.5, 0.75]), _x).astype(float)
    _q25 = _null(_gate(_sc, _pn, preregistered=True, block=25), "q95")
    _q1 = _null(_gate(_sc, _pn, preregistered=True, block=1), "q95")
    ok("confluence_gate: on persistent data the block null's q95 exceeds the iid null's by 25%+",
       _q25 > 1.25 * _q1, f"block 25: {_q25:.2f}, iid: {_q1:.2f}")

    # Gate 5's second test checks the slope itself (jackknife_slope). Its reference is the
    # definition in plain Python: c_i = z_i (p_i - mean p) over the periods whose score and
    # P&L are finite, in time order; the share of Σc that its 3 largest terms in its direction
    # carry; the smallest b_k / b over 20 contiguous blocks, split as np.array_split splits
    # (the first n mod 20 blocks one period longer), the means recomputed on what is kept.
    def _slope_ref(score, pnl, blocks=20):
        pairs = [(a, b) for a, b in zip(score, pnl) if math.isfinite(a) and math.isfinite(b)]

        def ols(ps):
            ms, mp = sum(a for a, _ in ps) / len(ps), sum(b for _, b in ps) / len(ps)
            return sum((a - ms) * (b - mp) for a, b in ps) / sum((a - ms) ** 2 for a, _ in ps)

        n = len(pairs)
        ms, mp = sum(a for a, _ in pairs) / n, sum(b for _, b in pairs) / n
        c = [(a - ms) * (b - mp) for a, b in pairs]
        d = 1 if sum(c) > 0 else -1
        sizes = [n // blocks + (k < n % blocks) for k in range(blocks)]
        starts = [sum(sizes[:k]) for k in range(blocks)]
        full = ols(pairs)
        return (sum(sorted(d * v for v in c)[-3:]) / abs(sum(c)),
                min(ols(pairs[:s] + pairs[s + m:]) / full for s, m in zip(starts, sizes)))

    def _slope(*a, **kw):
        try:
            return QF.jackknife_slope(*a, **kw)
        except Exception as e:  # a raise is a FAIL on the check that reads it, not a crash
            return {"error": repr(e)}

    def _bucket_alone(c):
        """the bucket test alone, read from the record"""
        return (c.get("bucket_mean", 0.0) > 0 and c.get("sign_survives_all") is True
                and c.get("top3_share", 1.0) < 0.5)

    def _gates_1_to_4(res):
        return all(res.get("gates", {}).get(k) is True for k in
                   ("1_dose_response", "2_beats_own_null", "3_out_of_sample", "4_replicates"))

    # 103 periods with a NaN in the P&L at period 40 and an inf in the score at period 77 leave
    # 101 finite periods, split into a first block of 6 and 19 blocks of 5.
    _g = _np.random.default_rng(21)
    _sc = _g.permutation(_np.repeat([0.0, 1.0, 2.0, 3.0], 26))[:103]
    _pn = 0.4 + 0.3 * (_sc - 1.5) + _g.normal(0, 1, 103)
    _pn[40], _sc[77] = float("nan"), float("inf")
    with _np.errstate(invalid="ignore"):   # the sizing's z-score of an inf score warns; not under test
        _c = _gate(_sc, _pn, preregistered=True).get("concentration", {})
    _ref = _slope_ref(_sc.tolist(), _pn.tolist())
    close("confluence_gate: slope_top3_share is the top-3 share of the slope's numerator over the finite periods (plain Python)",
          _c.get("slope_top3_share", float("nan")), _ref[0], 1e-9)
    close("confluence_gate: slope_block_min_ratio is the smallest b_k / b over 20 contiguous blocks of the finite periods (plain Python)",
          _c.get("slope_block_min_ratio", float("nan")), _ref[1], 1e-9)
    ok("confluence_gate: the record counts 20 slope blocks", _c.get("slope_blocks") == 20,
       str(_c.get("slope_blocks")))

    # Here three crash days of -60 sit in the score-0 bucket, the one a rising sizing weights
    # DOWN, on an effect of 0.03 per level (two crashes in the first 60% of the periods, one
    # after). Seed 0
    # passes gates 1-4 with its replication and passes the bucket test, and the crashes carry
    # 69% of the slope's numerator. Its verdict was UTILIZZABILE before the slope test, as on
    # 63 of 200 such datasets in the lead's bench (2026-10-25). With the score negated the
    # relation falls, the crashes sit in the top bucket, the bottom one is examined, and the
    # numbers are the same.
    def _crashes(seed, sign):
        g = _np.random.default_rng(seed)

        def one():
            sc = g.integers(0, 4, 2000).astype(float)
            pn = 0.1 + 0.03 * (sc - 1.5) + g.normal(0, 1, 2000)
            low = _np.flatnonzero(sc == 0)
            pn[list(g.choice(low[low < 1200], 2, replace=False))
               + list(g.choice(low[low >= 1200], 1))] -= 60.0
            return sign * sc, pn
        return one(), one()

    for _sign, _side in ((1.0, "top"), (-1.0, "bottom")):
        (_sc, _pn), _rep = _crashes(0, _sign)
        _cg = _gate(_sc, _pn, replication=_rep, preregistered=True)
        _c = _cg.get("concentration", {})
        ok(f"confluence_gate: three crash days in the down-weighted bucket ({_side} bucket examined) pass gates 1-4 and the bucket test, and fail gate 5 on the slope's top-3 share",
           _gates_1_to_4(_cg) and _c.get("bucket") == _side and _bucket_alone(_c)
           and _c.get("slope_top3_share", 0.0) >= 0.5 and _gate5(_cg) is False,
           f"gates {_cg.get('gates')}, bucket {_c.get('bucket')}, slope share {_c.get('slope_top3_share')}")

    # Here the whole effect lives in one 100-period window, periods 300-399 (the fourth of 20
    # blocks). Seed 0 passes gates 1-4 with its replication and passes the bucket test, the
    # slope's top-3 share is 0.10, and deleting that block leaves 18% of the slope.
    def _one_window(seed):
        g = _np.random.default_rng(seed)

        def one():
            sc = g.integers(0, 4, 2000).astype(float)
            pn = 0.05 + g.normal(0, 1, 2000)
            pn[300:400] += 0.9 * (sc[300:400] - 1.5)
            return sc, pn
        return one(), one()

    (_sc, _pn), _rep = _one_window(0)
    _cg = _gate(_sc, _pn, replication=_rep, preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: an effect living in one 100-period window passes gates 1-4 and the bucket test, and fails gate 5 because deleting that block leaves less than half the slope",
       _gates_1_to_4(_cg) and _bucket_alone(_c) and _c.get("slope_top3_share", 1.0) < 0.5
       and _c.get("slope_block_min_ratio", 1.0) < 0.5 and _gate5(_cg) is False,
       f"gates {_cg.get('gates')}, slope share {_c.get('slope_top3_share')}, block ratio {_c.get('slope_block_min_ratio')}")

    # The rest of the sample falls (-0.1 per level) and periods 300-399 rise steeply (+4 per
    # level), so the whole sample rises. On seed 1, deleting that block leaves a slope as large
    # as the full one, reversed (ratio -1.04). The sign has to survive, not only the size.
    _g = _np.random.default_rng(1)
    _sc = _g.integers(0, 4, 2000).astype(float)
    _pn = 1.0 - 0.1 * (_sc - 1.5) + _g.normal(0, 1, 2000)
    _pn[300:400] += 4.0 * (_sc[300:400] - 1.5)
    _cg = _gate(_sc, _pn, preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: a block deletion that reverses the slope fails gate 5 even when what is left is as large as the slope",
       _c.get("bucket") == "top" and _bucket_alone(_c) and _c.get("slope_top3_share", 1.0) < 0.5
       and _c.get("slope_block_min_ratio", 0.0) <= -0.5 and _gate5(_cg) is False,
       f"gate {_gate5(_cg)}, ratio {_c.get('slope_block_min_ratio')}, share {_c.get('slope_top3_share')}")

    # These checks put each bound at its edge, with one knob on the same data (seed 1, effect
    # 0.1 per level). Three crash days in the score-0 bucket put the slope's top-3 share at
    # 0.495 with crashes of 54 and at 0.503 with crashes of 56. An extra effect in periods
    # 300-399 puts the smallest block ratio at 0.504 with 1.9 per level and at 0.492 with 2.0,
    # a ratio that keeps the sign. The other bound and the bucket test pass on all four.
    def _edge(knob, x):
        g = _np.random.default_rng(1)
        sc = g.integers(0, 4, 2000).astype(float)
        pn = 0.1 + 0.1 * (sc - 1.5) + g.normal(0, 1, 2000)
        if knob == "crash":
            pn[_np.flatnonzero(sc == 0)[[100, 250, 400]]] -= x
        else:
            pn[300:400] += x * (sc[300:400] - 1.5)
        return sc, pn

    for _name, _knob, _x, _band, _want in (
            ("a slope top-3 share of 0.495 passes gate 5", "crash", 54.0, (0.49, 0.5), True),
            ("a slope top-3 share of 0.503 fails gate 5", "crash", 56.0, (0.5, 0.51), False),
            ("a smallest block ratio of 0.504 passes gate 5", "window", 1.9, (0.5, 0.51), True),
            ("a smallest block ratio of 0.492 fails gate 5, though no deletion flips the sign",
             "window", 2.0, (0.49, 0.5), False)):
        _sc, _pn = _edge(_knob, _x)
        _share, _ratio = _slope_ref(_sc.tolist(), _pn.tolist())
        _at, _other = (_share, _ratio >= 0.5) if _knob == "crash" else (_ratio, _share < 0.5)
        _cg = _gate(_sc, _pn, preregistered=True)
        ok("confluence_gate: " + _name,
           _band[0] <= _at < _band[1] and _other and _bucket_alone(_cg.get("concentration", {}))
           and _gate5(_cg) is _want,
           f"gate {_gate5(_cg)}, share {_share:.4f}, ratio {_ratio:.4f}")

    # These two put each bound at exactly 0.5, in exact arithmetic and in floating point alike,
    # since every value is a short binary fraction. The first fixture has 128 periods with the
    # score 0-3 repeating, a P&L of 0.25 + (549/1024)·(score - 1.5), and three crash days of
    # 18.3125 at score 0. Its crashes carry exactly half of the slope's numerator, which is not
    # less than half. The second has 160 periods in 20 blocks of 8, a P&L of
    # 1 + 0.25·(score - 1.5), and 20 times that slope in periods 40-47. Deleting that block
    # leaves exactly half of the slope, which is enough.
    import fractions as _fr
    _sc = _np.tile([0.0, 1.0, 2.0, 3.0], 32)
    _pn = 0.25 + 549 / 1024 * (_sc - 1.5)
    _pn[[8, 60, 100]] -= 18.3125
    _ref = _slope_ref([_fr.Fraction(v) for v in _sc], [_fr.Fraction(v) for v in _pn])
    _cg = _gate(_sc, _pn, preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: a slope top-3 share of exactly 0.5 fails gate 5, since the 3 periods must carry less than half",
       _ref[0] == _fr.Fraction(1, 2) and _c.get("slope_top3_share") == 0.5 and _ref[1] >= 0.5
       and _bucket_alone(_c) and _gate5(_cg) is False,
       f"gate {_gate5(_cg)}, exact share {_ref[0]}, recorded {_c.get('slope_top3_share')}")
    _sc = _np.tile([0.0, 1.0, 2.0, 3.0], 40)
    _pn = 1.0 + 0.25 * (_sc - 1.5)
    _pn[40:48] += 5.0 * (_sc[40:48] - 1.5)
    _ref = _slope_ref([_fr.Fraction(v) for v in _sc], [_fr.Fraction(v) for v in _pn])
    _cg = _gate(_sc, _pn, preregistered=True)
    _c = _cg.get("concentration", {})
    ok("confluence_gate: a smallest block ratio of exactly 0.5 passes gate 5, since half of the slope is enough",
       _ref[1] == _fr.Fraction(1, 2) and _c.get("slope_block_min_ratio") == 0.5 and _ref[0] < 0.5
       and _bucket_alone(_c) and _gate5(_cg) is True,
       f"gate {_gate5(_cg)}, exact ratio {_ref[1]}, recorded {_c.get('slope_block_min_ratio')}")

    # A score that varies only in periods 400-499 has no variance left once that block is
    # deleted, so no slope is left, its ratio is NaN, and gate 5 fails. The bucket test passes
    # and the slope's top-3 share is 0.11.
    _g = _np.random.default_rng(0)
    _sc = _np.zeros(2000)
    _sc[400:500] = _g.integers(1, 4, 100)
    _cg = _gate(_sc, 0.1 + 0.5 * _sc + _g.normal(0, 1, 2000), preregistered=True)
    _c = _cg.get("concentration", {})
    _r = _c.get("slope_block_min_ratio", 0.0)
    ok("confluence_gate: a score that varies inside one block only fails gate 5, since deleting that block leaves no slope (ratio NaN)",
       _r != _r and _bucket_alone(_c) and _c.get("slope_top3_share", 1.0) < 0.5
       and _gate5(_cg) is False,
       f"gate {_gate5(_cg)}, ratio {_r}, share {_c.get('slope_top3_share')}")

    # jackknife_slope needs 3 periods per block, so 60 finite periods for 20 blocks, counted
    # after the non-finite ones are dropped. Through confluence_gate, fewer than 60 finite
    # periods leave no dose-response slope either, so gate 5 is not run before the slope test
    # is reached.
    _g = _np.random.default_rng(31)
    _sc = _np.tile([0.0, 1.0, 2.0, 3.0], 16)[:61]
    _pn = 0.5 * _sc + _g.normal(0, 1, 61)
    _js = _slope(_sc[:60], _pn[:60])
    ok("jackknife_slope: 60 finite periods, 3 per block, give a result", _js.get("ok") is True, str(_js))
    _pn[[7, 30]] = [float("nan"), float("inf")]
    _js = _slope(_sc, _pn)
    ok("jackknife_slope: 61 periods of which 59 are finite are not a verdict, with a note naming 60",
       _js.get("ok") is False and "60" in _js.get("note", ""), str(_js))
    _cg = _gate(_sc, _pn, preregistered=True)
    ok("confluence_gate: fewer than 60 finite periods leave gate 5 not run, with a note",
       _gate5(_cg) is None and bool(_cg.get("concentration", {}).get("note")),
       str(_cg.get("concentration")))

    # A U-shaped P&L, exactly symmetric in the score, has a slope numerator of exactly 0 in
    # floating point. There is no slope to attribute, and the sizing weights nothing up.
    _cg = _gate(_np.tile([0.0, 1.0, 2.0, 3.0], 100), _np.tile([1.0, -1.0, -1.0, 1.0], 100),
                preregistered=True)
    ok("confluence_gate: a slope numerator of exactly 0 leaves gate 5 not run, with a note on the numerator",
       _gate5(_cg) is None and "numeratore" in _cg.get("concentration", {}).get("note", ""),
       str(_cg.get("concentration")))

    # The rate check uses clean effects in seven families, with a median |slope t| of 5 to 8.
    # Gate 5 failed 5 of 7000 such datasets, 1000 per family on the seeds from
    # 10000·(family + 1), and none of these 70 (2026-10-25). At that rate 70 datasets fail
    # twice or more with probability 0.12%; a gate 5 that failed 5% of clean effects would
    # with probability 87%. At the lead's weaker effects (median |t| 3 to 5, 1000 per family)
    # the slope test fails up to 4.8% of clean datasets; with normal noise and on persistent
    # data, where it fails most, those datasets' realized |t| has a median of 0.7 to 1.4 and
    # none passes gates 1-3. On the lead's seven clean fixtures of 200 datasets it turned no
    # UTILIZZABILE verdict into another.
    def _clean(kind, seed, n=2000):
        g = _np.random.default_rng(seed)
        if kind == "persistent":   # P&L AR(0.6), score the quartile of an AR(0.97) series
            pn, x = _ar(g, 0.6, n), _ar(g, 0.97, n)
            sc = _np.searchsorted(_np.quantile(x, [0.25, 0.5, 0.75]), x).astype(float)
            return sc, 0.1 + 0.15 * (sc - 1.5) + 0.8 * pn
        sc = g.integers(0, 4, n).astype(float)
        if kind == "normal":
            return sc, 0.1 + 0.1 * (sc - 1.5) + g.normal(0, 1, n)
        if kind == "falling":
            return sc, 0.1 - 0.1 * (sc - 1.5) + g.normal(0, 1, n)
        if kind == "t(3)":
            return sc, 0.1 + 0.12 * (sc - 1.5) + g.standard_t(3, n) / math.sqrt(3.0)
        if kind == "t(2.5)":
            return sc, 0.1 + 0.15 * (sc - 1.5) + g.standard_t(2.5, n) / math.sqrt(5.0)
        if kind == "variance rising with the score":
            sd = _np.array([1.0, 1.2, 1.6, 2.5])[sc.astype(int)]
            return sc, 0.1 + 0.2 * (sc - 1.5) + g.normal(0, 1, n) * sd
        if kind == "two-valued (TP/SL)":
            return sc, _np.where(g.random(n) < 0.8 + 0.05 * (sc - 1.5), 1.0, -3.5)
        raise ValueError(kind)

    _kinds = ("normal", "falling", "t(3)", "t(2.5)", "variance rising with the score",
              "two-valued (TP/SL)", "persistent")
    _v = [_gate5(_gate(*_clean(_k, 10_000 * (_i + 1) + _j), preregistered=True))
          for _i, _k in enumerate(_kinds) for _j in range(10)]
    ok("confluence_gate: clean effects in seven families pass gate 5 at least 69 times in 70",
       all(isinstance(v, bool) for v in _v) and sum(v is True for v in _v) >= 69,
       f"{sum(v is True for v in _v)} passes, {sum(v is None for v in _v)} not run")


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
    # real and identifiable in advance.
    _vp = QR.validate_power()
    ok("validate_power: a regime-limited edge beats the block null (P2)", _vp["regime_edge"]["P2_more_than_chance"])
    ok("validate_power: and is identifiable in advance (P3)", _vp["regime_edge"]["P3_identifiable_ex_ante"])
    ok("regime_report: a real, identifiable regime reads ROTAZIONE SOSTENUTA",
       _vp["regime_edge"]["verdict"].startswith("ROTAZIONE SOSTENUTA"), _vp["regime_edge"]["verdict"])
    # Pure noise fires P2 by chance: two window scans at p < 0.05 make that ~10% by construction
    # (3 of 40 streams, 2026-10-25), so one stream proves nothing either way. The check is the
    # rate over 10 streams, against a bound that a 10% rate breaks about 2 times in 1000.
    _fp = sum(QR.regime_report(_np.random.default_rng(500 + s).normal(0, 1000.0, 2800), n_null=100, seed=s)
              ["P2_more_than_chance"] for s in range(10))
    ok("regime_report: pure noise fires P2 in at most 4 of 10 streams", _fp <= 4, f"{_fp}/10")

# ══════════════════════════════════ summary
print(f"\n{'='*66}\nquant-rigor test suite\n{'='*66}")
for f in FAIL:
    print(f"  FAIL  {f}")
for s in SKIP:
    print(f"  SKIP  {s}")

print(f"\n  passed {len(PASS)} · failed {len(FAIL)} · sections skipped {len(SKIP)}")
print(f"  {'ALL GREEN' if not FAIL else 'SUITE RED — fix before shipping'}")
sys.exit(1 if FAIL else 0)
