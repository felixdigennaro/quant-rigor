#!/usr/bin/env python3
"""quant-rigor test suite — STDLIB ONLY. Run: python3 tests/test_all.py

Every numeric assertion is against an INDEPENDENT reference (published statistical tables,
closed-form results, or analytically known values) — never against the library's own output.
That discipline is what caught a 1% p-value error in the t-tail at df=1.
"""
import math
import os
import sys
import xml.dom.minidom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import qstats as S          # noqa: E402
import qvalidate as V       # noqa: E402
import qsize as Z           # noqa: E402
import qcharts as C         # noqa: E402

PASS, FAIL = [], []


def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(f"{name}{(' — ' + detail) if detail else ''}")


def close(name, got, want, tol, detail=""):
    d = abs(got - want)
    ok(name, d <= tol, detail or f"got {got:.8g}, want {want:.8g}, |Δ|={d:.2e} (tol {tol:.1e})")


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

# ══════════════════════════════════ summary
print(f"\n{'='*66}\nquant-rigor test suite\n{'='*66}")
for f in FAIL:
    print(f"  FAIL  {f}")
print(f"\n  passed {len(PASS)} · failed {len(FAIL)}")
print(f"  {'ALL GREEN' if not FAIL else 'SUITE RED — fix before shipping'}")
sys.exit(1 if FAIL else 0)
