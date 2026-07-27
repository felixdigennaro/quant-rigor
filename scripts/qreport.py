#!/usr/bin/env python3
"""qreport — assemble a strategy analysis into one self-contained HTML report. STDLIB ONLY.

Design rule: the report states what was measured, under which assumptions, with which uncertainty,
and what the evidence does NOT support. No persuasive language, no adjectives that a number could
replace, no omitted caveats. Every figure carries its standard error or its band. Sections that
would be empty are omitted rather than padded.

Output is a single .html with inlined SVG and CSS — no network, no fonts, no scripts. Suitable for
publishing as an artifact, printing to PDF, or attaching to a decision record.

    from qreport import Report
    r = Report("Trailing stop evaluation", subtitle="60 synthetic sessions")
    r.assumptions(data="…", periods=43, resolution="20 bars/session",
                  costs="fee 2.0 + slippage 2.5 on market exits (measured)",
                  not_modelled=["funding", "queue position", "gaps"])
    r.measurability(bar_ranges, {"fixed stop 60": 60, "trailing 20": 20})
    r.add_series("baseline", base); r.add_series("trail 20/40", variant)
    r.compare("baseline", "trail 20/40", n_variants=12)
    r.overfitting(perf_matrix=[...], n_variants=12)
    r.sizing(mean, sd, dd_tolerated=0.15)
    r.verdict("not decidable", "unmeasurable at this resolution; gap inside best-of-N headroom",
              path_to_decide="collect 1-minute data forward; pre-register the criterion")
    r.write("report.html")
"""
import html
import math
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import qstats as S        # noqa: E402
import qcharts as C       # noqa: E402

try:
    import qvalidate as V
except Exception:                                    # optional
    V = None
try:
    import qsize as Z
except Exception:
    Z = None

CSS = """
:root{--ink:#16181d;--soft:#5b6270;--faint:#8b929e;--line:#e3e6ec;--bg:#fff;--panel:#f7f8fa;
--pos:#1a7f5a;--neg:#b3261e;--warn:#8a5a00;--accent:#2f5fd0}
@media (prefers-color-scheme:dark){:root{--ink:#e8eaed;--soft:#a8afbb;--faint:#7d848f;
--line:#2a2f39;--bg:#0f1116;--panel:#161a21;--pos:#54c79a;--neg:#e5787d;--warn:#e0a33a;--accent:#8fb6ff}}
:root[data-theme="dark"]{--ink:#e8eaed;--soft:#a8afbb;--faint:#7d848f;--line:#2a2f39;--bg:#0f1116;
--panel:#161a21;--pos:#54c79a;--neg:#e5787d;--warn:#e0a33a;--accent:#8fb6ff}
:root[data-theme="light"]{--ink:#16181d;--soft:#5b6270;--faint:#8b929e;--line:#e3e6ec;--bg:#fff;
--panel:#f7f8fa;--pos:#1a7f5a;--neg:#b3261e;--warn:#8a5a00;--accent:#2f5fd0}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:920px;margin:0 auto;padding:40px 22px 72px}
h1{font-size:1.55rem;line-height:1.25;margin:0 0 6px;letter-spacing:-.01em}
h2{font-size:1.05rem;margin:38px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--line)}
h3{font-size:.93rem;margin:20px 0 6px;color:var(--soft);font-weight:600}
p{margin:8px 0}
.sub{color:var(--soft);font-size:.9rem;margin:0 0 4px}
.meta{color:var(--faint);font-size:.8rem;margin:0}
table{width:100%;border-collapse:collapse;margin:10px 0;font-size:.88rem}
th,td{text-align:right;padding:7px 9px;border-bottom:1px solid var(--line);
font-variant-numeric:tabular-nums}
th:first-child,td:first-child{text-align:left;font-variant-numeric:normal}
thead th{color:var(--soft);font-weight:600;font-size:.8rem;text-transform:uppercase;
letter-spacing:.04em;border-bottom:1px solid var(--line)}
tbody tr:hover{background:var(--panel)}
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
.pos{color:var(--pos)}.neg{color:var(--neg)}.warn{color:var(--warn)}
.box{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:12px 0}
.box.verdict{border-left:3px solid var(--accent)}
.k{color:var(--soft);font-size:.82rem;text-transform:uppercase;letter-spacing:.05em;margin:0 0 4px}
.v{font-size:1.02rem;font-weight:600}
ul{margin:6px 0;padding-left:20px}li{margin:3px 0}
figure{margin:16px 0}figure svg{width:100%;height:auto;display:block;
border:1px solid var(--line);border-radius:8px}
figcaption{color:var(--faint);font-size:.8rem;margin-top:6px}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.86em;
background:var(--panel);padding:1px 5px;border-radius:4px}
.foot{margin-top:44px;padding-top:14px;border-top:1px solid var(--line);
color:var(--faint);font-size:.78rem}
@media print{body{background:#fff;color:#000}.wrap{max-width:none}figure svg{border:1px solid #ccc}}
"""


def _f(v, d=2, sign=False):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "—"
    s = f"{v:+.{d}f}" if sign else f"{v:.{d}f}"
    return s


def _cls(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return "pos" if v > 0 else ("neg" if v < 0 else "")


def _esc(s):
    return html.escape(str(s))


class Report:
    def __init__(self, title, subtitle="", generated_by="quant-rigor", date=None):
        self.title, self.subtitle = title, subtitle
        self.generated_by, self.date = generated_by, date
        self.parts, self.series, self.order = [], {}, []
        self._verdict = None

    # ───────────────────────────────────────── content blocks
    def assumptions(self, data="", periods=None, resolution="", costs="",
                    not_modelled=None, extra=None):
        """State the frame before any number. A report without this section is not auditable."""
        rows = []
        if data:
            rows.append(("Data source / sample", data))
        if periods is not None:
            rows.append(("Independent periods (n)", str(periods)))
        if resolution:
            rows.append(("Resolution", resolution))
        if costs:
            rows.append(("Costs applied", costs))
        for k, v in (extra or {}).items():
            rows.append((k, str(v)))
        h = ["<h2>1 · Assumptions and data</h2>", "<div class='scroll'><table><tbody>"]
        for k, v in rows:
            h.append(f"<tr><td>{_esc(k)}</td><td style='text-align:left'>{_esc(v)}</td></tr>")
        h.append("</tbody></table></div>")
        if not_modelled:
            h.append("<h3>Not modelled — these costs and risks are absent from every figure below</h3><ul>")
            for x in not_modelled:
                h.append(f"<li>{_esc(x)}</li>")
            h.append("</ul>")
        self.parts.append("\n".join(h))
        return self

    def measurability(self, bar_ranges, rules, refinement=None):
        """rules: {name: characteristic distance}. Answers whether the data can measure the rule."""
        h = ["<h2>2 · Measurability</h2>",
             "<p>A rule is measurable only if its characteristic distance is large relative to the "
             "typical bar range. Ratios below 3 mean the result is an artefact of resolution rather "
             "than a property of the strategy.</p>",
             "<div class='scroll'><table><thead><tr><th>Rule</th><th>Distance</th>"
             "<th>Median bar range</th><th>R</th><th>Assessment</th></tr></thead><tbody>"]
        med = st.median(bar_ranges) if bar_ranges else float('nan')
        for name, dist in rules.items():
            m = S.measurability(bar_ranges, dist)
            cl = "" if m["ratio_R"] >= 3 else ("warn" if m["ratio_R"] >= 1 else "neg")
            h.append(f"<tr><td>{_esc(name)}</td><td>{_f(dist,1)}</td><td>{_f(med,1)}</td>"
                     f"<td class='{cl}'>{_f(m['ratio_R'])}</td>"
                     f"<td style='text-align:left' class='{cl}'>{_esc(m['verdict'])}</td></tr>")
        h.append("</tbody></table></div>")
        if refinement:
            h.append(f"<p>Path refinement factor measured at {_f(refinement.get('median_ratio'))} "
                     f"against a random-walk null of {_f(refinement.get('random_walk_null_sqrt_n'))}. "
                     f"A ratio near the null means measured path length does not converge as "
                     f"resolution improves, so path-dependent rules cannot be evaluated with this data.</p>")
        self.parts.append("\n".join(h))
        return self

    def add_series(self, name, values):
        self.series[name] = list(values)
        if name not in self.order:
            self.order.append(name)
        return self

    def breakdown(self, dated, title="3b · Year by year, month by month", split=None):
        """Temporal decomposition of one or more dated P&L series.

        ``dated`` maps a series name to a list of ``(date_string, value)`` pairs, dates as
        ``YYYY-MM-DD``. Emits, per series: a year table (n, total, mean, t, win rate, max
        drawdown of the year's own equity) and a year x month matrix of totals.

        The seasonality row is printed with its best-of-12 search band, because picking the
        best calendar month out of twelve is a 12-variant search: a month total inside that
        band carries no information and must not be read as a seasonal effect.
        """
        import collections
        import statistics as _st
        h = [f"<h2>{_esc(title)}</h2>"]
        for name, ser in dated.items():
            vals = [v for _, v in ser]
            by_y, by_ym, by_m = (collections.defaultdict(list) for _ in range(3))
            for d, v in ser:
                by_y[d[:4]].append(v)
                by_ym[(d[:4], int(d[5:7]))].append(v)
                by_m[int(d[5:7])].append(v)
            h.append(f"<h3>{_esc(name)} — by year</h3><div class='scroll'><table><thead><tr>"
                     "<th>Year</th><th>n</th><th>Total</th><th>Mean</th><th>t</th>"
                     "<th>Win %</th><th>Worst period</th><th>Max drawdown</th>"
                     + ("<th>Side</th>" if split else "") + "</tr></thead><tbody>")
            for y in sorted(by_y):
                v = by_y[y]
                d = S.describe(v)
                dd = S.drawdown_from_pnl(v)["max_dd"]
                wr = 100.0 * sum(1 for x in v if x > 0) / len(v)
                side = ("<td>" + ("in" if y < split[:4] else "out") + "</td>") if split else ""
                h.append(f"<tr><td>{y}</td><td>{len(v)}</td>"
                         f"<td class='{_cls(d['total'])}'>{_f(d['total'],0,True)}</td>"
                         f"<td class='{_cls(d['mean'])}'>{_f(d['mean'],1,True)}</td>"
                         f"<td>{_f(d['t'],2,True)}</td><td>{_f(wr,0)}%</td>"
                         f"<td class='neg'>{_f(min(v),0,True)}</td>"
                         f"<td class='neg'>{_f(dd,0,True)}</td>{side}</tr>")
            dtot = S.describe(vals)
            h.append(f"<tr class='tot'><td><b>all</b></td><td>{len(vals)}</td>"
                     f"<td class='{_cls(dtot['total'])}'><b>{_f(dtot['total'],0,True)}</b></td>"
                     f"<td class='{_cls(dtot['mean'])}'>{_f(dtot['mean'],1,True)}</td>"
                     f"<td>{_f(dtot['t'],2,True)}</td>"
                     f"<td>{_f(100.0*sum(1 for x in vals if x>0)/len(vals),0)}%</td>"
                     f"<td class='neg'>{_f(min(vals),0,True)}</td>"
                     f"<td class='neg'>{_f(S.drawdown_from_pnl(vals)['max_dd'],0,True)}</td>"
                     + ("<td></td>" if split else "") + "</tr>")
            h.append("</tbody></table></div>")
            M = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
            h.append(f"<h3>{_esc(name)} — month by month (total per month)</h3>"
                     "<div class='scroll'><table><thead><tr><th>Year</th>"
                     + "".join(f"<th>{m}</th>" for m in M) + "<th>Year</th></tr></thead><tbody>")
            for y in sorted(by_y):
                row, tot = [], 0.0
                for m in range(1, 13):
                    v = by_ym.get((y, m))
                    if v:
                        tot += sum(v)
                        row.append(f"<td class='{_cls(sum(v))}'>{_f(sum(v),0,True)}</td>")
                    else:
                        row.append("<td class='meta'>·</td>")
                h.append(f"<tr><td>{y}</td>" + "".join(row)
                         + f"<td class='{_cls(tot)}'><b>{_f(tot,0,True)}</b></td></tr>")
            n_months = len(by_ym)
            sd = _st.pstdev(vals) if len(vals) > 1 else 0.0
            n_avg = _st.mean([len(v) for v in by_m.values()]) if by_m else 0
            band = S.headroom(sd, n_avg, 12)
            cells = []
            for m in range(1, 13):
                v = by_m.get(m)
                cells.append(f"<td class='{_cls(sum(v))}'>{_f(sum(v),0,True)}</td>" if v
                             else "<td class='meta'>·</td>")
            h.append("<tr class='tot'><td><b>all yrs</b></td>" + "".join(cells)
                     + "<td></td></tr></tbody></table></div>")
            n_pos = sum(1 for k, v in by_ym.items() if sum(v) > 0)
            best = max(((sum(v), m) for m, v in by_m.items()), default=(0, 0))
            inside = abs(best[0]) < band
            h.append(f"<p class='meta'>{n_pos} of {n_months} calendar months positive "
                     f"({_f(100.0*n_pos/n_months,0)}%). Best calendar month: "
                     f"{M[best[1]-1]} at {_f(best[0],0,True)}; best-of-12 search band on "
                     f"{_f(n_avg,0)} periods is ±{_f(band,0)}, so this month is "
                     f"<b>{'inside the noise — not a seasonal effect' if inside else 'outside the search band'}</b>. "
                     f"Year and month totals are descriptive: the number of periods per cell is "
                     f"small and each cell carries its own sampling error.</p>")
        self.parts.append("\n".join(h))
        return self

    def results(self, title="3 · Results"):
        """Full metric set for every series. Same columns for all, no selective reporting."""
        if not self.series:
            return self
        cols = [("Total", lambda d: _f(d["total"], 1, True)),
                ("Mean/period", lambda d: _f(d["mean"], 2, True)),
                ("SE", lambda d: _f(d["se"], 2)),
                ("t", lambda d: _f(d["t"], 2, True)),
                ("p", lambda d: _f(d["p"], 3)),
                ("Win %", lambda d: _f(d["win_rate"], 0)),
                ("Worst", lambda d: _f(d["worst"], 1, True)),
                ("Skew", lambda d: _f(d["skew"], 2, True)),
                ("n", lambda d: str(d["n"]))]
        h = [f"<h2>{_esc(title)}</h2>", "<div class='scroll'><table><thead><tr><th>Series</th>"]
        h += [f"<th>{c}</th>" for c, _ in cols]
        h.append("</tr></thead><tbody>")
        for name in self.order:
            d = S.describe(self.series[name], name)
            h.append(f"<tr><td>{_esc(name)}</td>")
            for _, fn in cols:
                h.append(f"<td>{fn(d)}</td>")
            h.append("</tr>")
        h.append("</tbody></table></div>")
        h.append("<h3>Drawdown and dispersion</h3><div class='scroll'><table><thead><tr>"
                 "<th>Series</th><th>Max drawdown</th><th>Periods under water</th>"
                 "<th>Longest run</th><th>Empirical VaR 95%</th><th>CVaR</th>"
                 "<th>Top-3 of gains</th><th>Top-3 of losses</th></tr></thead><tbody>")
        for name in self.order:
            v = self.series[name]
            dd = S.drawdown(S.equity_from_returns(v))
            vc = S.var_cvar(v)
            cc = S.concentration(v)
            h.append(f"<tr><td>{_esc(name)}</td><td class='neg'>{_f(dd['max_dd'],1)}</td>"
                     f"<td>{dd['periods_under_water']}/{len(v)}</td>"
                     f"<td>{dd['longest_under_water']}</td>"
                     f"<td>{_f(vc['VaR_empirical'],1,True)}</td>"
                     f"<td>{_f(vc['CVaR_expected_shortfall'],1,True)}</td>"
                     f"<td>{_f(cc['gain_share_pct'],0)}%</td>"
                     f"<td>{_f(cc['loss_share_pct'],0)}%</td></tr>")
        h.append("</tbody></table></div>")
        h.append("<p class='meta'>t and p are two-sided against a zero mean. They describe the "
                 "sample, not the future, and assume the periods are independent.</p>")
        self.parts.append("\n".join(h))
        return self

    def compare(self, base_name, variant_name, n_variants=1, title="4 · Difference and uncertainty"):
        a, b = self.series[base_name], self.series[variant_name]
        pr = S.paired(a, b, (base_name, variant_name))
        bsa, bsb = S.bootstrap(a, reps=4000), S.bootstrap(b, reps=4000)
        jk = S.jackknife(pr["diffs"])
        hr = S.headroom(S.describe(a)["sd"], len(a), n_variants)
        gap = sum(b) - sum(a)
        inside = abs(gap) < hr if hr else False
        h = [f"<h2>{_esc(title)}</h2>",
             "<h3>Paired difference (variant − baseline, same periods)</h3>",
             "<div class='scroll'><table><tbody>",
             f"<tr><td>Mean difference per period</td><td class='{_cls(pr['mean'])}'>"
             f"{_f(pr['mean'],2,True)} ± {_f(pr['se'],2)} (SE)</td></tr>",
             f"<tr><td>t · p</td><td>{_f(pr['t'],2,True)} · {_f(pr['p'],3)}</td></tr>",
             f"<tr><td>Periods improved / worsened / unchanged</td>"
             f"<td>{pr['improved']} / {pr['worsened']} / {pr['unchanged']}</td></tr>",
             f"<tr><td>Jackknife: sign survives</td><td>{_esc(jk['sign_survives'])} "
             f"(worst deletion {_f(jk['worst_value'],0,True)})</td></tr>",
             "</tbody></table></div>",
             "<h3>Bootstrap band on the total (4,000 resamples)</h3>",
             "<div class='scroll'><table><thead><tr><th>Series</th><th>5th</th><th>Median</th>"
             "<th>95th</th><th>Resamples with opposite sign</th></tr></thead><tbody>"]
        for nm, bs in ((base_name, bsa), (variant_name, bsb)):
            h.append(f"<tr><td>{_esc(nm)}</td><td>{_f(bs['p05'],0,True)}</td>"
                     f"<td>{_f(bs['p50'],0,True)}</td><td>{_f(bs['p95'],0,True)}</td>"
                     f"<td>{_f(bs.get('share_opposite_sign_pct', bs.get('share_opposite_sign')),0)}%</td></tr>")
        h.append("</tbody></table></div>")
        if n_variants > 1:
            h.append("<h3>Search correction</h3>"
                     f"<p>{n_variants} configurations were evaluated. The expected best-of-N "
                     f"deviation under no true effect is ±{_f(hr,0)} "
                     f"(<code>σ·√n·√(2 ln N)</code>). The observed gap is {_f(gap,0,True)}, "
                     f"which is <b>{'inside' if inside else 'outside'}</b> that band.</p>")
            if inside:
                h.append("<p>A gap inside the band is not distinguishable from the effect of "
                         "searching. It supports neither adoption nor rejection on magnitude alone.</p>")
        self.parts.append("\n".join(h))
        return self

    def contraction_check(self, base_name, variant_name):
        """Explicit algebra check: does the variant only 'win' where the baseline lost?"""
        a, b = self.series[base_name], self.series[variant_name]
        idx = [i for i, v in enumerate(a) if v < 0]
        if not idx:
            return self
        an, bn = sum(a[i] for i in idx), sum(b[i] for i in idx)
        ap, bp = sum(v for v in a if v >= 0), sum(b[i] for i, v in enumerate(a) if v >= 0)
        h = ["<h3>Contraction check</h3>",
             "<div class='scroll'><table><thead><tr><th>Subset</th>"
             f"<th>{_esc(base_name)}</th><th>{_esc(variant_name)}</th><th>Difference</th>"
             "</tr></thead><tbody>",
             f"<tr><td>Periods where baseline lost ({len(idx)})</td><td>{_f(an,0,True)}</td>"
             f"<td>{_f(bn,0,True)}</td><td class='{_cls(bn-an)}'>{_f(bn-an,0,True)}</td></tr>",
             f"<tr><td>Periods where baseline gained ({len(a)-len(idx)})</td><td>{_f(ap,0,True)}</td>"
             f"<td>{_f(bp,0,True)}</td><td class='{_cls(bp-ap)}'>{_f(bp-ap,0,True)}</td></tr>",
             "</tbody></table></div>"]
        if bn > an and bp < ap:
            h.append("<p>The variant improves the losing subset and reduces the gaining subset. "
                     "This is the signature of a contraction: a rule that caps both tails "
                     "necessarily shows a higher result on any window where the baseline is negative, "
                     "independently of whether it has merit. Comparisons restricted to such a "
                     "window carry no information.</p>")
        self.parts.append("\n".join(h))
        return self

    def overfitting(self, series_name=None, perf_matrix=None, n_variants=1,
                    sd_of_trial_sharpes=None, title="5 · Overfitting diagnostics"):
        if V is None:
            return self
        h = [f"<h2>{_esc(title)}</h2>", "<div class='scroll'><table><tbody>"]
        if series_name:
            v = self.series[series_name]
            d = S.describe(v)
            hac = V.hac_t(v)
            h.append(f"<tr><td>t, naive vs HAC (Newey-West, {hac['lags']} lags)</td>"
                     f"<td>{_f(hac['t_naive'],2,True)} vs {_f(hac['t_hac'],2,True)} "
                     f"(SE ×{_f(hac['se_inflation'])})</td></tr>")
            h.append(f"<tr><td>Effective sample size</td>"
                     f"<td>{_f(V.effective_sample_size(v),1)} of {len(v)}</td></tr>")
            if d["sd"]:
                sr = d["mean"] / d["sd"]
                mt = V.min_track_record_length(sr, d["skew"], 3.0)
                h.append(f"<tr><td>Per-period Sharpe</td><td>{_f(sr,3,True)}</td></tr>")
                h.append(f"<tr><td>Minimum track record length (95%)</td>"
                         f"<td>{_f(mt,0)} periods, {len(v)} available "
                         f"{'(sufficient)' if len(v) >= mt else '(insufficient)'}</td></tr>")
                if n_variants > 1:
                    sdt = sd_of_trial_sharpes if sd_of_trial_sharpes else abs(sr) * 0.5
                    ds = V.deflated_sharpe(sr, len(v), n_variants, sdt, d["skew"], 3.0)
                    h.append(f"<tr><td>Deflated Sharpe (N={n_variants})</td>"
                             f"<td>hurdle SR₀ {_f(ds.get('sr0'),3,True)} · "
                             f"DSR {_f(ds.get('dsr'),3)} — {_esc(ds.get('verdict',''))}</td></tr>")
        if perf_matrix:
            pb = V.pbo(perf_matrix, n_splits=6, max_combos=60)
            h.append(f"<tr><td>Probability of backtest overfitting (PBO)</td>"
                     f"<td>{_f(pb.get('pbo'),2)} over {pb.get('n_partitions','—')} partitions — "
                     f"{_esc(pb.get('verdict',''))}</td></tr>")
        h.append("</tbody></table></div>")
        h.append("<p class='meta'>PBO measures the selection procedure, not a single "
                 "configuration. DSR corrects an observed Sharpe for the number of trials, "
                 "non-normality and sample length.</p>")
        self.parts.append("\n".join(h))
        return self

    def sizing(self, mean_per_period, sd_per_period, dd_tolerated=0.20, prob=0.10,
               title="6 · Sizing constraints"):
        if Z is None:
            return self
        s = Z.sizing_summary(mean_per_period, sd_per_period, dd_tolerated, prob)
        h = [f"<h2>{_esc(title)}</h2>", "<div class='scroll'><table><tbody>",
             f"<tr><td>Kelly ceiling f* (μ/σ²)</td><td>{_f(s['kelly_full'],3)}</td></tr>",
             f"<tr><td>Quarter Kelly</td><td>{_f(s['kelly_fractional'],3)}</td></tr>",
             f"<tr><td>Drawdown-constrained f (≤{dd_tolerated:.0%} at {prob:.0%} probability)</td>"
             f"<td>{_f(s['dd_constrained'],3)}</td></tr>",
             f"<tr><td>Binding constraint</td><td>{_esc(s['binding_constraint'])}</td></tr>",
             f"<tr><td>Resulting fraction</td><td>{_f(s['recommended_f'],3)}</td></tr>",
             f"<tr><td>Expected log-growth at that fraction</td>"
             f"<td>{_f(s['growth_at_recommended'],5,True)} per period</td></tr>",
             "</tbody></table></div>",
             "<p class='meta'>Kelly is an upper bound derived from estimated moments, not a target; "
             "above it, additional size lowers median wealth while raising the mean. The drawdown "
             "constraint uses a barrier approximation and assumes the estimated moments hold.</p>"]
        self.parts.append("\n".join(h))
        return self

    def charts(self, param_grid=None, headroom=None, title="7 · Charts"):
        if not self.series:
            return self
        h = [f"<h2>{_esc(title)}</h2>"]
        figs = [(C.equity_chart({k: self.series[k] for k in self.order},
                                title="Cumulative result by period"),
                 "Cumulative result. A step indicates the outcome depends on a single period; "
                 "a slope indicates a persistent effect."),
               (C.drawdown_chart({k: self.series[k] for k in self.order}),
                "Drawdown depth and duration. Equal maxima with different durations are not "
                "equivalent outcomes."),
               (C.dist_chart(self.series[self.order[0]],
                             title=f"Distribution — {self.order[0]}"),
                "Per-period distribution with a normal fit overlaid. Departures indicate that "
                "Sharpe- and VaR-based statements understate tail risk.")]
        if len(self.order) >= 2:
            pr = S.paired(self.series[self.order[0]], self.series[self.order[1]])
            figs.insert(2, (C.contribution_chart(pr["diffs"],
                            title=f"Per-period difference: {self.order[1]} − {self.order[0]}"),
                            "Per-period paired difference, sorted. Reveals whether the difference "
                            "is distributed or concentrated in a few periods."))
        if param_grid:
            params, metric = param_grid
            figs.append((C.surface_chart(params, metric, headroom=headroom,
                                         title="Metric by parameter value"),
                         "Metric against parameter. An isolated peak indicates a fit to noise; a "
                         "broad plateau indicates robustness. The shaded band, when present, is the "
                         "best-of-N noise range."))
        for svg, cap in figs:
            body = svg.split("?>", 1)[-1]
            h.append(f"<figure>{body}<figcaption>{_esc(cap)}</figcaption></figure>")
        self.parts.append("\n".join(h))
        return self

    def verdict(self, status, reason, path_to_decide="", does_not_support=None,
                default_action="Change nothing."):
        """status: 'decidable-positive' | 'decidable-negative' | 'not decidable'."""
        self._verdict = dict(status=status, reason=reason, path=path_to_decide,
                             limits=does_not_support or [], action=default_action)
        return self

    def _verdict_html(self):
        v = self._verdict
        if not v:
            return ""
        h = ["<h2>8 · Verdict</h2>", "<div class='box verdict'>",
             f"<p class='k'>Status</p><p class='v'>{_esc(v['status'])}</p>",
             f"<p>{_esc(v['reason'])}</p>"]
        if v["path"]:
            h.append(f"<p class='k' style='margin-top:12px'>Path to a decision</p>"
                     f"<p>{_esc(v['path'])}</p>")
        h.append(f"<p class='k' style='margin-top:12px'>Default action</p>"
                 f"<p>{_esc(v['action'])} The burden of proof falls on the modification, "
                 f"not on the current configuration.</p></div>")
        if v["limits"]:
            h.append("<h3>What this evidence does not support</h3><ul>")
            for x in v["limits"]:
                h.append(f"<li>{_esc(x)}</li>")
            h.append("</ul>")
        return "\n".join(h)

    # ───────────────────────────────────────── output
    def to_html(self):
        head = [f"<div class='wrap'><h1>{_esc(self.title)}</h1>"]
        if self.subtitle:
            head.append(f"<p class='sub'>{_esc(self.subtitle)}</p>")
        stamp = f" · {_esc(self.date)}" if self.date else ""
        head.append(f"<p class='meta'>Prepared with {_esc(self.generated_by)}{stamp}. "
                    f"All figures computed from the stated sample; none are forecasts.</p>")
        foot = ("<div class='foot'><p>Research and educational material. Not investment advice, "
                "not a recommendation, and not an offer to buy or sell any instrument. Backtested "
                "and simulated results are hypothetical, do not reflect actual trading, and are not "
                "indicative of future results. Trading involves risk of substantial loss.</p></div>")
        return ("<!doctype html><html><head><meta charset='utf-8'>"
                "<meta name='viewport' content='width=device-width,initial-scale=1'>"
                f"<title>{_esc(self.title)}</title><style>{CSS}</style></head><body>"
                + "\n".join(head) + "\n" + "\n".join(self.parts) + "\n"
                + self._verdict_html() + foot + "</div></body></html>")

    def to_markdown(self):
        """Plain-text sibling for commit messages, issues and decision records."""
        L = [f"# {self.title}"]
        if self.subtitle:
            L.append(f"_{self.subtitle}_")
        for name in self.order:
            d = S.describe(self.series[name], name)
            L.append(f"\n**{name}** — total {_f(d['total'],1,True)} · "
                     f"mean {_f(d['mean'],2,True)} ± {_f(d['se'],2)} · t {_f(d['t'],2,True)} · "
                     f"p {_f(d['p'],3)} · n {d['n']} · worst {_f(d['worst'],1,True)}")
        if self._verdict:
            v = self._verdict
            L.append(f"\n**Verdict: {v['status']}** — {v['reason']}")
            if v["path"]:
                L.append(f"Path to decide: {v['path']}")
            L.append(f"Default action: {v['action']}")
        return "\n".join(L)

    def write(self, path):
        with open(path, "w") as f:
            f.write(self.to_html())
        return path
