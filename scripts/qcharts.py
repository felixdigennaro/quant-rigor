#!/usr/bin/env python3
"""qcharts — the five required strategy charts as standalone SVG. STDLIB ONLY.

No matplotlib needed. Output is a self-contained .svg (openable in a browser, embeddable in an
artifact, publishable). Dark-theme by default, light-theme safe via currentColor where possible.

    import sys; sys.path.insert(0,'<skill>/scripts')
    from qcharts import equity_chart, drawdown_chart, contribution_chart, surface_chart, dist_chart, write

    write("equity.svg", equity_chart({"baseline": base, "TP+100": variant}))

CLI demo:
    python3 qcharts.py --a 1 -2 3 --b 2 -1 4 --out /tmp/chart
"""
import argparse
import math
import statistics as st

W, H = 900, 380
PAD = dict(l=64, r=18, t=34, b=42)
INK = "#e8eaed"
SOFT = "#9aa0a6"
GRID = "#2a2f3a"
BG = "#0f1116"
SERIES = ["#8fb6ff", "#5fdfe6", "#f0b357", "#e5787d", "#a78bfa"]


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _hdr(title, sub=""):
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif">',
         f'<rect width="{W}" height="{H}" fill="{BG}"/>',
         f'<text x="{PAD["l"]}" y="20" fill="{INK}" font-size="14" font-weight="600">{_esc(title)}</text>']
    if sub:
        s.append(f'<text x="{PAD["l"]}" y="{H-12}" fill="{SOFT}" font-size="10.5">{_esc(sub)}</text>')
    return s


def _axes(xmin, xmax, ymin, ymax, xlab="", ylab="", yfmt="{:.0f}"):
    """Returns (svg_parts, X(), Y()) mapping data → pixels."""
    x0, x1 = PAD["l"], W - PAD["r"]
    y0, y1 = H - PAD["b"], PAD["t"]
    if xmax == xmin:
        xmax = xmin + 1
    if ymax == ymin:
        ymax = ymin + 1
    X = lambda v: x0 + (v - xmin) / (xmax - xmin) * (x1 - x0)
    Y = lambda v: y0 + (v - ymin) / (ymax - ymin) * (y1 - y0)
    p = []
    for i in range(5):
        v = ymin + i * (ymax - ymin) / 4
        y = Y(v)
        p.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x1}" y2="{y:.1f}" stroke="{GRID}" stroke-width="1"/>')
        p.append(f'<text x="{x0-8}" y="{y+3.5:.1f}" fill="{SOFT}" font-size="10" text-anchor="end">'
                 f'{yfmt.format(v)}</text>')
    if ymin < 0 < ymax:
        p.append(f'<line x1="{x0}" y1="{Y(0):.1f}" x2="{x1}" y2="{Y(0):.1f}" stroke="{SOFT}" '
                 f'stroke-width="1" stroke-dasharray="3 3" opacity=".6"/>')
    if xlab:
        p.append(f'<text x="{(x0+x1)/2}" y="{H-24}" fill="{SOFT}" font-size="10" text-anchor="middle">{_esc(xlab)}</text>')
    if ylab:
        p.append(f'<text x="14" y="{(y0+y1)/2}" fill="{SOFT}" font-size="10" text-anchor="middle" '
                 f'transform="rotate(-90 14 {(y0+y1)/2})">{_esc(ylab)}</text>')
    return p, X, Y


def _legend(names):
    out = []
    x = PAD["l"]
    for i, n in enumerate(names):
        c = SERIES[i % len(SERIES)]
        out.append(f'<rect x="{x}" y="{PAD["t"]-16}" width="10" height="3" rx="1.5" fill="{c}"/>')
        out.append(f'<text x="{x+15}" y="{PAD["t"]-11}" fill="{SOFT}" font-size="10.5">{_esc(n)}</text>')
        x += 22 + 7 * len(str(n))
    return out


def _cum(x):
    out, a = [], 0.0
    for v in x:
        a += v
        out.append(a)
    return out

# ------------------------------------------------------------------ 1. equity


def equity_chart(series: dict, title="Equity curve (cumulative)", sub="", labels=None):
    """series: {name: per-period list}. Exposes: slope (real edge) vs step (one lucky period)."""
    cums = {k: _cum(v) for k, v in series.items()}
    n = max(len(v) for v in cums.values())
    lo = min(0, min(min(v) for v in cums.values()))
    hi = max(max(v) for v in cums.values())
    pad = (hi - lo) * .08 or 1
    s = _hdr(title, sub or f"n={n} periods · cumulative units")
    ax, X, Y = _axes(0, n - 1, lo - pad, hi + pad, "period", "cumulative")
    s += ax + _legend(list(cums.keys()))
    for i, (name, cum) in enumerate(cums.items()):
        c = SERIES[i % len(SERIES)]
        pts = " ".join(f"{X(j):.1f},{Y(v):.1f}" for j, v in enumerate(cum))
        s.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="2.4" '
                 f'stroke-linejoin="round" stroke-linecap="round"/>')
        s.append(f'<circle cx="{X(len(cum)-1):.1f}" cy="{Y(cum[-1]):.1f}" r="3.5" fill="{c}"/>')
        s.append(f'<text x="{X(len(cum)-1)-6:.1f}" y="{Y(cum[-1])-9:.1f}" fill="{c}" font-size="11" '
                 f'font-weight="600" text-anchor="end">{cum[-1]:+.0f}</text>')
    s.append("</svg>")
    return "\n".join(s)

# ------------------------------------------------------------------ 2. drawdown


def drawdown_chart(series: dict, title="Underwater curve (drawdown)", sub=""):
    """Depth AND duration. Two strategies with equal max-DD differ if one stays under water longer."""
    dds = {}
    for k, v in series.items():
        cum = _cum(v)
        peak, d = cum[0], []
        for c in cum:
            peak = max(peak, c)
            d.append(c - peak)
        dds[k] = d
    n = max(len(v) for v in dds.values())
    lo = min(min(v) for v in dds.values())
    s = _hdr(title, sub or "0 = at high-water mark; deeper and longer are both bad")
    ax, X, Y = _axes(0, n - 1, lo * 1.12 or -1, 0, "period", "drawdown")
    s += ax + _legend(list(dds.keys()))
    for i, (name, d) in enumerate(dds.items()):
        c = SERIES[i % len(SERIES)]
        pts = " ".join(f"{X(j):.1f},{Y(v):.1f}" for j, v in enumerate(d))
        area = f"{X(0):.1f},{Y(0):.1f} " + pts + f" {X(len(d)-1):.1f},{Y(0):.1f}"
        s.append(f'<polygon points="{area}" fill="{c}" opacity=".13"/>')
        s.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="2"/>')
        under = sum(1 for v in d if v < 0)
        s.append(f'<text x="{W-PAD["r"]}" y="{PAD["t"]-11+i*13}" fill="{c}" font-size="10.5" '
                 f'text-anchor="end">{_esc(name)}: max {min(d):+.0f} · under water {under}/{len(d)}</text>')
    s.append("</svg>")
    return "\n".join(s)

# ------------------------------------------------------------------ 3. contribution


def contribution_chart(values, labels=None, title="Per-period contribution (sorted)", sub=""):
    """THE cheapest, most informative chart. Exposes concentration: does it live in 3 periods?"""
    idx = sorted(range(len(values)), key=lambda i: values[i], reverse=True)
    vals = [values[i] for i in idx]
    labs = [(labels[i] if labels else str(i)) for i in idx]
    tot = sum(values)
    top3 = sum(sorted(values, key=abs, reverse=True)[:3])
    share = 100 * top3 / tot if tot else float("nan")
    s = _hdr(title, sub or f"total {tot:+.0f} · top-3 periods = {share:.0f}% of total · "
                          f"{sum(1 for v in values if v>0)} up / {sum(1 for v in values if v<0)} down")
    lo, hi = min(0, min(vals)), max(0, max(vals))
    pad = (hi - lo) * .1 or 1
    ax, X, Y = _axes(-0.6, len(vals) - 0.4, lo - pad, hi + pad, "", "contribution")
    s += ax
    bw = max(4.0, (W - PAD["l"] - PAD["r"]) / max(1, len(vals)) * 0.66)
    for j, v in enumerate(vals):
        c = SERIES[1] if v >= 0 else SERIES[3]
        y = Y(max(v, 0))
        h = abs(Y(v) - Y(0))
        s.append(f'<rect x="{X(j)-bw/2:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{max(h,0.8):.1f}" '
                 f'rx="2" fill="{c}" opacity=".9"/>')
        if len(vals) <= 24:
            s.append(f'<text x="{X(j):.1f}" y="{H-PAD["b"]+13:.1f}" fill="{SOFT}" font-size="8.5" '
                     f'text-anchor="middle" transform="rotate(35 {X(j):.1f} {H-PAD["b"]+13:.1f})">'
                     f'{_esc(labs[j])}</text>')
    s.append("</svg>")
    return "\n".join(s)

# ------------------------------------------------------------------ 4. parameter surface


def surface_chart(params, metric, headroom=None, title="Parameter surface", sub="", xlab="parameter"):
    """Sharp isolated peak = overfit. Broad plateau = robust (adopt the MIDDLE).
    Pass `headroom` (best-of-N noise band) to show what variation is meaningless."""
    lo, hi = min(min(metric), 0), max(metric)
    pad = (hi - lo) * .12 or 1
    best = max(range(len(metric)), key=lambda i: metric[i])
    s = _hdr(title, sub or ("peak at " + str(params[best]) +
                            (f" · noise headroom ±{headroom:.0f} — variation inside it is luck" if headroom else "")))
    ax, X, Y = _axes(min(params), max(params), lo - pad, hi + pad, xlab, "metric")
    s += ax
    if headroom:
        top = metric[best]
        y1, y2 = Y(top), Y(max(lo - pad, top - 2 * headroom))
        s.append(f'<rect x="{PAD["l"]}" y="{min(y1,y2):.1f}" width="{W-PAD["l"]-PAD["r"]}" '
                 f'height="{abs(y2-y1):.1f}" fill="{SOFT}" opacity=".10"/>')
        s.append(f'<text x="{W-PAD["r"]-6}" y="{min(y1,y2)+13:.1f}" fill="{SOFT}" font-size="9.5" '
                 f'text-anchor="end">within noise of the peak</text>')
    pts = " ".join(f"{X(p):.1f},{Y(m):.1f}" for p, m in zip(params, metric))
    s.append(f'<polyline points="{pts}" fill="none" stroke="{SERIES[0]}" stroke-width="2.4"/>')
    for p, m in zip(params, metric):
        s.append(f'<circle cx="{X(p):.1f}" cy="{Y(m):.1f}" r="3" fill="{SERIES[0]}"/>')
    s.append(f'<circle cx="{X(params[best]):.1f}" cy="{Y(metric[best]):.1f}" r="5.5" fill="none" '
             f'stroke="{SERIES[2]}" stroke-width="2"/>')
    s.append("</svg>")
    return "\n".join(s)

# ------------------------------------------------------------------ 5. distribution


def dist_chart(values, bins=14, title="Return distribution", sub=""):
    """Where the tails actually are. Negative skew + fat left tail = short gamma (see 04-finance §27)."""
    n = len(values)
    m, sd = st.mean(values), (st.stdev(values) if n > 1 else 1)
    lo, hi = min(values), max(values)
    if hi <= lo:                      # degenerate (constant series): build a symmetric window
        pad = abs(lo) * 0.05 or 1.0   # so bins and the axis agree instead of running off-canvas
        lo, hi = lo - pad, hi + pad
    w = (hi - lo) / bins
    counts = [0] * bins
    for v in values:
        counts[min(bins - 1, int((v - lo) / w))] += 1
    sk = 0.0
    if n > 2 and sd:
        sk = (n / ((n - 1) * (n - 2))) * sum(((v - m) / sd) ** 3 for v in values)
    srt = sorted(values)
    var95 = srt[max(0, int(0.05 * n))]
    s = _hdr(title, sub or f"n={n} · mean {m:+.1f} · sd {sd:.1f} · skew {sk:+.2f} · "
                          f"empirical VaR95 {var95:+.1f} · worst {srt[0]:+.1f}")
    ax, X, Y = _axes(lo - w, hi + w, 0, max(counts) * 1.15 or 1, "outcome", "count")
    s += ax
    for i, c in enumerate(counts):
        x = lo + i * w
        bw = (X(x + w) - X(x)) * .88
        col = SERIES[3] if x + w / 2 < 0 else SERIES[1]
        s.append(f'<rect x="{X(x):.1f}" y="{Y(c):.1f}" width="{bw:.1f}" '
                 f'height="{abs(Y(0)-Y(c)):.1f}" fill="{col}" opacity=".85" rx="1.5"/>')
    # normal overlay (to show what Sharpe/VaR assume)
    if sd:
        pts = []
        for k in range(80):
            xv = lo - w + k * (hi - lo + 2 * w) / 79
            d = math.exp(-((xv - m) ** 2) / (2 * sd ** 2)) / (sd * math.sqrt(2 * math.pi))
            pts.append(f"{X(xv):.1f},{Y(d*n*w):.1f}")
        s.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{SERIES[2]}" '
                 f'stroke-width="1.6" stroke-dasharray="4 3" opacity=".85"/>')
        s.append(f'<text x="{W-PAD["r"]}" y="{PAD["t"]-11}" fill="{SERIES[2]}" font-size="10" '
                 f'text-anchor="end">- - normal fit (what Sharpe/VaR assume)</text>')
    s.append(f'<line x1="{X(var95):.1f}" y1="{PAD["t"]}" x2="{X(var95):.1f}" y2="{H-PAD["b"]}" '
             f'stroke="{SERIES[3]}" stroke-width="1.4" stroke-dasharray="2 3"/>')
    s.append("</svg>")
    return "\n".join(s)

# ------------------------------------------------------------------ helpers


def write(path, svg):
    with open(path, "w") as f:
        f.write(svg)
    return path


def sparkline(values, width=48):
    """ASCII sparkline for terminal output — no file needed."""
    blocks = "▁▂▃▄▅▆▇█"
    if not values:
        return ""
    step = max(1, len(values) // width)
    v = [values[i] for i in range(0, len(values), step)]
    lo, hi = min(v), max(v)
    rng = (hi - lo) or 1
    return "".join(blocks[min(7, int((x - lo) / rng * 7.999))] for x in v)


def _main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", nargs="+", type=float, required=True)
    ap.add_argument("--b", nargs="+", type=float)
    ap.add_argument("--out", default="/tmp/qchart")
    x = ap.parse_args()
    series = {"baseline": x.a}
    if x.b:
        series["variant"] = x.b
    write(f"{x.out}-equity.svg", equity_chart(series))
    write(f"{x.out}-dd.svg", drawdown_chart(series))
    write(f"{x.out}-contrib.svg", contribution_chart(x.a))
    write(f"{x.out}-dist.svg", dist_chart(x.a))
    print(f"wrote {x.out}-{{equity,dd,contrib,dist}}.svg")
    print("baseline sparkline:", sparkline(_cum(x.a)))


if __name__ == "__main__":
    _main()
