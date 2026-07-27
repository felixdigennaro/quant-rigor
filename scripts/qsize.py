#!/usr/bin/env python3
"""qsize — sizing, capacity and execution reality. STDLIB ONLY (no numpy/scipy).

An edge can be real and still lose money if it is sized or executed badly. These are the
constraints that decide how much of an edge reaches the account.

Reference notes: references/06-sizing-execution.md (formulas 38-44).
"""
import math
import statistics as st

# ---------------------------------------------------------------- 38. Kelly


def kelly_continuous(mean_per_period, sd_per_period):
    """f* = μ/σ² — the growth-optimal fraction. A CEILING, not a target.
    Above f*, extra size lowers MEDIAN wealth while raising the mean."""
    if sd_per_period == 0:
        return float('inf')
    return mean_per_period / (sd_per_period ** 2)


def kelly_discrete(p, win, loss):
    """f* = p/loss - (1-p)/win for a bet paying `win` on success, losing `loss` on failure
    (both positive, expressed as fractions of the stake)."""
    if win <= 0 or loss <= 0:
        return float('nan')
    return p / loss - (1 - p) / win


def fractional_kelly(f_star, fraction=0.25):
    """Quarter-to-half Kelly is the practical range: σ is estimated, and the drawdown of full
    Kelly is intolerable (expected max DD ≈ f*, i.e. ~100% at full Kelly)."""
    return f_star * fraction


def growth_rate(f, mean_per_period, sd_per_period):
    """Expected log-growth per period at fraction f: g ≈ f·μ − f²σ²/2.
    Note g(2f*) = 0 — twice Kelly grows at zero while carrying enormous variance."""
    return f * mean_per_period - 0.5 * (f * sd_per_period) ** 2

# ---------------------------------------------------------------- 39. drawdown-constrained


def prob_hit_drawdown(dd, f, mean_per_period, sd_per_period):
    """P(ever touching a relative drawdown `dd`) for a GBM-like process:
    ≈ exp(−2μ_f·D/σ_f²) with μ_f = f·μ, σ_f = f·σ. dd as a positive fraction (0.2 = 20%)."""
    mu, sg = f * mean_per_period, f * sd_per_period
    if mu <= 0 or sg == 0:
        return 1.0
    return min(1.0, math.exp(-2.0 * mu * dd / (sg ** 2)))


def dd_constrained_fraction(dd_tolerated, mean_per_period, sd_per_period, prob=0.10):
    """Largest f whose probability of ever hitting `dd_tolerated` stays ≤ prob.
    Inverting the barrier formula: f ≤ −2·μ·D / (σ²·ln prob).
    This is usually a much tighter (and more useful) constraint than Kelly."""
    if mean_per_period <= 0 or sd_per_period == 0 or not 0 < prob < 1:
        return 0.0
    f = (-2.0 * mean_per_period * dd_tolerated) / (sd_per_period ** 2 * math.log(prob))
    return max(0.0, f)


def sizing_summary(mean_per_period, sd_per_period, dd_tolerated=0.20, prob=0.10, fraction=0.25):
    if sd_per_period <= 0 or not math.isfinite(sd_per_period):
        return {"kelly_full": float('inf'), "kelly_fractional": float('inf'),
                "dd_constrained": 0.0, "binding_constraint": "undefined (zero variance)",
                "recommended_f": 0.0, "growth_at_recommended": 0.0}
    ks = kelly_continuous(mean_per_period, sd_per_period)
    dd = dd_constrained_fraction(dd_tolerated, mean_per_period, sd_per_period, prob)
    binding = "drawdown" if dd < fractional_kelly(ks, fraction) else "Kelly fraction"
    return {"kelly_full": ks, "kelly_fractional": fractional_kelly(ks, fraction),
            "dd_constrained": dd, "binding_constraint": binding,
            "recommended_f": min(fractional_kelly(ks, fraction), dd),
            "growth_at_recommended": growth_rate(min(fractional_kelly(ks, fraction), dd),
                                                 mean_per_period, sd_per_period)}

# ---------------------------------------------------------------- 40. vol targeting


def ewma_vol(returns, lam=0.94, seed_window=20):
    """EWMA volatility: σ²_t = λσ²_{t−1} + (1−λ)r²_t. λ=0.94 ≈ RiskMetrics daily."""
    if not returns:
        return float('nan')
    seed = returns[:min(seed_window, len(returns))]
    var = st.pvariance(seed) if len(seed) > 1 else (seed[0] ** 2 if seed else 0.0)
    for r in returns[len(seed):]:
        var = lam * var + (1 - lam) * r * r
    return math.sqrt(max(var, 0.0))


def vol_target_size(sigma_target, sigma_estimate, cap=3.0):
    """f_t = σ_target/σ̂_t, capped. The cap matters: a low realised vol (holiday, closed market,
    illiquidity) otherwise levers you up right before the gap."""
    if sigma_estimate <= 0:
        return 0.0
    return min(cap, sigma_target / sigma_estimate)

# ---------------------------------------------------------------- 41-42. impact & capacity


def impact_sqrt(size, daily_volume, sigma_daily, Y=0.8):
    """Square-root law: impact ≈ Y·σ·√(Q/V). Cost grows like √size, so doubling size does NOT
    double cost — but it does not stay flat either, which is what a linear cost model assumes."""
    if daily_volume <= 0:
        return float('inf')
    return Y * sigma_daily * math.sqrt(max(0.0, size) / daily_volume)


def capacity(alpha_per_unit, daily_volume, sigma_daily, Y=0.8):
    """Size at which expected alpha per unit is fully eaten by impact:
    solve alpha = Y·σ·√(Q/V) → Q* = V·(alpha/(Y·σ))².
    Caveat that kills naive estimates: V evaporates exactly when you need to trade."""
    if sigma_daily <= 0 or Y <= 0:
        return float('inf')
    return daily_volume * (alpha_per_unit / (Y * sigma_daily)) ** 2


def turnover_drag(gross_return, turnover, cost_per_turn):
    """net = gross − turnover × cost_per_turn, with the survivability ratio.
    An edge needs gross ≥ 2-3× total cost to survive the worst regimes."""
    cost = turnover * cost_per_turn
    net = gross_return - cost
    return {"gross": gross_return, "cost": cost, "net": net,
            "cost_share_pct": (100 * cost / gross_return) if gross_return else float('inf'),
            "survivability": (gross_return / cost) if cost else float('inf'),
            "verdict": ("fragile — cost is a large share of edge" if cost and gross_return / cost < 2
                        else "acceptable")}


def cost_ratio(cost_per_round_trip, target):
    """Cost as a fraction of the TARGET — the number that reveals whether a parameter change
    is a tweak or a different economic proposition."""
    if target == 0:
        return float('inf')
    r = cost_per_round_trip / abs(target)
    return {"ratio_pct": 100 * r,
            "verdict": ("fine (<10%)" if r < 0.10 else
                        "material (10-25%)" if r < 0.25 else "dominant (>25%) — rethink")}

# ---------------------------------------------------------------- 44. fills & latency


def effective_rule_distance(nominal_distance, poll_seconds, bar_range_per_second):
    """A path-dependent rule re-emitted every `poll_seconds` lags the true running extreme,
    so it behaves WIDER than its nominal distance. Returns the effective distance to test."""
    slack = bar_range_per_second * max(0.0, poll_seconds)
    return {"nominal": nominal_distance, "effective": nominal_distance + slack,
            "slack": slack,
            "note": "test the effective distance — it is what actually executes"}


def adverse_fill_adjustment(fill_prob_touch, adverse_share=0.5):
    """Limit orders fill preferentially when the price keeps going against you (adverse
    selection). A backtest assuming fill-on-touch overstates limit-order strategies.
    Returns the fill probability to assume for an unfavourable-continuation scenario."""
    return {"assumed_on_touch": fill_prob_touch,
            "adverse_adjusted": fill_prob_touch * (1 - adverse_share),
            "note": "if the edge disappears under this, it was a fill artefact"}


def report(mean_per_period, sd_per_period, dd_tolerated=0.20, gross=None,
           turnover=None, cost_per_turn=None, label="strategy"):
    s = sizing_summary(mean_per_period, sd_per_period, dd_tolerated)
    L = [f"=== sizing · {label} ===",
         f"  Kelly full {s['kelly_full']:.3f} · quarter-Kelly {s['kelly_fractional']:.3f}",
         f"  drawdown-constrained (≤{dd_tolerated:.0%} with 10% prob): {s['dd_constrained']:.3f}",
         f"  binding constraint: {s['binding_constraint']} → recommended f = {s['recommended_f']:.3f}",
         f"  expected log-growth at that f: {s['growth_at_recommended']:+.5f}/period"]
    if gross is not None and turnover is not None and cost_per_turn is not None:
        t = turnover_drag(gross, turnover, cost_per_turn)
        L.append(f"  costs: gross {t['gross']:+.4f} − cost {t['cost']:.4f} = net {t['net']:+.4f} "
                 f"· survivability {t['survivability']:.1f}× → {t['verdict']}")
    return "\n".join(L)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mean", type=float, required=True)
    ap.add_argument("--sd", type=float, required=True)
    ap.add_argument("--dd", type=float, default=0.20)
    a = ap.parse_args()
    print(report(a.mean, a.sd, a.dd))
