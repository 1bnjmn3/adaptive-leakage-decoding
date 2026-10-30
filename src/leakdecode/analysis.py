"""Turn the sweep into fits and an (extrapolated) qubit-saving estimate."""

from __future__ import annotations

import math
from collections import defaultdict

from .circuits import num_qubits_rotated
from .stats import distance_for_target, fit_exponential_suppression

# (kind, value). "per_round": logical error per round = value.
# "window": the LCD paper's target - failure probability value for one
# d x d x d block (d rounds), approximated as d * eps_round(d) = value.
TARGETS = (("per_round", 1e-6), ("per_round", 1e-9), ("per_round", 1e-12), ("window", 1e-6))


def _distance(a, b, kind, t):
    if kind == "per_round":
        return distance_for_target(a, b, t)
    # Solve log10(d) + a + b d = log10(t). The left side rises until its peak at
    # d* = -1/(b ln 10) and falls after it, so for shallow slopes it can cross the
    # target twice. The distance we want is the LARGER root (beyond it the window
    # stays under target), found by bisection on the decreasing branch [d*, 1000].
    if b >= 0:
        return math.inf
    f = lambda d: math.log10(d) + a + b * d - math.log10(t)
    lo, hi = max(1.0, -1.0 / (b * math.log(10))), 1000.0
    if f(lo) <= 0:
        return 1.0  # target met at every d >= 1
    if f(hi) > 0:
        return math.inf
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fits_from_sweep(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(float(r["p_leak"]), r["decoder"])].append(r)
    out = []
    for (pl, dec), rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: int(r["distance"]))
        f = fit_exponential_suppression(
            [int(r["distance"]) for r in rs], [float(r["ler_per_round"]) for r in rs],
            [float(r["ler_per_round_lo"]) for r in rs], [float(r["ler_per_round_hi"]) for r in rs])
        a, b, lam = f if f else (math.nan, math.nan, math.nan)
        out.append(dict(p_leak=pl, decoder=dec, a=a, b=b, Lambda=lam,
                        n_points_fitted=sum(1 for r in rs if int(r["errors"]) > 0),
                        distances=" ".join(r["distance"] for r in rs)))
    return out


def savings_from_fits(fits, targets=TARGETS):
    by = {(f["p_leak"], f["decoder"]): f for f in fits}
    out = []
    for pl in sorted({f["p_leak"] for f in fits}):
        s, a = by.get((pl, "pymatching_static")), by.get((pl, "pymatching_adaptive"))
        if not s or not a:
            continue
        for kind, t in targets:
            row = dict(p_leak=pl, target_kind=kind, target=t)
            for name, f in (("static", s), ("adaptive", a)):
                d = _distance(f["a"], f["b"], kind, t) if not math.isnan(f["a"]) else math.nan
                row[f"d_{name}"] = d
                # round up to the next odd distance a real device could use
                d_odd = (math.ceil((d - 1) / 2) * 2 + 1) if math.isfinite(d) else math.nan
                row[f"d_{name}_odd"] = max(3, d_odd) if math.isfinite(d_odd) else math.nan
                row[f"qubits_{name}"] = 2 * d * d - 1 if math.isfinite(d) else math.nan
                row[f"qubits_{name}_odd"] = (num_qubits_rotated(int(row[f"d_{name}_odd"]))
                                             if math.isfinite(d_odd) else math.nan)
            row["qubit_ratio"] = (row["qubits_static"] / row["qubits_adaptive"]
                                  if math.isfinite(row["qubits_static"]) and math.isfinite(row["qubits_adaptive"]) else math.nan)
            row["qubit_ratio_odd"] = (row["qubits_static_odd"] / row["qubits_adaptive_odd"]
                                      if math.isfinite(row["qubits_static_odd"]) and math.isfinite(row["qubits_adaptive_odd"]) else math.nan)
            out.append(row)
    return out
