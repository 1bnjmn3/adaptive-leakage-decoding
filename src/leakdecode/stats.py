"""Small statistics helpers (no scipy needed)."""

from __future__ import annotations

import math

import numpy as np

Z_95 = 1.959963984540054


def wilson_interval(k: int, n: int, z: float = Z_95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n (95% by default)."""
    if n == 0:
        return (0.0, 1.0)
    phat = k / n
    denom = 1 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def per_round(p_shot: float, rounds: int) -> float:
    """Convert a whole-experiment logical error probability to a per-round rate.

    If each round independently flips the logical with probability e, after r
    rounds the flip probability is P = (1 - (1 - 2e)^r) / 2. Inverting:
    e = (1 - (1 - 2P)^(1/r)) / 2. (Valid for P < 1/2.)
    """
    p_shot = min(max(p_shot, 0.0), 0.5)
    return 0.5 * (1 - (1 - 2 * p_shot) ** (1.0 / rounds))


def fit_exponential_suppression(distances, eps, eps_lo, eps_hi):
    """Fit log10(eps) = a + b*d by weighted least squares.

    Weights come from the width of the Wilson interval in log space. Returns
    (a, b, Lambda) where Lambda = 10^(-2b) is the suppression factor per
    distance step of 2. Points with zero observed errors are skipped.
    """
    d = np.asarray(distances, float)
    e = np.asarray(eps, float)
    lo = np.asarray(eps_lo, float)
    hi = np.asarray(eps_hi, float)
    ok = (e > 0) & (lo > 0)
    if ok.sum() < 2:
        return None
    y = np.log10(e[ok])
    sigma = (np.log10(hi[ok]) - np.log10(lo[ok])) / (2 * Z_95)
    w = 1 / np.maximum(sigma, 1e-6) ** 2
    A = np.vstack([np.ones(ok.sum()), d[ok]]).T
    W = np.diag(w)
    a, b = np.linalg.solve(A.T @ W @ A, A.T @ W @ y)
    return float(a), float(b), float(10 ** (-2 * b))


def distance_for_target(a: float, b: float, target: float) -> float:
    """Continuous code distance d solving a + b*d = log10(target)."""
    if b >= 0:
        return math.inf
    return (math.log10(target) - a) / b
