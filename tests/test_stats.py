"""Wilson intervals, per-round rates, and the exponential fit."""

from __future__ import annotations

import math

import pytest

from leakdecode.stats import Z_95, fit_exponential_suppression, per_round, wilson_interval


@pytest.mark.parametrize("k,n", [(0, 1), (0, 10), (1, 2), (1, 10), (5, 20), (9, 10), (10, 10), (3, 3), (2, 100)])
def test_wilson_interval_contains_the_point_estimate(k: int, n: int):
    lo, hi = wilson_interval(k, n)
    phat = k / n
    assert 0.0 <= lo <= hi <= 1.0
    # One ulp of slack: when k == n the closed form is 1, and the float
    # evaluation can land just under 1 before the clip.
    assert lo <= phat + 1e-9
    assert phat <= hi + 1e-9


def test_wilson_interval_known_values():
    assert wilson_interval(0, 0) == (0.0, 1.0)

    # Zero successes: lower endpoint is 0 and upper = z^2 / (n + z^2).
    # The subtraction centre - half is 0 in real arithmetic; allow one ulp.
    n = 25
    expected_hi = Z_95 * Z_95 / (n + Z_95 * Z_95)
    lo, hi = wilson_interval(0, n)
    assert math.isclose(lo, 0.0, abs_tol=1e-12)
    assert math.isclose(hi, expected_hi, rel_tol=1e-12)

    lo_full, hi_full = wilson_interval(n, n)
    assert math.isclose(hi_full, 1.0, rel_tol=1e-12, abs_tol=1e-12)
    assert lo_full > 0.0


@pytest.mark.parametrize("p_shot", [0.0, 0.25, 0.5])
def test_per_round_of_one_round_is_the_shot_rate(p_shot: float):
    assert per_round(p_shot, 1) == p_shot


def test_per_round_inverts_the_shot_formula():
    def shot_rate(error_per_round: float, rounds: int) -> float:
        return (1.0 - (1.0 - 2.0 * error_per_round) ** rounds) / 2.0

    for error_per_round, rounds in [(0.0, 4), (1e-4, 3), (0.02, 5), (0.1, 7), (0.25, 1), (0.5, 6)]:
        recovered = per_round(shot_rate(error_per_round, rounds), rounds)
        assert math.isclose(recovered, error_per_round, rel_tol=1e-12, abs_tol=1e-15)

    # Values outside [0, 1/2] are clipped before the inverse is applied.
    assert per_round(1.0, 4) == per_round(0.5, 4)
    assert per_round(-0.2, 4) == 0.0


def test_fit_exponential_suppression_recovers_exact_line():
    a_true = -1.7
    b_true = -0.35
    distances = [3, 5, 7, 9]
    eps = [10 ** (a_true + b_true * d) for d in distances]
    lo = [value / 3 for value in eps]
    hi = [value * 3 for value in eps]
    fitted = fit_exponential_suppression(distances, eps, lo, hi)
    assert fitted is not None
    a, b, lam = fitted
    assert a == pytest.approx(a_true, abs=1e-8)
    assert b == pytest.approx(b_true, abs=1e-8)
    assert lam == pytest.approx(10 ** (-2 * b_true), rel=1e-8)

    # A point with no observed errors is left out; the other points still fix the line.
    eps_with_zero = list(eps)
    lo_with_zero = list(lo)
    hi_with_zero = list(hi)
    eps_with_zero[1] = 0.0
    lo_with_zero[1] = 0.0
    hi_with_zero[1] = 0.0
    fitted_skip = fit_exponential_suppression(distances, eps_with_zero, lo_with_zero, hi_with_zero)
    assert fitted_skip is not None
    a_skip, b_skip, _lam_skip = fitted_skip
    assert a_skip == pytest.approx(a_true, abs=1e-8)
    assert b_skip == pytest.approx(b_true, abs=1e-8)


def test_fit_exponential_suppression_needs_two_positive_points():
    assert fit_exponential_suppression([3, 5], [0.0, 0.01], [0.0, 0.001], [0.0, 0.05]) is None
    assert fit_exponential_suppression([3], [0.01], [0.001], [0.05]) is None
