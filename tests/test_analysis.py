"""Qubit-saving ratios from a pair of exponential fits."""

from __future__ import annotations

import math

from leakdecode.analysis import TARGETS, savings_from_fits
from leakdecode.stats import distance_for_target


def _fit(decoder: str, a: float = -2.0, b: float = -0.5, p_leak: float = 0.01) -> dict:
    return dict(p_leak=p_leak, decoder=decoder, a=a, b=b, Lambda=10 ** (-2 * b))


def test_savings_from_identical_fits_has_unit_ratio():
    fits = [_fit("pymatching_static"), _fit("pymatching_adaptive")]
    rows = savings_from_fits(fits)
    assert len(rows) == len(TARGETS)
    for row in rows:
        assert row["qubit_ratio"] == 1.0
        assert row["qubit_ratio_odd"] == 1.0
        assert row["d_static"] == row["d_adaptive"]
        assert row["d_static_odd"] == row["d_adaptive_odd"]
        assert row["qubits_static"] == row["qubits_adaptive"]
        assert math.isfinite(row["d_static"])

    per_round_only = savings_from_fits(fits, targets=(("per_round", 1e-6),))
    assert len(per_round_only) == 1
    expected_d = distance_for_target(-2.0, -0.5, 1e-6)
    assert per_round_only[0]["d_static"] == expected_d
    assert per_round_only[0]["qubit_ratio"] == 1.0
    assert per_round_only[0]["target_kind"] == "per_round"
    assert per_round_only[0]["target"] == 1e-6


def test_steeper_adaptive_fit_uses_fewer_qubits():
    fits = [
        _fit("pymatching_static", b=-0.2),
        _fit("pymatching_adaptive", b=-0.5),
    ]
    rows = savings_from_fits(fits, targets=(("per_round", 1e-9),))
    assert len(rows) == 1
    assert rows[0]["d_adaptive"] < rows[0]["d_static"]
    assert rows[0]["qubit_ratio"] > 1.0


def test_savings_skips_a_leak_rate_missing_a_decoder():
    assert savings_from_fits([_fit("pymatching_static")]) == []


def test_window_distance_takes_the_larger_root():
    """Shallow slope: d*eps(d) crosses the target twice; we need the larger root."""
    import math
    from leakdecode.analysis import _distance
    a, b, t = -5.9, -0.15, 1e-6
    d = _distance(a, b, "window", t)
    assert 5.5 < d < 6.0
    assert abs(math.log10(d) + a + b * d - math.log10(t)) < 1e-9
    # beyond the returned root the window stays under target
    assert all(math.log10(x) + a + b * x < math.log10(t) for x in (d + 0.1, d + 5, d + 50))
