"""PyMatching baseline and herald-adaptive decoder."""

from __future__ import annotations

import math

import numpy as np
import pytest

from leakdecode.circuits import build_memory
from leakdecode.decoders import AdaptiveMatcher, BaselineMatcher
from leakdecode.sim import sample


def edge_snapshot(matching) -> tuple:
    """Canonical list of ``(u, v, weight, probability, fault ids)``."""
    rows = []
    for u, v, data in matching.edges():
        rows.append(
            (
                int(u),
                None if v is None else int(v),
                float(data["weight"]),
                float(data["error_probability"]),
                frozenset(int(fault) for fault in data["fault_ids"]),
            )
        )
    return tuple(sorted(rows, key=lambda row: (row[0], -1 if row[1] is None else row[1], row[2], row[3])))


@pytest.mark.parametrize("distance", [3, 5])
def test_zero_noise_has_no_detections_and_no_logical_errors(distance: int):
    shots = 32
    mem = build_memory(distance, distance, 0.0, 0.0)
    dets, obs, heralds = sample(mem, shots, seed=0)
    assert int(dets.sum()) == 0
    assert int(obs.sum()) == 0
    assert int(heralds.sum()) == 0

    baseline = BaselineMatcher(mem)
    adaptive = AdaptiveMatcher(mem)
    pred_base = baseline.decode_batch(dets)
    pred_adapt = adaptive.decode_batch(dets, heralds)
    assert pred_base.shape == (shots,)
    assert pred_adapt.shape == (shots,)
    assert np.array_equal(pred_base, obs)
    assert np.array_equal(pred_adapt, obs)
    assert int(pred_base.sum()) == 0
    assert int(pred_adapt.sum()) == 0


def test_zero_leak_predictions_agree():
    shots = 1500
    mem = build_memory(3, 3, 2e-3, 0.0)
    dets, _obs, heralds = sample(mem, shots, seed=42)
    assert int(heralds.sum()) == 0
    assert int(dets.sum()) > 0

    pred_base = BaselineMatcher(mem).decode_batch(dets)
    pred_adapt = AdaptiveMatcher(mem).decode_batch(dets, heralds)
    assert np.array_equal(pred_base, pred_adapt)


def test_pure_erasure_adaptive_is_correct_on_light_patterns():
    """Pure heralded erasure, no other circuit noise.

    A logical operator uses at least ``d`` decoding-graph edges, and one herald
    contributes at most two edges. Fewer than ``d - 1`` fired heralds is the
    safe count asked for here: with ``d = 3`` that is at most one herald, so
    at most two edges. Fewer than ``d`` distinct erased edges is the same
    guarantee stated directly on the graph. The adaptive decoder has to get
    every such shot right.

    Across the whole fixed sample the adaptive decoder should also make
    strictly fewer logical errors than the static decoder, and the adaptive
    error fraction should stay small. Shots with many heralds can be
    uncorrectable, so that comparison is statistical.
    """
    distance = 3
    shots = 4000
    p_leak = 0.02
    mem = build_memory(distance, distance, 0.0, p_leak)
    dets, obs, heralds = sample(mem, shots, seed=123)
    baseline = BaselineMatcher(mem)
    adaptive = AdaptiveMatcher(mem)
    pred_base = baseline.decode_batch(dets)
    pred_adapt = adaptive.decode_batch(dets, heralds)

    static_errors = int((pred_base != obs).sum())
    adapt_errors = int((pred_adapt != obs).sum())
    assert static_errors > 0
    assert adapt_errors < static_errors
    assert adapt_errors / shots <= 0.02

    light_herald_shots = 0
    light_herald_with_fire = 0
    light_edge_shots = 0
    for shot in range(shots):
        fired = np.flatnonzero(heralds[shot])
        erased = {
            edge
            for herald_index in fired
            for edge in adaptive.herald_map.edges[int(herald_index)]
        }
        if len(fired) < distance - 1:
            light_herald_shots += 1
            if len(fired):
                light_herald_with_fire += 1
            assert int(pred_adapt[shot]) == int(obs[shot])
        if len(erased) < distance:
            light_edge_shots += 1
            assert int(pred_adapt[shot]) == int(obs[shot])

    assert light_herald_shots > 0
    assert light_herald_with_fire > 0
    assert light_edge_shots > 0


def test_adaptive_within_two_percent_of_static_under_circuit_noise():
    """Moderate circuit noise and leakage. Adaptive may not exceed static by 2%."""
    shots = 3000
    mem = build_memory(3, 3, 2e-3, 1e-2)
    dets, obs, heralds = sample(mem, shots, seed=12345)
    assert int(dets.sum()) > 0

    pred_base = BaselineMatcher(mem).decode_batch(dets, heralds)
    pred_adapt = AdaptiveMatcher(mem).decode_batch(dets, heralds)
    static_errors = int((pred_base != obs).sum())
    adapt_errors = int((pred_adapt != obs).sum())
    assert adapt_errors / shots <= static_errors / shots + 0.02
    assert adapt_errors <= static_errors + math.ceil(0.02 * shots)


def test_decode_one_restores_static_weights():
    shots = 50
    mem = build_memory(3, 3, 2e-3, 0.05)
    dets, _obs, heralds = sample(mem, shots, seed=11)
    fired = heralds.any(axis=1)
    assert int(fired.sum()) > 0
    assert int((~fired).sum()) > 0

    adaptive = AdaptiveMatcher(mem)
    before = edge_snapshot(adaptive.matching)
    pred_batch = adaptive.decode_batch(dets, heralds)
    assert edge_snapshot(adaptive.matching) == before

    for shot in range(shots):
        first = adaptive.decode_one(dets[shot], heralds[shot])
        assert edge_snapshot(adaptive.matching) == before
        second = adaptive.decode_one(dets[shot], heralds[shot])
        assert first == second
        assert first == int(pred_batch[shot])
        assert edge_snapshot(adaptive.matching) == before

    # A heralded shot must actually take the reweight path.
    heralded = int(np.flatnonzero(fired)[0])
    assert adaptive.erased_edges(heralds[heralded]).size > 0

    fresh = AdaptiveMatcher(mem)
    assert edge_snapshot(fresh.matching) == before

    baseline = BaselineMatcher(mem)
    base_before = edge_snapshot(baseline.matching)
    baseline.decode_batch(dets)
    assert edge_snapshot(baseline.matching) == base_before
