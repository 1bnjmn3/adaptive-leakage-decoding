"""Circuit layout: herald counts, measurement record, detectors, and the DEM."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pytest
import stim

from leakdecode.circuits import build_memory, num_qubits_rotated
from leakdecode.sim import sample


def plain_memory(distance: int, rounds: int, p: float) -> stim.Circuit:
    """Stim rotated Z-memory with the study's circuit-level noise, and no heralds."""
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


def herald_record_indices(circuit: stim.Circuit) -> list[int]:
    """Absolute measurement-record indices written by HERALDED_ERASE."""
    indices: list[int] = []
    cursor = 0
    for inst in circuit.flattened():
        if inst.name == "HERALDED_ERASE":
            for _target in inst.targets_copy():
                indices.append(cursor)
                cursor += 1
        else:
            cursor += inst.num_measurements
    assert cursor == circuit.num_measurements
    return indices


def error_mechanisms(circuit: stim.Circuit) -> Counter:
    """Multiset of DEM error mechanisms, ignoring instruction order.

    Each mechanism is ``(probability, frozenset of components)``. A component
    is the frozenset of detector and observable targets between ``^`` separators,
    so the order of targets and of decomposition components does not matter.
    """
    dem = circuit.detector_error_model(
        decompose_errors=True,
        approximate_disjoint_errors=True,
    )
    counts: Counter = Counter()
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        probability = inst.args_copy()[0]
        components: list[frozenset] = []
        current: list[tuple[str, int]] = []
        for target in inst.targets_copy():
            if target.is_separator():
                components.append(frozenset(current))
                current = []
            elif target.is_relative_detector_id():
                current.append(("D", int(target.val)))
            elif target.is_logical_observable_id():
                current.append(("L", int(target.val)))
            else:
                raise AssertionError(f"unexpected DEM target {target}")
        components.append(frozenset(current))
        counts[(probability, frozenset(components))] += 1
    return counts


@pytest.mark.parametrize("distance", [0, 1, 2, 4, 6, -3])
def test_invalid_distance_raises(distance: int):
    with pytest.raises(ValueError):
        build_memory(distance, 3, 0.0, 0.0)


@pytest.mark.parametrize("p,p_leak", [(-0.01, 0.0), (1.01, 0.0), (0.0, -0.01), (0.0, 1.01)])
def test_probability_outside_unit_interval_raises(p: float, p_leak: float):
    with pytest.raises(ValueError):
        build_memory(3, 3, p, p_leak)


def test_unknown_leak_location_raises():
    with pytest.raises(ValueError):
        build_memory(3, 3, 0.0, 0.0, leak_on=("data", "bulk"))


@pytest.mark.parametrize("distance", [3, 5])
def test_herald_counts_match_the_patch(distance: int):
    rounds = distance
    mem = build_memory(distance, rounds, 2e-3, 1e-3)
    data_per_round = distance * distance
    ancilla_per_round = distance * distance - 1

    assert mem.data_qubits == sorted(mem.data_qubits)
    assert len(mem.data_qubits) == data_per_round
    assert len(mem.ancilla_qubits) == ancilla_per_round
    assert mem.num_physical_qubits == num_qubits_rotated(distance)
    assert mem.num_heralds == (data_per_round + ancilla_per_round) * rounds

    for round_index in range(rounds):
        in_round = mem.herald_round == round_index
        data = in_round & (mem.herald_kind == 0)
        ancilla = in_round & (mem.herald_kind == 1)
        assert int(data.sum()) == data_per_round
        assert int(ancilla.sum()) == ancilla_per_round
        assert sorted(mem.herald_qubit[data].tolist()) == mem.data_qubits
        assert sorted(mem.herald_qubit[ancilla].tolist()) == mem.ancilla_qubits

    assert set(np.unique(mem.herald_kind)).issubset({0, 1})


@pytest.mark.parametrize("distance", [3, 5])
def test_herald_meas_index_points_at_heralded_erase(distance: int):
    mem = build_memory(distance, distance, 2e-3, 0.25)
    assert herald_record_indices(mem.circuit) == mem.herald_meas_index.tolist()
    assert mem.herald_meas_index.min() >= 0
    assert mem.herald_meas_index.max() < mem.circuit.num_measurements
    # One new record per herald, inserted in circuit order.
    assert np.all(np.diff(mem.herald_meas_index) > 0)


@pytest.mark.parametrize("distance", [3, 5])
@pytest.mark.parametrize("p", [0.0, 2e-3])
@pytest.mark.parametrize("p_leak,all_fired", [(0.0, False), (1.0, True)])
def test_herald_bits_are_deterministic_at_the_extremes(distance, p, p_leak, all_fired):
    mem = build_memory(distance, distance, p, p_leak)
    _dets, _obs, heralds = sample(mem, 16, seed=1)
    if all_fired:
        assert bool(heralds.all())
    else:
        assert not bool(heralds.any())


@pytest.mark.parametrize("distance", [3, 5])
@pytest.mark.parametrize("p_leak", [0.0, 1e-3, 0.2])
def test_detector_count_matches_the_plain_circuit(distance: int, p_leak: float):
    p = 2e-3
    rounds = distance
    mem = build_memory(distance, rounds, p, p_leak)
    plain = plain_memory(distance, rounds, p)
    assert mem.circuit.num_detectors == plain.num_detectors
    assert mem.circuit.num_observables == plain.num_observables
    assert mem.circuit.get_detector_coordinates() == plain.get_detector_coordinates()


@pytest.mark.parametrize("distance", [3, 5])
@pytest.mark.parametrize("p", [0.0, 2e-3])
def test_dem_at_zero_leak_matches_the_plain_circuit(distance: int, p: float):
    mem = build_memory(distance, distance, p, 0.0)
    plain = plain_memory(distance, distance, p)
    assert error_mechanisms(mem.circuit) == error_mechanisms(plain)


def test_with_leak_rate_keeps_the_herald_layout():
    mem = build_memory(3, 3, 2e-3, 0.0)
    shifted = mem.with_leak_rate(0.01)
    assert shifted.p_leak == 0.01
    assert shifted.p == mem.p
    assert shifted.distance == mem.distance
    assert np.array_equal(shifted.herald_meas_index, mem.herald_meas_index)
    assert shifted.circuit.num_detectors == mem.circuit.num_detectors
    assert shifted.num_heralds == mem.num_heralds
