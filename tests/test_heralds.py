"""Herald signatures: graph-like edges, final-round X edges, and erasure span."""

from __future__ import annotations

import numpy as np
import pytest

from leakdecode.circuits import build_memory
from leakdecode.heralds import BOUNDARY, Edge, herald_signatures
from leakdecode.sim import sample


def _detectors(edge: Edge) -> list[int]:
    found = [edge.u]
    if edge.v != BOUNDARY:
        found.append(edge.v)
    return found


def syndrome_in_span(syndrome: np.ndarray, edges: list[Edge], num_detectors: int) -> bool:
    """True when ``syndrome`` is a GF(2) sum of the given graph edges.

    An edge ``(u, v)`` is the vector with those two detectors set. A boundary
    edge sets one detector. Gaussian elimination decides membership.
    """
    rows: list[int] = []
    for edge in edges:
        mask = 1 << int(edge.u)
        if edge.v != BOUNDARY:
            mask |= 1 << int(edge.v)
        rows.append(mask)
    target = 0
    for index, bit in enumerate(syndrome):
        if bit:
            target |= 1 << int(index)
    used = [False] * len(rows)
    for col in range(num_detectors):
        pivot = None
        for i, row in enumerate(rows):
            if not used[i] and (row >> col) & 1:
                pivot = i
                break
        if pivot is None:
            if (target >> col) & 1:
                return False
            continue
        used[pivot] = True
        pivot_row = rows[pivot]
        for j in range(len(rows)):
            if j != pivot and (rows[j] >> col) & 1:
                rows[j] ^= pivot_row
        if (target >> col) & 1:
            target ^= pivot_row
    return target == 0


def observable_in_span(syndrome: np.ndarray, observable: int, edges: list[Edge], num_detectors: int) -> bool:
    """Syndrome plus the logical bit, with each edge's ``obs`` flag as an extra column."""
    obs_col = num_detectors
    rows: list[int] = []
    for edge in edges:
        mask = 1 << int(edge.u)
        if edge.v != BOUNDARY:
            mask |= 1 << int(edge.v)
        if edge.obs:
            mask |= 1 << obs_col
        rows.append(mask)
    target = (1 << obs_col) if observable else 0
    for index, bit in enumerate(syndrome):
        if bit:
            target |= 1 << int(index)
    used = [False] * len(rows)
    for col in range(num_detectors + 1):
        pivot = None
        for i, row in enumerate(rows):
            if not used[i] and (row >> col) & 1:
                pivot = i
                break
        if pivot is None:
            if (target >> col) & 1:
                return False
            continue
        used[pivot] = True
        pivot_row = rows[pivot]
        for j in range(len(rows)):
            if j != pivot and (rows[j] >> col) & 1:
                rows[j] ^= pivot_row
        if (target >> col) & 1:
            target ^= pivot_row
    return target == 0


def test_span_checker_accepts_known_vectors():
    edge_bulk = Edge(0, 1, 0)
    edge_boundary = Edge(0, BOUNDARY, 1)
    assert syndrome_in_span(np.array([0, 0, 0]), [], 3)
    assert syndrome_in_span(np.array([1, 1, 0]), [edge_bulk], 3)
    assert not syndrome_in_span(np.array([1, 0, 0]), [edge_bulk], 3)
    assert syndrome_in_span(np.array([1, 0]), [edge_boundary], 2)
    assert syndrome_in_span(np.array([1, 0, 0]), [edge_bulk, Edge(1, BOUNDARY, 0)], 3)
    assert not syndrome_in_span(np.array([0, 0, 1]), [edge_bulk, edge_boundary], 3)
    assert observable_in_span(np.array([1, 0]), 1, [edge_boundary], 2)
    assert not observable_in_span(np.array([1, 0]), 0, [edge_boundary], 2)


@pytest.mark.parametrize("distance", [3, 5])
def test_every_herald_edge_is_graph_like(distance: int):
    rounds = distance
    mem = build_memory(distance, rounds, 0.0, 0.0)
    herald_map = herald_signatures(mem)
    num_detectors = mem.circuit.num_detectors
    coords = mem.circuit.get_detector_coordinates()

    assert len(herald_map.edges) == mem.num_heralds

    logical_qubits_by_round: list[set[int]] = []
    for round_index in range(rounds):
        data_index = np.flatnonzero((mem.herald_kind == 0) & (mem.herald_round == round_index))
        logical: set[int] = set()
        for herald_index in data_index:
            edges = herald_map.edges[int(herald_index)]
            # Round 0 still has data in |0>, so a Z at the herald is trivial and
            # only the X edge remains. Later rounds keep both sectors.
            if round_index == 0:
                assert len(edges) == 1
            else:
                assert len(edges) == 2
                x_edge, z_edge = edges
                assert z_edge.obs == 0
                assert set(_detectors(x_edge)).isdisjoint(_detectors(z_edge))
            if any(edge.obs for edge in edges):
                assert edges[0].obs == 1
                logical.add(int(mem.herald_qubit[int(herald_index)]))
        logical_qubits_by_round.append(logical)

        ancilla_index = np.flatnonzero((mem.herald_kind == 1) & (mem.herald_round == round_index))
        assert len(ancilla_index) == distance * distance - 1
        for herald_index in ancilla_index:
            edges = herald_map.edges[int(herald_index)]
            # The ancilla is measured immediately, so one Pauli sector flips that
            # record and the other commutes with the measurement.
            assert len(edges) == 1
            assert edges[0].obs == 0

    assert logical_qubits_by_round
    assert all(qubits == logical_qubits_by_round[0] for qubits in logical_qubits_by_round)
    assert len(logical_qubits_by_round[0]) == distance

    # Final-round data X edges: the first of the two edges. They sit on the last
    # ancilla round (the partner that would compare to the final data measurement
    # cancels), and the d qubits of the logical Z string give a boundary edge.
    final = np.flatnonzero((mem.herald_kind == 0) & (mem.herald_round == rounds - 1))
    assert len(final) == distance * distance
    final_logical = 0
    for herald_index in final:
        x_edge, z_edge = herald_map.edges[int(herald_index)]
        for edge in (x_edge, z_edge):
            times = {coords[det][2] for det in _detectors(edge)}
            assert times == {float(rounds - 1)}
        assert z_edge.obs == 0
        if x_edge.obs:
            final_logical += 1
            assert x_edge.v == BOUNDARY
            assert int(mem.herald_qubit[int(herald_index)]) in logical_qubits_by_round[0]
    assert final_logical == distance

    for edges in herald_map.edges:
        assert 1 <= len(edges) <= 2
        for edge in edges:
            detectors = _detectors(edge)
            assert 1 <= len(detectors) <= 2
            assert len(set(detectors)) == len(detectors)
            assert all(0 <= det < num_detectors for det in detectors)
            assert edge.obs in (0, 1)
            if edge.v != BOUNDARY:
                assert edge.u < edge.v


def test_pure_erasure_syndrome_lies_in_the_herald_span():
    """With circuit noise off, every detection event comes from a fired herald."""
    distance = 3
    mem = build_memory(distance, distance, 0.0, 0.05)
    herald_map = herald_signatures(mem)
    shots = 300
    dets, obs, heralds = sample(mem, shots, seed=5)
    num_detectors = mem.circuit.num_detectors
    assert int(dets.sum()) > 0

    for shot in range(shots):
        fired = np.flatnonzero(heralds[shot])
        edges: list[Edge] = []
        for herald_index in fired:
            edges.extend(herald_map.edges[int(herald_index)])
        assert syndrome_in_span(dets[shot], edges, num_detectors)
        assert observable_in_span(dets[shot], int(obs[shot]), edges, num_detectors)
