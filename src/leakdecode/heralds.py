"""Map every herald bit to the decoding-graph edges it can affect.

An erasure at location k applies I, X, Y or Z with probability 1/4 each.
For a Z-basis memory the decoding graph splits into two independent pieces:

* the "Z-detector" graph, which sees X errors (and the X half of Y), and
* the "X-detector" graph, which sees Z errors (and the Z half of Y).

So erasure k can only flip (at most) one edge in each piece:
``edge_X(k)`` = the detectors flipped by an X at k, and
``edge_Z(k)`` = the detectors flipped by a Z at k.
Given that the herald fired, each of those edges flips with probability exactly
1/2 (X or Y flips ``edge_X``; Z or Y flips ``edge_Z``).

We find the signatures exactly, with no guessing about geometry, by running a
noiseless Stim ``FlipSimulator`` with one batch instance per (herald, Pauli),
injecting that Pauli right after the HERALDED_ERASE instruction, and reading
back which detectors and which logical observable flipped.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import stim

from .circuits import LeakyMemory, build_memory

BOUNDARY = -1


@dataclass(frozen=True)
class Edge:
    u: int
    v: int  # BOUNDARY (-1) for a boundary edge
    obs: int  # 1 if this edge flips logical observable 0

    @staticmethod
    def from_dets(dets: np.ndarray, obs: int) -> "Edge | None":
        dets = sorted(int(x) for x in dets)
        if len(dets) == 0:
            if obs:
                raise ValueError("undetectable logical error: code distance too small?")
            return None
        if len(dets) == 1:
            return Edge(dets[0], BOUNDARY, int(obs))
        if len(dets) == 2:
            return Edge(dets[0], dets[1], int(obs))
        raise ValueError(f"erasure signature is not graph-like: {dets}")


@dataclass
class HeraldMap:
    """For herald k: ``edges[k]`` is a tuple of 0, 1 or 2 graph edges."""

    edges: list[tuple[Edge, ...]]

    def unique_edges(self) -> list[Edge]:
        seen: dict[Edge, None] = {}
        for es in self.edges:
            for e in es:
                seen.setdefault(e, None)
        return list(seen)


def herald_signatures(mem: LeakyMemory) -> HeraldMap:
    """Exact detector/observable signature of X and Z at every herald location."""
    # Same layout, but with every noise channel switched off.
    clean = build_memory(mem.distance, mem.rounds, 0.0, 0.0, mem.leak_on).circuit
    k_total = mem.num_heralds
    batch = 2 * k_total  # instance 2k: X at herald k; instance 2k+1: Z at herald k
    nq = clean.num_qubits
    sim = stim.FlipSimulator(
        batch_size=batch, disable_stabilizer_randomization=True, num_qubits=nq
    )
    k = 0
    for inst in clean:  # circuit is already flat
        sim.do(inst)
        if inst.name == "HERALDED_ERASE":
            mask_x = np.zeros((nq, batch), dtype=np.bool_)
            mask_z = np.zeros((nq, batch), dtype=np.bool_)
            for t in inst.targets_copy():
                mask_x[t.value, 2 * k] = True
                mask_z[t.value, 2 * k + 1] = True
                k += 1
            sim.broadcast_pauli_errors(pauli="X", mask=mask_x)
            sim.broadcast_pauli_errors(pauli="Z", mask=mask_z)
    assert k == k_total
    det = sim.get_detector_flips()  # (num_detectors, batch)
    obs = sim.get_observable_flips()  # (num_observables, batch)
    out: list[tuple[Edge, ...]] = []
    for k in range(k_total):
        es = []
        for col in (2 * k, 2 * k + 1):
            e = Edge.from_dets(np.flatnonzero(det[:, col]), int(obs[0, col]))
            if e is not None:
                es.append(e)
        out.append(tuple(es))
    return HeraldMap(out)
