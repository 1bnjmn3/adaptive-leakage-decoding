"""Decoders compared in this project.

* :class:`BaselineMatcher` - PyMatching on the detector error model of the
  *leaky* circuit, ignoring heralds. Stim can only build this DEM with
  ``approximate_disjoint_errors=True``; it then folds each erasure into
  graph edges with the erasure's marginal flip probability (p_leak/2 per edge
  it can flip), dropping the X/Z-edge correlation. So
  this decoder knows leakage happens "on average" but not *where* in a shot.

* :class:`AdaptiveMatcher` - PyMatching that reads the herald bits of each
  shot. Its static graph is the DEM *without* erasures (if a herald did not
  fire, no erasure happened there). For every herald that fired, the (at most
  two) edges that erasure can flip are set to error probability 1/2, i.e.
  weight log((1-0.5)/0.5) = 0: the matching may use them for free. This is the
  standard "erasure-aware MWPM" step, the simplest form of the adaptivity the
  Local Clustering Decoder paper does in hardware.

Why p = 1/2 exactly: given a herald, the qubit got I, X, Y or Z with 1/4 each.
The Z-detector-graph edge flips for X or Y (prob 1/2); the X-detector-graph edge
flips for Z or Y (prob 1/2). Any other mechanism on the same edge, with
probability q, combines as q(1-1/2) + (1-q)/2 = 1/2, so 1/2 is exact.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pymatching
import stim

from .circuits import LeakyMemory
from .heralds import BOUNDARY, Edge, HeraldMap, herald_signatures

# Probability given to a herald edge that is absent from the static graph
# (only happens when the circuit noise p is 0). Large weight = "never use
# unless heralded".
_ABSENT_EDGE_P = 1e-12


def dem_of(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    return circuit.detector_error_model(
        decompose_errors=True, approximate_disjoint_errors=True
    )


def edge_weight(p: float) -> float:
    """MWPM edge weight: minus log-likelihood ratio log((1-p)/p)."""
    return math.log((1 - p) / p)


class BaselineMatcher:
    name = "pymatching_static"

    def __init__(self, mem: LeakyMemory):
        self.mem = mem
        self.matching = pymatching.Matching.from_detector_error_model(dem_of(mem.circuit))

    def decode_batch(self, dets: np.ndarray, heralds: np.ndarray | None = None) -> np.ndarray:
        return self.matching.decode_batch(dets)[:, 0].astype(np.uint8)


class AdaptiveMatcher:
    name = "pymatching_adaptive"

    def __init__(self, mem: LeakyMemory, herald_map: HeraldMap | None = None):
        self.mem = mem
        self.herald_map = herald_map or herald_signatures(mem)
        clean = mem.with_leak_rate(0.0).circuit
        m = pymatching.Matching.from_detector_error_model(dem_of(clean))
        m.ensure_num_fault_ids(1)

        # Unique herald edges, their original data, and herald -> edge ids.
        self.edges: list[Edge] = self.herald_map.unique_edges()
        index = {e: i for i, e in enumerate(self.edges)}
        self.orig: list[tuple[set, float, float]] = []
        for e in self.edges:
            if e.v == BOUNDARY:
                if not m.has_boundary_edge(e.u):
                    m.add_boundary_edge(
                        e.u, fault_ids={0} if e.obs else set(),
                        weight=edge_weight(_ABSENT_EDGE_P), error_probability=_ABSENT_EDGE_P,
                    )
                d = m.get_boundary_edge_data(e.u)
            else:
                if not m.has_edge(e.u, e.v):
                    m.add_edge(
                        e.u, e.v, fault_ids={0} if e.obs else set(),
                        weight=edge_weight(_ABSENT_EDGE_P), error_probability=_ABSENT_EDGE_P,
                    )
                d = m.get_edge_data(e.u, e.v)
            self.orig.append((set(d["fault_ids"]), d["weight"], d["error_probability"]))
        self.herald_edges: list[np.ndarray] = [
            np.array([index[e] for e in es], dtype=np.int64) for es in self.herald_map.edges
        ]
        self.matching = m

    def _set(self, i: int, fault_ids: set, weight: float, p: float) -> None:
        e = self.edges[i]
        if e.v == BOUNDARY:
            self.matching.add_boundary_edge(
                e.u, fault_ids=fault_ids, weight=weight, error_probability=p,
                merge_strategy="replace",
            )
        else:
            self.matching.add_edge(
                e.u, e.v, fault_ids=fault_ids, weight=weight, error_probability=p,
                merge_strategy="replace",
            )

    def erased_edges(self, herald_row: np.ndarray) -> np.ndarray:
        fired = np.flatnonzero(herald_row)
        if len(fired) == 0:
            return np.zeros(0, dtype=np.int64)
        return np.unique(np.concatenate([self.herald_edges[k] for k in fired]))

    def decode_one(self, det_row: np.ndarray, herald_row: np.ndarray) -> int:
        ids = self.erased_edges(herald_row)
        for i in ids:
            self._set(int(i), self.orig[i][0], 0.0, 0.5)
        try:
            pred = int(self.matching.decode(det_row)[0])
        finally:
            for i in ids:
                fid, w, p = self.orig[i]
                self._set(int(i), fid, w, p)
        return pred

    def decode_batch(self, dets: np.ndarray, heralds: np.ndarray) -> np.ndarray:
        out = np.zeros(dets.shape[0], dtype=np.uint8)
        any_h = heralds.any(axis=1)
        quiet = np.flatnonzero(~any_h)
        if len(quiet):
            # No herald fired: static graph is exactly right, use the fast batch path.
            out[quiet] = self.matching.decode_batch(dets[quiet])[:, 0]
        for s in np.flatnonzero(any_h):
            out[s] = self.decode_one(dets[s], heralds[s])
        return out
