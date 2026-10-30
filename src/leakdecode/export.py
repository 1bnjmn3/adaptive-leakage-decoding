"""Export a decoding problem in the flat file format read by the C++ decoder.

See hw/INTERFACE.md. The graph is the adaptive decoder's static graph (DEM of
the circuit with p_leak = 0, plus every herald edge), so the same file serves
both "UF ignoring heralds" (no erasure file) and "UF with erasures".
Weights are integers: round(WEIGHT_SCALE * log((1-p)/p)), clipped at >= 0.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .decoders import AdaptiveMatcher
from .heralds import BOUNDARY

WEIGHT_SCALE = 10


def graph_edges(adapt: AdaptiveMatcher):
    """List of (u, v, weight_int, obs) for the static adaptive graph, and a map
    from herald-edge index (AdaptiveMatcher.edges) to exported edge index."""
    rows = []
    index = {}
    for u, v, data in adapt.matching.edges():
        vv = BOUNDARY if v is None else v
        a, b = (u, vv) if vv == BOUNDARY or u < vv else (vv, u)
        p = data.get("error_probability", -1.0)
        if p is None or not (0 < p < 1):
            p = math.exp(-data["weight"]) / (1 + math.exp(-data["weight"]))
        w = max(0, int(round(WEIGHT_SCALE * math.log((1 - p) / p))))
        obs = 1 if 0 in set(data.get("fault_ids", set())) else 0
        index[(a, b)] = len(rows)
        rows.append((a, b, w, obs))
    herald_to_row = np.array(
        [index[(e.u, e.v) if e.v == BOUNDARY or e.u < e.v else (e.v, e.u)] for e in adapt.edges],
        dtype=np.int64,
    )
    return rows, herald_to_row


def write_problem(adapt: AdaptiveMatcher, dets: np.ndarray, heralds: np.ndarray, outdir: Path):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows, herald_to_row = graph_edges(adapt)
    nd = adapt.mem.circuit.num_detectors
    with open(outdir / "graph.txt", "w") as f:
        f.write(f"{nd} {len(rows)}\n")
        for a, b, w, o in rows:
            f.write(f"{a} {b} {w} {o}\n")
    np.ascontiguousarray(dets.astype(np.uint8)).tofile(outdir / "syndromes.u8")
    er = np.zeros((dets.shape[0], len(rows)), dtype=np.uint8)
    for s in np.flatnonzero(heralds.any(axis=1)):
        ids = adapt.erased_edges(heralds[s])
        er[s, herald_to_row[ids]] = 1
    er.tofile(outdir / "erasures.u8")
    return outdir
