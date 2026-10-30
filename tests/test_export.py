"""Flat graph export read by the C++ union-find decoder."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from leakdecode.circuits import build_memory
from leakdecode.decoders import AdaptiveMatcher
from leakdecode.export import graph_edges, write_problem
from leakdecode.heralds import BOUNDARY


def _normalised_ends(u: int, v: int) -> tuple[int, int]:
    if v == BOUNDARY or u < v:
        return u, v
    return v, u


def test_graph_file_lists_every_herald_edge(tmp_path: Path):
    mem = build_memory(3, 3, 2e-3, 1e-2)
    adaptive = AdaptiveMatcher(mem)
    num_detectors = mem.circuit.num_detectors
    shots = 2
    dets = np.zeros((shots, num_detectors), dtype=np.uint8)
    heralds = np.zeros((shots, mem.num_heralds), dtype=np.uint8)
    heralds[1] = 1

    outdir = write_problem(adaptive, dets, heralds, tmp_path / "problem")
    rows, herald_to_row = graph_edges(adaptive)

    lines = (outdir / "graph.txt").read_text().splitlines()
    header_detectors, header_edges = (int(piece) for piece in lines[0].split())
    assert header_detectors == num_detectors
    assert header_edges == len(rows) == len(lines) - 1
    assert len(lines[0].split()) == 2

    parsed = []
    for line in lines[1:]:
        u, v, weight, obs = (int(piece) for piece in line.split())
        assert 0 <= u < num_detectors
        assert v == BOUNDARY or 0 <= v < num_detectors
        assert weight >= 0
        assert obs in (0, 1)
        parsed.append((u, v, weight, obs))
    assert parsed == [(a, b, w, o) for a, b, w, o in rows]

    assert len(herald_to_row) == len(adaptive.edges)
    for edge, row_index in zip(adaptive.edges, herald_to_row):
        row_index = int(row_index)
        assert 0 <= row_index < len(parsed)
        u, v, _weight, obs = parsed[row_index]
        assert (u, v) == _normalised_ends(edge.u, edge.v)
        assert obs == edge.obs

    erasures = np.fromfile(outdir / "erasures.u8", dtype=np.uint8).reshape(shots, header_edges)
    syndromes = np.fromfile(outdir / "syndromes.u8", dtype=np.uint8)
    assert syndromes.size == shots * num_detectors
    assert int(erasures[0].sum()) == 0
    assert int(erasures[1, herald_to_row].min()) == 1
