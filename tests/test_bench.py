"""End-to-end quick benchmark: CSVs and figures, without the C++ decoder."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from leakdecode.bench import main

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
PNG_FILES = (
    "ler_vs_leak.png",
    "ler_vs_distance.png",
    "qubit_savings.png",
    "decode_timing.png",
)


def _rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        return fieldnames, list(reader)


@pytest.mark.slow
def test_quick_bench_writes_tables_and_figures(tmp_path: Path):
    out = tmp_path / "r"
    rc = main(["--quick", "--out", str(out), "--no-hw", "--workers", "2"])
    assert rc == 0

    sweep_fields, sweep = _rows(out / "ler_sweep.csv")
    for name in (
        "distance",
        "rounds",
        "p",
        "p_leak",
        "decoder",
        "shots",
        "errors",
        "ler_per_round",
        "ler_per_round_lo",
        "ler_per_round_hi",
    ):
        assert name in sweep_fields
    assert sweep
    assert {row["decoder"] for row in sweep} == {"pymatching_static", "pymatching_adaptive"}
    assert {int(row["distance"]) for row in sweep} == {3, 5}

    fit_fields, fits = _rows(out / "fits.csv")
    for name in ("p_leak", "decoder", "a", "b", "Lambda"):
        assert name in fit_fields
    assert fits
    assert {row["decoder"] for row in fits} == {"pymatching_static", "pymatching_adaptive"}

    save_fields, savings = _rows(out / "qubit_savings.csv")
    for name in ("p_leak", "target_kind", "target", "qubit_ratio", "qubit_ratio_odd", "d_static", "d_adaptive"):
        assert name in save_fields
    assert savings
    assert {"per_round", "window"} <= {row["target_kind"] for row in savings}

    time_fields, timing = _rows(out / "decode_timing.csv")
    for name in ("distance", "decoder", "method", "us_per_shot", "us_per_round"):
        assert name in time_fields
    assert timing
    assert {row["decoder"] for row in timing} >= {"pymatching_static", "pymatching_adaptive"}

    for name in PNG_FILES:
        blob = (out / name).read_bytes()
        assert blob.startswith(PNG_MAGIC)
        assert len(blob) > len(PNG_MAGIC)
