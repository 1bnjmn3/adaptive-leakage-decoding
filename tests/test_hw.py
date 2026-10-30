"""The C++ union-find decoder (hw/) against PyMatching on the same shots.

Skipped if no C++ toolchain is available to build hw/.
"""

from __future__ import annotations

import subprocess

import numpy as np
import pytest

from leakdecode import build_memory
from leakdecode.decoders import AdaptiveMatcher, BaselineMatcher
from leakdecode.export import write_problem
from leakdecode.sim import sample


@pytest.fixture(scope="module")
def uf_binary():
    from leakdecode.hwbench import build_hw
    try:
        return build_hw(log=lambda *_: None)
    except Exception as e:  # no compiler / cmake
        pytest.skip(f"cannot build hw/: {e}")


def _run(binary, pdir, erasures):
    cmd = [str(binary), str(pdir / "graph.txt"), str(pdir / "syndromes.u8"), str(pdir / "pred.u8")]
    if erasures:
        cmd += ["--erasures", str(pdir / "erasures.u8")]
    subprocess.run(cmd, check=True, capture_output=True)
    return np.fromfile(pdir / "pred.u8", dtype=np.uint8)


def test_cpp_unit_tests_pass(uf_binary):
    t = uf_binary.parent / "test_uf"
    if not t.exists():
        pytest.skip("test_uf not built")
    assert subprocess.run([str(t)], capture_output=True).returncode == 0


@pytest.mark.parametrize("distance", [3, 5])
def test_uf_agrees_with_pymatching(uf_binary, tmp_path, distance):
    mem = build_memory(distance, distance, 2e-3, 1e-2)
    base, adapt = BaselineMatcher(mem), AdaptiveMatcher(mem)
    dets, obs, her = sample(mem, 3000, 99)
    pdir = write_problem(adapt, dets, her, tmp_path)
    uf_static = _run(uf_binary, pdir, erasures=False)
    uf_erase = _run(uf_binary, pdir, erasures=True)
    assert (uf_static == base.decode_batch(dets)).mean() > 0.98
    assert (uf_erase == adapt.decode_batch(dets, her)).mean() > 0.98
    # using the heralds must not make the union-find worse
    assert (uf_erase != obs).sum() <= (uf_static != obs).sum()


def test_uf_erasure_only_noise_is_decoded_well(uf_binary, tmp_path):
    """No Pauli noise, only heralded erasures: the peeling decoder should almost never fail."""
    mem = build_memory(5, 5, 0.0, 2e-2)
    adapt = AdaptiveMatcher(mem)
    dets, obs, her = sample(mem, 2000, 5)
    pdir = write_problem(adapt, dets, her, tmp_path)
    pred = _run(uf_binary, pdir, erasures=True)
    assert (pred != obs).mean() < 0.01
    assert not np.any(pred[~dets.any(axis=1)])  # quiet shots predict no flip
