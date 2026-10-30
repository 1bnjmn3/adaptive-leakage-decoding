"""Hardware add-on stage: build and run the C++ union-find decoder (hw/).

For each distance we sample one batch of shots, export it in the flat format of
hw/INTERFACE.md, and decode it four ways on exactly the same syndromes:
PyMatching static, PyMatching adaptive, C++ UF ignoring heralds, C++ UF with
heralded edges pre-grown as erasures. We report logical error rates, how often
UF agrees with PyMatching, and the UF's time and memory per round.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .circuits import build_memory
from .decoders import AdaptiveMatcher, BaselineMatcher
from .export import write_problem
from .sim import sample
from .stats import per_round, wilson_interval

REPO = Path(__file__).resolve().parents[2]
HW = REPO / "hw"


def _binary() -> Path:
    for cand in (HW / "build" / "uf_decode", HW / "uf_decode"):
        if cand.exists():
            return cand
    raise FileNotFoundError("uf_decode binary not found")


def build_hw(log=print) -> Path:
    """Build hw/ with make if present, else CMake+Ninja. Returns the binary path."""
    cxx = os.environ.get("CXX") or str(HW / "tools" / "cxx")
    if shutil.which("make"):
        cmd = ["make", "-C", str(HW), f"CXX={cxx}"]
    else:
        exe_dir = Path(os.sys.executable).parent
        cmake = shutil.which("cmake") or str(exe_dir / "cmake")
        ninja = shutil.which("ninja") or str(exe_dir / "ninja")
        bdir = HW / "build"
        subprocess.run([cmake, "-S", str(HW), "-B", str(bdir), "-G", "Ninja",
                        f"-DCMAKE_MAKE_PROGRAM={ninja}", f"-DCMAKE_CXX_COMPILER={cxx}",
                        "-DCMAKE_BUILD_TYPE=Release"], check=True, capture_output=True, text=True)
        cmd = [cmake, "--build", str(bdir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"hw build failed: {r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    log("  built C++ union-find decoder")
    return _binary()


def _run(binary, d: Path, erasures: bool, repeat: int) -> tuple[np.ndarray, dict]:
    cmd = [str(binary), str(d / "graph.txt"), str(d / "syndromes.u8"), str(d / "pred.u8")]
    if erasures:
        cmd += ["--erasures", str(d / "erasures.u8")]
    cmd += ["--repeat", str(repeat)]
    r = subprocess.run(cmd, capture_output=True, text=True, check=True)
    line = [ln for ln in r.stdout.splitlines() if ln.strip().startswith("{")][-1]
    return np.fromfile(d / "pred.u8", dtype=np.uint8), json.loads(line)


def run_hw(distances, p, p_leak, leak_on, shots, out: Path, log=print, repeat=5):
    binary = build_hw(log)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for dist in distances:
            mem = build_memory(dist, dist, p, p_leak, leak_on)
            base, adapt = BaselineMatcher(mem), AdaptiveMatcher(mem)
            dets, obs, her = sample(mem, shots, 777 + dist)
            pdir = write_problem(adapt, dets, her, Path(tmp) / f"d{dist}")
            p_static = base.decode_batch(dets)
            p_adapt = adapt.decode_batch(dets, her)
            for mode, use_er, ref_name, ref in (
                ("ignore heralds", False, "pymatching_static", p_static),
                ("erasure-aware", True, "pymatching_adaptive", p_adapt),
            ):
                pred, info = _run(binary, pdir, use_er, repeat)
                assert len(pred) == shots
                k = int((pred != obs).sum())
                lo, hi = wilson_interval(k, shots)
                k_ref = int((ref != obs).sum())
                rows.append(dict(
                    distance=dist, rounds=dist, p=p, p_leak=p_leak, mode=mode, shots=shots,
                    uf_errors=k, uf_ler_per_round=per_round(k / shots, dist),
                    uf_ler_per_round_lo=per_round(lo, dist), uf_ler_per_round_hi=per_round(hi, dist),
                    reference=ref_name, ref_errors=k_ref, ref_ler_per_round=per_round(k_ref / shots, dist),
                    agreement_with_reference=float((pred == ref).mean()),
                    num_detectors=info["num_detectors"], num_edges=info["num_edges"],
                    us_per_shot=info["ns_per_shot"] / 1000, us_per_round=info["ns_per_shot"] / 1000 / dist,
                    max_us_per_round=info["max_ns_per_shot"] / 1000 / dist,
                    static_bytes=info["static_bytes"],
                ))
                log(f"  hw d={dist} {mode}: UF errors {k}/{shots} vs {ref_name} {k_ref}; "
                    f"agree {rows[-1]['agreement_with_reference']:.4f}; {rows[-1]['us_per_round']:.3f} us/round")
    path = out / "hw_uf.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return path
