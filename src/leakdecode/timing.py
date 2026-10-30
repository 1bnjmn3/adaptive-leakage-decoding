"""Single-threaded decoding-time measurements (microseconds per round)."""

from __future__ import annotations

import time

import numpy as np

from .circuits import build_memory
from .decoders import AdaptiveMatcher, BaselineMatcher
from .sim import sample


def time_decoders(distances, p, p_leak, leak_on, shots, seed=999, repeats=3):
    rows = []
    for d in distances:
        mem = build_memory(d, d, p, p_leak, leak_on)
        base, adapt = BaselineMatcher(mem), AdaptiveMatcher(mem)
        dets, obs, her = sample(mem, shots, seed + d)
        # warm-up (builds PyMatching's internal graph)
        base.decode_batch(dets[:10]); adapt.decode_batch(dets[:10], her[:10])

        def best(fn):
            ts = []
            for _ in range(repeats):
                t0 = time.perf_counter(); fn(); ts.append(time.perf_counter() - t0)
            return min(ts)

        t_batch = best(lambda: base.decode_batch(dets))
        t_loop = best(lambda: [base.matching.decode(dets[s]) for s in range(shots)])
        t_adapt = best(lambda: adapt.decode_batch(dets, her))
        for name, t, how in (
            ("pymatching_static", t_batch, "decode_batch (C++ loop)"),
            ("pymatching_static", t_loop, "python loop of decode()"),
            ("pymatching_adaptive", t_adapt, "python loop, reweight+decode+restore"),
        ):
            rows.append(dict(distance=d, rounds=d, p=p, p_leak=p_leak, decoder=name, method=how,
                             shots=shots, us_per_shot=1e6 * t / shots,
                             us_per_round=1e6 * t / shots / d,
                             frac_shots_with_herald=float(her.any(axis=1).mean())))
    return rows
