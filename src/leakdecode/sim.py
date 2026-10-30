"""Sampling + decoding loop, parallelised over processes."""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, asdict
from multiprocessing import get_context

import numpy as np

from .circuits import build_memory
from .decoders import AdaptiveMatcher, BaselineMatcher
from .stats import per_round, wilson_interval

DECODERS = ("pymatching_static", "pymatching_adaptive")


def sample(mem, shots: int, seed: int):
    """Return (detectors, observable0, heralds) as uint8/bool arrays."""
    sampler = mem.circuit.compile_sampler(seed=seed)
    meas = sampler.sample(shots)
    dets, obs = mem.circuit.compile_m2d_converter().convert(
        measurements=meas, separate_observables=True
    )
    heralds = meas[:, mem.herald_meas_index]
    return dets, obs[:, 0].astype(np.uint8), heralds


_CACHE: dict = {}


def _decoders_for(d, rounds, p, p_leak, leak_on):
    key = (d, rounds, p, p_leak, tuple(leak_on))
    if key not in _CACHE:
        mem = build_memory(d, rounds, p, p_leak, leak_on)
        _CACHE.clear()  # keep memory small: one point per worker at a time
        _CACHE[key] = (mem, BaselineMatcher(mem), AdaptiveMatcher(mem))
    return _CACHE[key]


def run_chunk(args):
    d, rounds, p, p_leak, leak_on, shots, seed = args
    mem, base, adapt = _decoders_for(d, rounds, p, p_leak, leak_on)
    dets, obs, her = sample(mem, shots, seed)
    t0 = time.perf_counter()
    pb = base.decode_batch(dets)
    t1 = time.perf_counter()
    pa = adapt.decode_batch(dets, her)
    t2 = time.perf_counter()
    return {
        "key": (d, p_leak),
        "shots": shots,
        "err_static": int((pb != obs).sum()),
        "err_adaptive": int((pa != obs).sum()),
        "heralded_shots": int(her.any(axis=1).sum()),
        "heralds_fired": int(her.sum()),
        "t_static": t1 - t0,
        "t_adaptive": t2 - t1,
    }


@dataclass
class PointResult:
    distance: int
    rounds: int
    p: float
    p_leak: float
    decoder: str
    shots: int
    errors: int
    p_shot: float
    p_shot_lo: float
    p_shot_hi: float
    ler_per_round: float
    ler_per_round_lo: float
    ler_per_round_hi: float
    mean_heralds_per_shot: float
    frac_shots_with_herald: float
    decode_us_per_shot_batch: float


def run_sweep(
    distances,
    p: float,
    leak_rates,
    leak_on,
    max_shots: dict,
    max_errors: int,
    chunk: int,
    workers: int,
    seed: int = 12345,
    log=print,
) -> list[PointResult]:
    """Sample every (d, p_leak) point until both decoders saw ``max_errors``
    logical errors or ``max_shots[d]`` shots were used."""
    points = [(d, pl) for d in distances for pl in leak_rates]
    acc = {pt: dict(shots=0, err_static=0, err_adaptive=0, heralded_shots=0,
                    heralds_fired=0, t_static=0.0, t_adaptive=0.0) for pt in points}
    next_seed = seed
    ctx = get_context("fork")
    with ctx.Pool(workers) as pool:
        while True:
            todo = []
            for (d, pl) in points:
                a = acc[(d, pl)]
                done = a["shots"] >= max_shots[d] or min(a["err_static"], a["err_adaptive"]) >= max_errors
                if done:
                    continue
                remaining = max_shots[d] - a["shots"]
                n_chunks = max(1, min(workers, -(-remaining // chunk)))
                for _ in range(n_chunks):
                    n = min(chunk, remaining)
                    if n <= 0:
                        break
                    remaining -= n
                    todo.append((d, d, p, pl, tuple(leak_on), n, next_seed))
                    next_seed += 1
            if not todo:
                break
            # group by point so each worker reuses its cached decoders
            for r in pool.imap_unordered(run_chunk, todo):
                a = acc[r["key"]]
                for k in a:
                    a[k] += r[k]
            log("  progress: " + ", ".join(
                f"d={d} pL={pl:g}: {acc[(d, pl)]['shots']} shots "
                f"({acc[(d, pl)]['err_static']}/{acc[(d, pl)]['err_adaptive']} err)"
                for (d, pl) in points))
    out = []
    for (d, pl) in points:
        a = acc[(d, pl)]
        for dec, ek, tk in (("pymatching_static", "err_static", "t_static"),
                            ("pymatching_adaptive", "err_adaptive", "t_adaptive")):
            n, k = a["shots"], a[ek]
            lo, hi = wilson_interval(k, n)
            out.append(PointResult(
                distance=d, rounds=d, p=p, p_leak=pl, decoder=dec, shots=n, errors=k,
                p_shot=k / n, p_shot_lo=lo, p_shot_hi=hi,
                ler_per_round=per_round(k / n, d),
                ler_per_round_lo=per_round(lo, d), ler_per_round_hi=per_round(hi, d),
                mean_heralds_per_shot=a["heralds_fired"] / n,
                frac_shots_with_herald=a["heralded_shots"] / n,
                decode_us_per_shot_batch=1e6 * a[tk] / n,
            ))
    return out
