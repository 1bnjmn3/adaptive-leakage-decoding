"""``leakdecode-bench``: regenerate every results CSV and plot from scratch.

    leakdecode-bench            # full run (minutes)
    leakdecode-bench --quick    # small CI run (seconds)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import shutil
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pymatching
import stim

from . import plots
from .analysis import fits_from_sweep, savings_from_fits
from .sim import run_sweep
from .timing import time_decoders

P_DEFAULT = 2e-3
FULL = dict(distances=(3, 5, 7), leak_rates=(0.0, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2),
            max_shots={3: 1_000_000, 5: 600_000, 7: 600_000}, max_errors=500, chunk=25_000,
            timing_shots=4000, hw_shots=20_000, hw_distances=(3, 5, 7, 9, 11))
QUICK = dict(distances=(3, 5), leak_rates=(0.0, 5e-3, 2e-2),
             max_shots={3: 20_000, 5: 10_000}, max_errors=50, chunk=5_000,
             timing_shots=300, hw_shots=2_000, hw_distances=(3, 5))


def write_csv(path, rows):
    rows = list(rows)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="small, fast run for CI")
    ap.add_argument("--out", default="results", help="output directory (deleted and recreated)")
    ap.add_argument("--p", type=float, default=P_DEFAULT, help="circuit-level depolarising noise strength")
    ap.add_argument("--leak-on", default="data,ancilla", help="comma list from {data,ancilla}")
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--no-hw", action="store_true", help="skip the C++ union-find stage")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args(argv)
    cfg = QUICK if args.quick else FULL
    leak_on = tuple(s for s in args.leak_on.split(",") if s)
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    t_start = time.time()
    log_f = open(out / "bench_log.txt", "w")

    def log(msg):
        print(msg, flush=True)
        log_f.write(msg + "\n"); log_f.flush()

    log(f"leakdecode-bench ({'quick' if args.quick else 'full'}) p={args.p} leak_on={leak_on} workers={args.workers}")

    # 1. logical error rate sweep
    t0 = time.time()
    res = run_sweep(cfg["distances"], args.p, cfg["leak_rates"], leak_on, cfg["max_shots"],
                    cfg["max_errors"], cfg["chunk"], args.workers, seed=args.seed, log=log)
    sweep_rows = [asdict(r) for r in res]
    write_csv(out / "ler_sweep.csv", sweep_rows)
    t_sweep = time.time() - t0
    log(f"sweep done in {t_sweep:.1f}s")

    # 2. fits + qubit savings (string-typed like the CSV, so plots use the same path)
    str_rows = [{k: str(v) for k, v in r.items()} for r in sweep_rows]
    fits = fits_from_sweep(str_rows)
    write_csv(out / "fits.csv", fits)
    write_csv(out / "qubit_savings.csv", savings_from_fits(fits))

    # 3. timing
    t0 = time.time()
    tim = time_decoders(cfg["distances"], args.p, 5e-3, leak_on, cfg["timing_shots"])
    write_csv(out / "decode_timing.csv", tim)
    log(f"timing done in {time.time() - t0:.1f}s")

    # 4. C++ union-find (hardware add-on)
    hw_csv = None
    hw_status = "skipped (--no-hw)"
    if not args.no_hw:
        try:
            from .hwbench import run_hw
            hw_csv = run_hw(cfg["hw_distances"], args.p, 5e-3, leak_on, cfg["hw_shots"], out, log)
            hw_status = "ok"
        except Exception as e:  # report, don't hide
            hw_status = f"failed: {type(e).__name__}: {e}"
            log(f"HW stage {hw_status}")

    # 5. plots (from the CSVs on disk)
    plots.plot_ler_vs_leak(out / "ler_sweep.csv", out / "ler_vs_leak.png", args.p)
    plots.plot_ler_vs_distance(out / "ler_sweep.csv", out / "fits.csv", out / "ler_vs_distance.png")
    plots.plot_qubit_savings(out / "qubit_savings.csv", out / "qubit_savings.png")
    plots.plot_timing(out / "decode_timing.csv", out / "decode_timing.png", hw_csv)
    if hw_csv:
        plots.plot_hw_scaling(hw_csv, out / "hw_uf_scaling.png")

    total = time.time() - t_start
    info = dict(mode="quick" if args.quick else "full", p=args.p, leak_on=leak_on, config={k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()}) for k, v in cfg.items()},
                workers=args.workers, seed=args.seed, sweep_seconds=t_sweep, total_seconds=total,
                hw_status=hw_status, python=sys.version.split()[0], stim=stim.__version__,
                pymatching=pymatching.__version__, numpy=np.__version__, machine=platform.machine(),
                cpu_count=os.cpu_count(), finished=time.strftime("%Y-%m-%d %H:%M:%S %Z"))
    (out / "run_info.json").write_text(json.dumps(info, indent=2))
    log(f"TOTAL runtime {total:.1f}s; outputs in {out}/")
    log_f.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
