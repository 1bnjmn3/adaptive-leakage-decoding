"""All figures are drawn from the CSV files in results/ (never from memory)."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

COLORS = {3: "tab:blue", 5: "tab:orange", 7: "tab:green", 9: "tab:red", 11: "tab:purple"}
STYLE = {"pymatching_static": ("o", "-", "static (ignores heralds)"),
         "pymatching_adaptive": ("s", "--", "adaptive (uses heralds)")}


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def plot_ler_vs_leak(sweep_csv, out_png, p):
    rows = read_csv(sweep_csv)
    groups = defaultdict(list)
    for r in rows:
        groups[(int(r["distance"]), r["decoder"])].append(r)
    fig, ax = plt.subplots(figsize=(7, 5))
    for (d, dec), rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: float(r["p_leak"]))
        x = np.array([float(r["p_leak"]) for r in rs])
        y = np.array([float(r["ler_per_round"]) for r in rs])
        lo = np.array([float(r["ler_per_round_lo"]) for r in rs])
        hi = np.array([float(r["ler_per_round_hi"]) for r in rs])
        ok = y > 0
        m, ls, lab = STYLE[dec]
        # p_leak = 0 cannot be shown on a log axis: draw it at a fixed small x
        xplot = np.where(x > 0, x, x[x > 0].min() / 3 if (x > 0).any() else 1e-4)
        ax.errorbar(xplot[ok], y[ok], yerr=[y[ok] - lo[ok], hi[ok] - y[ok]], marker=m, ls=ls,
                    color=COLORS.get(d, "k"), capsize=2, label=f"d={d}, {lab}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("heralded-erasure (leakage) rate p_L per location per round\n(leftmost points: p_L = 0)")
    ax.set_ylabel("logical error rate per round")
    ax.set_title(f"Rotated surface-code memory, rounds = d, circuit noise p = {p:g}\n95% Wilson intervals")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)


def plot_ler_vs_distance(sweep_csv, fits_csv, out_png):
    rows = read_csv(sweep_csv)
    fits = {(float(r["p_leak"]), r["decoder"]): r for r in read_csv(fits_csv)}
    leaks = sorted({float(r["p_leak"]) for r in rows})
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")
    for i, pl in enumerate(leaks):
        c = cmap(i / max(1, len(leaks) - 1))
        for dec in STYLE:
            rs = sorted([r for r in rows if float(r["p_leak"]) == pl and r["decoder"] == dec],
                        key=lambda r: int(r["distance"]))
            d = np.array([int(r["distance"]) for r in rs])
            y = np.array([float(r["ler_per_round"]) for r in rs])
            lo = np.array([float(r["ler_per_round_lo"]) for r in rs])
            hi = np.array([float(r["ler_per_round_hi"]) for r in rs])
            ok = y > 0
            m, ls, lab = STYLE[dec]
            ax.errorbar(d[ok], y[ok], yerr=[y[ok] - lo[ok], hi[ok] - y[ok]], marker=m, ls="none",
                        color=c, capsize=2, label=f"p_L={pl:g} {dec.split('_')[1]}")
            f = fits.get((pl, dec))
            if f and f["a"] not in ("", "nan"):
                dd = np.linspace(d.min(), d.max(), 20)
                ax.plot(dd, 10 ** (float(f["a"]) + float(f["b"]) * dd), ls=ls, color=c, lw=1)
    ax.set_yscale("log"); ax.set_xlabel("code distance d"); ax.set_ylabel("logical error rate per round")
    ax.set_title("Suppression with distance\n(lines: weighted fit log10 eps = a + b d; solid static, dashed adaptive)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)


def plot_qubit_savings(savings_csv, out_png):
    rows = read_csv(savings_csv)
    targets = sorted({(r["target_kind"], float(r["target"])) for r in rows}, key=lambda x: (x[0] != "window", -x[1]))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for kind, t in targets:
        rs = sorted([r for r in rows if r["target_kind"] == kind and float(r["target"]) == t
                     and r["qubit_ratio"] not in ("", "nan", "inf")],
                    key=lambda r: float(r["p_leak"]))
        if not rs:
            continue
        ax.plot([float(r["p_leak"]) for r in rs], [float(r["qubit_ratio"]) for r in rs], marker="o",
                ls="-" if kind == "window" else ":",
                label=(f"target {t:g} per d-round window (LCD paper's)" if kind == "window" else f"target {t:g} per round"))
        if kind == "window":
            odd = [r for r in rs if r.get("qubit_ratio_odd") not in (None, "", "nan", "inf")]
            ax.plot([float(r["p_leak"]) for r in odd], [float(r["qubit_ratio_odd"]) for r in odd], marker="s",
                    mfc="none", ls="none", color="k", label=f"same window target, odd d (real patch sizes)")
    ax.axhline(1, color="k", lw=0.8)
    ax.set_xscale("log"); ax.set_xlabel("heralded-erasure rate p_L")
    ax.set_ylabel("qubits(static) / qubits(adaptive)")
    ax.set_title("Qubit saving from using heralds (EXTRAPOLATED from d = 3, 5, 7 fits)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)


def plot_timing(timing_csv, out_png, hw_csv=None):
    rows = read_csv(timing_csv)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    series = defaultdict(list)
    for r in rows:
        series[f"{r['decoder']}: {r['method']}"].append((int(r["distance"]), float(r["us_per_round"])))
    if hw_csv and Path(hw_csv).exists():
        for r in read_csv(hw_csv):
            series[f"C++ union-find ({r['mode']})"].append((int(r["distance"]), float(r["us_per_round"])))
    for lab, pts in series.items():
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", label=lab)
    ax.axhline(1.0, color="r", ls=":", label="1 us/round target (LCD paper)")
    ax.set_yscale("log"); ax.set_xlabel("code distance d"); ax.set_ylabel("decode time per round (us), 1 CPU core")
    ax.set_title("Decoding time per round (software on this machine, not hardware)")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)


def plot_hw_scaling(hw_csv, out_png):
    rows = read_csv(hw_csv)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for mode in sorted({r["mode"] for r in rows}):
        rs = sorted([r for r in rows if r["mode"] == mode], key=lambda r: int(r["distance"]))
        d = [int(r["distance"]) for r in rs]
        axs[0].plot(d, [float(r["us_per_round"]) for r in rs], marker="o", label=f"mean ({mode})")
        axs[0].plot(d, [float(r["max_us_per_round"]) for r in rs], marker="x", ls=":", label=f"worst shot ({mode})")
        axs[1].plot(d, [int(r["static_bytes"]) / 1024 for r in rs], marker="o", label=mode)
    axs[0].axhline(1.0, color="r", ls=":", label="1 us/round")
    axs[0].set_yscale("log"); axs[0].set_xlabel("d"); axs[0].set_ylabel("us per round"); axs[0].legend(fontsize=7)
    axs[0].set_title("C++ union-find on one CPU core (not FPGA)"); axs[0].grid(True, which="both", alpha=0.3)
    axs[1].set_xlabel("d"); axs[1].set_ylabel("decoder state (KiB)"); axs[1].set_title("Static memory (whole d-round window)")
    axs[1].grid(True, alpha=0.3); axs[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)
