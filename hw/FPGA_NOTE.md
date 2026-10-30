# Mapping the union-find decoder onto an FPGA (one-page sketch)

This note sketches how the C++ decoder in `hw/src/` could map onto an FPGA. **No RTL was written, nothing was synthesised,
and no FPGA latency is claimed.** All the measured numbers below are CPU numbers from `results/hw_uf.csv`, taken on one
core of a cloud x86-64 VM.

## What the C++ decoder does

It is an unweighted Delfosse–Zémor union-find decoder over the detector graph of one d-round window:

* Heralded edges (the `erasures` input) start fully grown, and their endpoints are merged before any growth.
* Each growth step adds one half-edge to every boundary edge of every *odd* cluster. An edge counter goes 0 → 1 → 2,
  and at 2 the edge is fully grown and its two clusters merge. Merging is union by size, the parities are XORed and the
  boundary flags are ORed.
* Peeling builds a spanning forest of the fully grown edges, rooted first at the boundary. It processes leaves first and
  flips the edges the syndrome needs. Edges with `obs = 1` toggle the predicted logical bit.
* Every array is allocated once in the constructor, and `decode()` does no allocation. Only touched nodes and edges are
  reset after each shot.

## State and memory

Per node the decoder keeps parent, size, parity, a boundary flag, a defect bit and peeling bookkeeping. Per edge it keeps
a 2-bit growth counter. It also keeps read-only tables for the endpoints, the observable bit and the adjacency.
The C++ `static_bytes` counts all of it, using 32-bit indices:

| d | detectors in window | edges | C++ decoder state |
|---|---|---|---|
| 3 | 24 | 78 | 3.4 KiB |
| 5 | 120 | 502 | 20.0 KiB |
| 7 | 336 | 1558 | 60.5 KiB |
| 9 | 720 | 3534 | 135.2 KiB |
| 11 | 1320 | 6718 | 254.9 KiB |

The state grows as about d³ (d² detectors per round × d rounds). On an FPGA, the indices would shrink to
⌈log₂ N⌉ bits and the growth counter is 2 bits, so the dynamic state is a few bits per node and per edge. That fits in
on-chip BRAM or registers at these sizes.

## Mapping

* **Heralds.** A ROM maps each herald bit to its edge set; `heralds.py` computes these edges offline. A fired herald sets
  those growth counters to 2 before the first growth step. This costs one cycle and no extra logic in the growth loop.
* **Growth.** One processing element per node (or per small tile of nodes) increments its odd cluster's boundary edges in
  parallel, so a growth step takes a constant number of cycles. Merges are the hard part, because many merges into one root
  are serial. Helios (Liyanage, Wu, Tagare, Zhong, arXiv:2406.08491) distributes union-find across processing elements on
  an FPGA. The LCD paper (arXiv:2411.10343) reports an adaptive clustering decoder under 1 µs per round up to d = 17 on an
  AMD VU19P at 285 MHz. Neither is a measurement of this design.
* **Streaming.** The C++ decodes a whole d-round block. A real decoder would slide a window, committing the oldest
  rounds while new ones arrive. That is not implemented here.

## Timing against the 1 µs-per-round target (CPU, not FPGA)

At p_L = 0.005, with 20,000 shots and the per-shot minimum of 5 repeats, using the erasure-aware mode:

| d | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|
| mean µs/round | 0.07 | 0.23 | 0.47 | 0.80 | 1.21 |
| worst shot µs/round | 1.1 | 2.4 | 3.5 | 5.4 | 5.6 |

The mean stays under 1 µs per round up to d = 9 and goes above it at d = 11. The worst shot misses the budget from
d = 5 up. A real-time decoder has to bound the *worst* case, because a slow shot builds a backlog.

## Risks

* The number of growth steps is bounded only by the window diameter, so a rare large cluster sets the tail latency.
* Unweighted growth gives up accuracy compared with matching (see the README table).
* `HERALDED_ERASE` is a one-location, perfectly flagged event. Real leakage persists and spreads, so it would need herald
  edges held across several rounds and a model of late detection. This sketch does not cover that.
