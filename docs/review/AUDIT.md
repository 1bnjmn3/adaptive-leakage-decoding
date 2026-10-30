# Independent audit, 30 Oct 2026

Checker's log. I did not rely on the build report or the docs; every check below was rerun or recomputed.

## Reproduction
* Copied `results/` to a backup, deleted it, ran `leakdecode-bench` (exit 0, 477 s total, sweep 345.5 s, 6 workers).
  All 12 CSVs and plots were regenerated.
* `ler_sweep.csv`: every column is identical except `decode_us_per_shot_batch` (a wall-clock timer).
  `fits.csv` and `qubit_savings.csv` are byte-identical. `hw_uf.csv`: every error count, reference count, agreement,
  detector/edge count and `static_bytes` is identical; only the three timing columns move (by at most about 5%).
* `decode_timing.csv` moved by up to about 6% (e.g. d = 7 adaptive 93 → 90 µs/round, batch 0.53 → 0.50).
* Because the Monte Carlo output reproduced exactly and the docs quote the original run's timings, I put the original
  `results/` back unchanged. The rerun's timing evidence is kept in `docs/review/audit_rerun_2026-10-30/`.
* `hw/` clean CMake build in a fresh directory: ctest 1/1 passed. `pytest -q`: 70 passed.

## Numbers checked against the CSVs
README d = 7 table (rates, error counts, ratios, intervals), Λ table, qubit-saving table (continuous, odd d, 10⁻⁹),
C++ table (errors, agreement, KiB), timing table, minimum error count 159; LEARN's quoted rates, Λ, savings and DEM
values (recomputed with Stim); FPGA_NOTE detector/edge/KiB and timing tables. All match, except for the items fixed below.

## Fixes
1. `LEARN.md` §5 and `docs/learn/01_physics.md`: the DEM was shown with `error(0.00505015) D1 D5`. The real DEM line is
   `error(0.00478582) D1 D5`. 0.00505015 is the edge PyMatching stores after it folds in the D1–D5 part of
   `D1 D5 ^ D4` (0.00478582 ⊕ 0.00026688). I corrected the line and said where 0.00505015 comes from.
2. `LEARN.md` §4 and `docs/learn/01_physics.md`: the erasure criterion ("leaves some representative … partly outside")
   was imprecise. It now says: correctable exactly when no nontrivial logical error is supported inside the erased set.
3. `LEARN.md` code walkthrough and `docs/learn/02_code_walkthrough.md`: the text said the C++ sources are not in the tree.
   They are in `hw/src/`. Corrected.
4. Same two files: the text said `plot_qubit_savings` does not draw `qubit_ratio_odd`. `plots.py:95-97` draws it as
   hollow squares. Corrected. Also `LEARN.md` §9 and `01_physics.md`: "the plot uses continuous d" now mentions the
   squares.
5. `LEARN.md` Q&A and `README.md` "How to run": CI was described as running. The workflow exists but has never run
   (no git repository or remote). Reworded as configured but not yet run.
6. `docs/LCD_PAPER_NOTES.md`: the Table 4 / Table 1 / IV.6 references and the 0.930 / 1.127 µs values are from arXiv v1.
   Added that v2 numbers it as Table 3 (0.676 / 0.820 µs) with the noise model in Table 4, IV.7.

None of these changes affects results, so the bench was not rerun again.

## Not changed (recommendations)
* README quotes 490 s / 8.2 min for the documented run; my rerun took 477 s. Fine as is; "~8 min" would be safer.
* README rounds 0.99985 to 99.98% (d = 11); 99.99% is the half-up value. Trivial.
* C++ timing is the per-shot minimum of 5 repeats with warm caches. The "worst shot" is the maximum of those minimums,
  and per-round time is block time divided by d (no streaming). This is optimistic for tail latency; say so if asked.
* The CI workflow uses `make`, which was never tested locally (no make on the box). Action tags v7 exist.
* LEARN's "agree to within about 10–20%" for the savings curves is borderline: at p_L = 2×10⁻³ the odd-d window
  ratio (1.00) and the continuous one (1.23) differ by about 19–23%.
