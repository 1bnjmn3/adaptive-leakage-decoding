# Physics and statistics review

> **Follow-up (after this review):** this review was written before the final full run, so its note that `results/` is absent is out of date. Changes made in response:
> * The window-target solver now bisects only on the decreasing branch beyond the peak of d·ε(d), so it returns the larger root (`analysis._distance`). There is a test for the reviewer's two-root example, and the full-run savings are unchanged.
> * The odd-d qubit ratio for the window target is now drawn on `qubit_savings.png` and given in the README.
> * The `BaselineMatcher` docstring now describes Stim's disjoint-error DEM correctly.
> * The fixed-n Wilson assumption, the lack of uncertainty on the fit, and the "herald bits inside MWPM" framing are listed as limits in the README.
>
> None of the runs in `results/` has a zero-error point (the minimum is 159 errors), so dropping zero-count points does not affect them.

Audit of `src/leakdecode/*.py` for the rotated Z-memory with `HERALDED_ERASE`. Stim 1.16.0, PyMatching 2.4.0. Checks used `/workspace/builder/riverlane/.venv/bin/python` on one core. No `leakdecode-bench` run, and `results/` is absent, so this note has no study logical-error rates, timings, or qubit savings. Counts below are from these checks. Worked examples of the fitter state their inputs.

## 1. Record lookback — OK

`build_memory` rewrites every `rec` target after inserting herald measurements. For the rebuilt circuit at `p_leak = 0`:

- Detector arguments, detector coordinates, and the observable point at the same measurements as Stim's generator. Checked for `(d, rounds)` = `(3, 1)`, `(3, 3)`, `(3, 4)`, `(5, 5)`, `(7, 7)`.
- No detector reads a herald bit.
- `herald_meas_index` is exactly the measurement-record positions of `HERALDED_ERASE`, and `herald_qubit` matches those targets.
- At `p = 0`, `p_leak = 0`, 32 shots at `d = 3, 5, 7`: every detector and the observable are 0. Some raw measurements are 1. X-stabiliser outcomes of a `|0⟩` memory are random even with the noise off; the detectors are differences, and those differences are 0.
- The detector error model of the rebuilt circuit equals the model of the **flattened** generator. At `d = 3` and `d = 5` the error lists match term by term (286 and 1953 errors). At `d = 7` both flattened models have 6092 errors and `==` is true. Their PyMatching graphs have the same 1558 edges, the same probabilities (maximum absolute difference 0), and the same fault ids.

Comparing the rebuilt model with the **unflattened** generator fails at `d = 7` (6602 errors against 6092). That gap is already there between Stim's `REPEAT` circuit and Stim's own `flattened()` circuit. The memory experiment uses the rebuilt flat circuit, and that circuit matches the flattened generator.

## 2. Where the erasures sit — OK

With `p = 2e-3` the herald instructions sit in the noise layer Stim uses for this generator:

- Data: previous op `DEPOLARIZE1` on all `d²` data qubits, next op `H` on the X ancillas. One such instruction per round, immediately where `before_round_data_depolarization` is applied. The two Pauli channels are on the same qubits with only that layer between them, so they commute.
- Ancilla: previous op `X_ERROR` on the measure qubits, next op `MR`. That `X_ERROR` is `before_measure_flip_probability`. The basis-change `H` is earlier in the round, so the erasure is in the Z measurement basis, immediately before `MR`.

With `p = 0` Stim drops the noise ops. The same slots are then `TICK` then the herald then `H`, and `TICK` then the herald then `MR`.

Counts match the claim. Every `MR` in a circuit has the same targets. For `d = 7`, rounds `0..6` each have 49 data heralds and 48 ancilla heralds (679 in total). Nothing is inserted after the final data `M`. The model therefore has no heralded erasure on the final data measurement, and none on the CNOTs. That matches the comment in `circuits.py`.

At `p = 0`, `p_leak = 0.2`, `d = 3`, 400 shots: the 51 herald columns have mean 0.203, column means from 0.155 to 0.260. Those columns are the erasure record.

## 3. Herald edges — OK

`herald_signatures` injects X and Z with a noiseless `FlipSimulator` just after each `HERALDED_ERASE`. Independent checks:

- For every herald at `d = 3` (51) and `d = 5` (245), a fresh circuit with `X_ERROR(1)`, `Y_ERROR(1)`, or `Z_ERROR(1)` placed immediately after that herald's instruction was sampled with `compile_detector_sampler`. Detector and observable flips matched the library edges. Mismatches: 0. The same patterns appear if the circuit is first built at `p = 2e-3` and the Pauli noise is then set to probability 0, so the signatures belong to the circuit that is actually decoded.
- The same X/Y/Z circuits, 6 heralds at `d = 3` (indices 0, 4, 9, 20, 42, 50), run in `TableauSimulator` and converted with `compile_m2d_converter`: 0 mismatches out of 18. `TableauSimulator` does execute `HERALDED_ERASE(0)`; the comparison above replaces the random Pauli by a definite gate so the two engines see one circuit.
- Real `HERALDED_ERASE(1)` on one qubit and 0 on the others, `p = 0`, 800 shots, `d = 3`. Every shot's `(detectors, observable)` was one of the four Pauli patterns. Unknown shots: 0.

| herald | kind | round | shot classes (800) |
| --- | --- | --- | --- |
| k=0, qubit 1 | data | 0 | I 405, X 395 |
| k=42, qubit 19 | data | 2 | I 201, X 193, Z 209, Y 197 |
| k=50, qubit 25 | ancilla | 2 | I 410, X 390 |

Round-2 data is split four ways, as a uniform `{I, X, Y, Z}` channel should be. Round-0 data and every ancilla show two classes of about half the shots. That is the same channel, with a simpler action:

- Every ancilla `Z` signature is empty (`d = 3, 5, 7`, every round). `MR` is a Z measurement, so Z does nothing and Y acts as X. `P(edge flips | herald) = 1/2` still.
- Every round-0 data `Z` signature is empty (9, 25, and 49 qubits). Those qubits have just been reset to `|0⟩`, so Z does not flip a later Z measurement, and it flips every later X-stabiliser measurement equally, which cancels in the detectors. Y acts as X. Again the surviving edge flips with probability 1/2.

X and Z detector sets are disjoint for every herald at `d = 3` and `d = 5`, and for the full `d = 7` list the X and Z supports each have length at most 2 (`max` 2). A Y often lights 3 or 4 detectors because it is the XOR of those two edges. `Edge.from_dets` only has to accept X and Z, and it does: `herald_signatures` completes for `d = 7` (679 heralds) and matches a second `FlipSimulator` pass with 0 mismatches. No herald is an undetectable logical.

## 4. Adaptive reweighting — OK

Given a herald, the qubit is I, X, Y, or Z with probability 1/4 each. The X-edge (Z detectors) flips on X or Y, the Z-edge on Z or Y. Each marginal is 1/2. For this channel the joint is independent as well: `P(both) = P(Y) = 1/4 = (1/2)·(1/2)`. Setting both edges to probability 1/2 is the exact conditional model. There is no Y-correlation left to ignore for that herald.

If any other fault flips the same edge with probability `q`, independent of the herald,

`P(flip) = (1/2)(1 − q) + (1/2)q = 1/2`.

The weight `log((1 − 1/2)/(1/2))` is 0 whatever `q` was. Replacing the whole merged edge by `p = 1/2` is the right update. Two or three heralds on one edge still leave probability 1/2, so writing it once is enough. At `d = 3`, 14 edges are shared by two or more heralds (maximum 3). At `d = 5`, 40 edges, maximum 3.

The static graph is the DEM at `p_leak = 0` with the same circuit `p`. That is the right prior on the shots where a given herald did not fire. Checked consequences:

- At `p = 2e-3`, every herald edge is already in that graph (`absent = 0` at `d = 3, 5, 7`). The `p = 1e-12` fallback is not used in the study's noise.
- For every herald edge the merged `fault_ids` bit equals `e.obs`. Mismatches: 0 at `d = 3` (51 edges), `d = 5` (297), `d = 7` (879). No pair `(u, v)` carries two different observable bits. At `d = 3` the `p_leak = 0` DEM has 78 endpoint groups and 0 groups whose graphlike pieces disagree on the observable. PyMatching keeps the first fault id when it merges parallel edges; here there is only one label to keep, and it is the herald's label.
- `add_edge(..., merge_strategy="replace")` sets `weight = 0`, `error_probability = 0.5`, and the original fault ids. Restoring the saved triple puts the original weight and probability back (checked on a boundary edge at `d = 3, 5, 7`).

Single-herald decode, graph built at `p = 2e-3`, syndrome taken from a pure I, X, Y, or Z at that herald (`p = 0` circuit): **0 failures out of 204** at `d = 3` and **0 out of 980** at `d = 5`. Five two-X pairs at `d = 3`, `(0,1)`, `(0,42)`, `(9,50)`, `(20,21)`, `(4,46)`, also matched the observable. `d = 7` was checked for fault ids and graphlike signatures, not for an exhaustive decode.

Low caveat. Minimum-weight matching returns one lowest-weight correction. Several weight-0 corrections can differ by a logical operator when many heralds fire. The decoder then picks one class rather than summing them. A single herald never did this in the tests above (every all-identity syndrome decoded to 0). This is the usual MWPM gap on erasures, including the Delfosse–Zémor picture, and both Python decoders share it.

## 5. Wilson intervals and the per-round rate — OK, with two small edges

`wilson_interval` matches the usual score interval. Checked against an independent implementation for `(k, n) = (0, 100), (1, 10), (10, 100), (500, 100000)`.

`per_round` inverts `P = (1 − (1 − 2e)^r) / 2`, which is the probability of an odd number of independent round-flips. Round trips at `(e, r) = (1e-4, 3), (0.02, 5), (0.1, 7)` come back with absolute error about `1e-17`. Using it on a memory with `rounds = d` is the standard effective rate. It folds reset, bulk, and the final measurement into one number per round. `sim.run_sweep` passes `rounds = d`, which matches the circuits it builds.

Applying `per_round` to the Wilson endpoints is valid while the map is monotone, which it is on `[0, 1/2]`.

**ISSUE low.** `per_round` clamps its input to `[0, 1/2]` (`stats.py`). `per_round(0.9, 3)` returns 0.5. A Wilson upper bound above 1/2 is reported as a per-round upper bound of 1/2, so the interval is censored. At the study's noise this is unlikely. The clamp should be left visible, or the per-round interval should be left undefined when the shot-level bound exceeds 1/2.

**ISSUE low.** `sim.run_sweep` stops when both decoders have reached `max_errors`, or the shot cap is hit (`sim.py`, the `min(err_static, err_adaptive) >= max_errors` test). `n` is random, and a wave of chunks can run past the stopping point. The Wilson interval is then computed as if `n` had been fixed. With `max_errors = 500` the bias of a negative-binomial style estimate is of order `1/500`. State that the interval assumes a fixed sample size, or stop on a fixed shot budget.

For `k = 0` the floating-point lower endpoint is about `1e-18` rather than a clean 0 (`0/100` gave `3.47e-18`). The point estimate is still 0, and the fit drops that point because the point estimate is 0. No change needed.

## 6. The fit and the qubit saving — ISSUE med

`fit_exponential_suppression` fits `log10(eps) = a + b d` by weighted least squares. `Lambda = 10^{-2b}` equals `eps(d) / eps(d+2)`, which is also the `Λ` in `eps ∝ Λ^{-d/2}`. A straight line `a = -2`, `b = -0.3` is recovered (`Lambda` matches `10^{0.6}` to the last digits the fit prints). The algebra of a line that is really there is fine.

The distance solver for a per-round target is the straight inverse. For `a = 0`, `b = -1`, target `1e-5` it returns `5.0`, and the odd-distance rounding stays on 5.

**The fit throws away every distance with no observed errors** (`stats.py`: keep points with `e > 0` and `lo > 0`; for integer counts those are the points with `k ≥ 1`). The upper bound from a clean large distance never enters the fit. Worked example, not a memory run: distances 3, 5, 7; errors 200, 150, 0; shots 20 000, 20 000, 200 000. The mask keeps 2 of the 3 points. The line predicts a per-round rate at `d = 7` of `6.79e-4`. The unused per-round Wilson upper bound on that point is `2.74e-6`. The prediction sits above the data that were dropped.

Weights are `1/sigma²` with `sigma = (log10(hi) − log10(lo)) / (2 Z_95)`. That treats the transformed Wilson interval as a symmetric Gaussian. For 5 errors in 1000 shots the two log10 gaps are 0.369 and 0.367, so the shape is mild once the count is not tiny. The zero-count rule is the part that changes the fit.

**The plotted saving is not a device qubit count.** `analysis.savings_from_fits` sets `qubits = 2 d² − 1` with the real-valued `d` from the line. `plots.plot_qubit_savings` draws `qubit_ratio`, the ratio of those real-valued counts. `qubit_ratio_odd` is computed and then left out of the figure.

Worked example, per-round target `1e-6`, both slopes `b = -0.5`, intercepts `-3.45` and `-3.5`:

| | static | adaptive |
| --- | --- | --- |
| continuous `d` | 5.1 | 5.0 |
| next odd distance | 7 | 5 |
| continuous qubits | 51.02 | 49 |
| odd-patch qubits (`2d²−1`) | 97 | 49 |

`qubit_ratio` is 1.04. `qubit_ratio_odd` is 1.98. One hundredth of a distance past an odd integer moves a real patch from `d` to `d+2`, and the quadratic count jumps. `num_qubits_rotated` itself is right: the circuits have 17, 49, and 97 qubits at `d = 3, 5, 7`.

The window target solves `d · eps(d) = target`, i.e. `log10(d) + a + b d = log10(target)`. For a target of `1e-6` the step from `d·e` to the exact shot error `(1−(1−2e)^d)/2` is negligible. The solver is a bisection from 1 to 1000. The comment says the left-hand side is monotone for `b < 0`. The derivative `1/(d ln 10) + b` is positive at small `d` when `b > -1/ln(10) ≈ -0.434`, i.e. when `10^{-2b}` is below about 7.4. Then `d·eps(d)` rises and falls, and two distances can hit the target.

Examples that do have two roots, and where the larger root is the distance above which the window stays under the target:

- `a = -6.2`, `b = -0.05` (`10^{-2b} ≈ 1.26`): roots near 2 and 23.5. Solver returns 23.37.
- `a = -5.9`, `b = -0.15` (`10^{-2b} ≈ 2.0`): roots near 1.25 and 5.75. Solver returns 5.71.

A short sweep of shallow slopes, with the peak set just above the target, did not produce a wrong root. The routine still implements an assumption the derivative does not obey. Suggested fix: define the answer as the larger root, and find it by a scan plus a local bisection. If `b ≥ -1/ln(10)`, say so on the row.

Further limits, all of which sit on top of those code choices:

- The configured sweep has three odd distances, 3, 5, and 7. A line through two or three points has no lack-of-fit check. The per-round rate at `d = 3` still includes the time boundaries.
- `a` and `b` are returned without an error bar, and the qubit ratio inherits none.
- If the static slope is non-negative the distance becomes infinity and the ratio becomes NaN, so the plot drops the point (`plots.py` skips `nan` and `inf`). A run in which only the adaptive decoder is below threshold disappears from the figure.
- `d_odd` is at least 3. The continuous qubit count is not clamped. Both decoders can then sit on distance 3 in the odd columns while the continuous ratio still claims a saving.

Suggested fix for the fit: keep zero-error distances as upper bounds (a binomial or Poisson likelihood on the raw `k` and `n`, or a censored least-squares term). Report `qubit_ratio_odd` in the figure, with a range from the fit covariance or from a bootstrap of the shots. Keep the continuous `d` as a sketch of the line only.

## 7. Is the static decoder a fair baseline? — OK

`BaselineMatcher` builds PyMatching on the DEM of the circuit that still contains `HERALDED_ERASE(p_leak)`, and it never reads the herald array (`decode_batch` ignores it). No detector contains a herald bit (section 1), so the syndrome does not smuggle the locations in. Both decoders see the same shots. The static decoder is told the average erasure rate. The adaptive decoder is told the locations and uses the `p_leak = 0` graph plus the `p = 1/2` update. That is the right pair for the question "how much do the herald bits buy, once the rate is already known?"

It is a strong control. A graph built at `p_leak = 0` and then run on leaky shots, with the heralds ignored, leaves the erasures out of the prior. That control fails more often and makes the adaptive gain look larger. The study uses the stronger one.

Two small qualifications.

**ISSUE low.** The comment on `BaselineMatcher` says Stim folds each erasure into independent X, Y, and Z faults at `p_leak/4`. Stim will not build a DEM for `HERALDED_ERASE` unless `approximate_disjoint_errors=True` (it raises otherwise). What it emits is not three independent `p_leak/4` terms. A one-qubit check, `HERALDED_ERASE(0.2)` then a Z measurement, gives a single `error(0.1) D0 L0`, i.e. the marginal `p_leak/2` of X or Y. On the `d = 3` memory with only erasures, edge probabilities in the DEM differ from the exact XOR of per-herald `p_leak/2` contributions by at most `2.5e-5` at `p_leak = 0.01`, `9.8e-5` at `p_leak = 0.02`, and `8.1e-3` at `p_leak = 0.2` (edge probabilities about 0.010, 0.020, and 0.18 on the worst edges). With circuit noise `p = 2e-3` as well, the same comparison at `p_leak = 0.02` has maximum absolute error `9.6e-5` and median about `1e-18`. Fault ids do not change. The full bench's largest configured `p_leak` is `2e-2` (`bench.py`), so the weight error is small next to the Monte Carlo noise. Suggested fix: describe the baseline as "Stim's disjoint DEM of the leaky circuit", and if exact marginals are wanted, set each edge to the XOR of its `p_leak = 0` probability and the exact erasure marginal instead of calling `detector_error_model` on the leaky circuit.

The herald-blind channel also correlates the X-edge and the Z-edge of one erasure (`P(both) = p_leak/4`, while the product of the marginals is `(p_leak/2)²`). Conditional on the herald those edges are independent, which is why section 4 is exact. MWPM never uses that cross-graph correlation, on the baseline or on the adaptive decoder. Part of the adaptive gain at large `p_leak` is that the conditional model matches independent matching, while the blind model does not. For a Z-memory both edges can carry the Z observable, so a correlated decoder could in principle use the blind correlation. That is a property of choosing MWPM on both sides, and it is the same decoder family. It does not make the baseline unfair. It does mean the gap is "herald bits, inside MWPM", not "herald bits against the best herald-blind decoder".

## Summary

| # | Item | Verdict |
| --- | --- | --- |
| 1 | `rec` remapping, DEM at `p_leak = 0`, noiseless detectors | OK |
| 2 | Data erasure at round start, ancilla erasure before `MR` | OK |
| 3 | Herald edges against `X_ERROR(1)`, Tableau, and real `HERALDED_ERASE(1)` | OK |
| 4 | `p = 1/2` reweight, `p_leak = 0` graph, merged fault ids | OK (degeneracy caveat) |
| 5 | Wilson interval and per-round inversion | OK; low issues on the clamp and on optional stopping |
| 6 | Log-linear fit and qubit extrapolation | ISSUE med |
| 7 | Static decoder that knows `p_leak` but not the locations | OK; low issue on the DEM approximation and the comment |
