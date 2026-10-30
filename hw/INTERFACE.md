# Interface between the Python package and the C++ union-find decoder

The Python side (`leakdecode.export`) writes, for one (distance, p, p_leak) point:

## graph.txt  (decoding graph, single logical observable)
```
<num_detectors> <num_edges>
<u> <v> <weight_int> <obs>     # one line per edge, num_edges lines
```
* `u` is a detector index in [0, num_detectors).
* `v` is a detector index, or `-1` for the virtual boundary.
* `weight_int` is a non-negative integer weight (round(log((1-p)/p) scaled), may be ignored by an unweighted UF).
* `obs` is 0 or 1: whether flipping this edge flips logical observable 0.
* Edge index = line order (0-based).

## syndromes.u8
Raw bytes, row-major, shape (num_shots, num_detectors), each byte 0 or 1.
num_shots = filesize / num_detectors.

## erasures.u8 (optional)
Raw bytes, shape (num_shots, num_edges), each byte 0 or 1; 1 = this edge is heralded (erased) in this shot.

## Output: predictions.u8
Raw bytes, length num_shots, each 0/1 = predicted flip of observable 0.

## stdout
One JSON line: {"shots":N,"num_detectors":D,"num_edges":E,"total_decode_ns":T,"ns_per_shot":X,"max_ns_per_shot":M,"static_bytes":B}
