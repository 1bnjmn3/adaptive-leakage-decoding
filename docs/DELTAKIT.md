# Deltakit without cloud access

Checked on 30 Oct 2026. PyPI lists `deltakit` 0.10.2 (`pip install deltakit`, Python `>=3.11,<3.16`). The copy in `/tmp/dkvenv` is that release: `deltakit-decode` 0.9.2, `deltakit-explorer` 0.9.4, `deltakit-circuit` 0.9.2, `deltakit-stim` 0.2.7, PyMatching 2.4.0. Docs used: the v0.10.2 [Setup](https://deltakit-docs.riverlane.com/en/stable/setup.html) and [Getting Started](https://deltakit-docs.riverlane.com/en/stable/guide/getting_started.html) pages, plus the [leakage-aware LCD notebook](https://deltakit-docs.riverlane.com/en/stable/examples/notebooks/simulation/leakage-aware_decoding.html).

## Local decoder

Getting Started says Deltakit exposes one decoder for local use: PyMatching, class `deltakit.decode.PyMatchingDecoder`. Build it with `PyMatchingDecoder.construct_decoder_and_stim_circuit` from a `deltakit.circuit.Circuit`. Both decode methods take a syndrome and nothing else:

```text
decode_batch_to_logical_flip(self, syndrome_batch: np.ndarray) -> np.ndarray
```

`deltakit.decode` also has shot-loop helpers and noise sources. One of those, `UniformErasureNoise`, draws erased edges for a Monte Carlo model. No decoder in `deltakit_decode` reads that list.

## Cloud decoders

`LCDecoder`, `MWPMDecoder`, `CCDecoder`, `BeliefMatchingDecoder`, `BPOSDecoder` and `ACDecoder` are re-exported from `deltakit.decode` and defined in `deltakit_explorer._cloud_decoders`. Each one requires a `deltakit.explorer.Client`. With `client` left out, construction raises:

```text
NotImplementedError: Currently, a `client` must be provided to instantiate this class.
```

That message was the same for all six classes, tried on the distance-3 circuit below. Their batch method is `decode_batch_to_logical_flip(syndrome_batch, leakage_batch=None)`, and the base class forwards `leakage_batch` to `Client.decode` as `leakage_flags`. The client docstring says leakage flags may be used by the Local Clustering Decoder (LCD, arXiv:2411.10343). The docs notebook for leakage-aware decoding is `Client.get_instance()` and `DecoderType.LCD`.

A token is a personal access token from the Deltakit site. Setup says: open <https://deltakit.riverlane.com/dashboard/token>, generate a token, then once:

```python
from deltakit.explorer import Client
Client.set_token("<your token>")
```

`Client.set_token` stores `DELTAKIT_TOKEN` in the process environment and in `~/.config/deltakit-explorer/.env` on Linux. Setting that variable yourself is the other way in. `Client.get_instance()` then uses `https://deltakit.riverlane.com/proxy` unless `DELTAKIT_SERVER` is set. With no session, the dashboard URL is a sign-in / create-account page. No account was created and no token was requested for this note. If neither the file nor the variable is set, `Client.get_instance()` raises locally, before any decode call. The installed text has a typo in the host (`rivelane`, missing an e); the docs and the README use `riverlane.com`:

```text
RuntimeError: Token could not be found neither in (/home/box/.config/deltakit-explorer/.env) nor environment variable (DELTAKIT_TOKEN). Please obtain your token at https://deltakit.rivelane.com/dashboard/token and use `Client.set_token` function to register it.
```

## Leakage and heralds

`PyMatchingDecoder` has no herald or leakage argument. On the distance-3 decoder below, an extra array and a `leakage_batch=` keyword both raised `TypeError`.

This study inserts Stim `HERALDED_ERASE`. `Circuit.from_stim_circuit` does not implement that gate, so the local constructor cannot be given these circuits:

```text
InstructionNotImplemented: Parsing of 'HERALDED_ERASE(0.001) 0' is not implemented yet.
```

Deltakit's own `LEAKAGE` and `HERALD_LEAKAGE_EVENT` do parse, and the local `deltakit_stim` sampler writes the herald into the measurement record (`LEAKAGE(1)` produced a herald bit of 1 on four shots; `LEAKAGE(0)` produced 0). `simulate_with_stim` without a client still returns `leakage=None`. Building `PyMatchingDecoder` on a small reset–`LEAKAGE`–measure circuit succeeded only in the sense that the call returned: the detector error model had a detector and an observable and no `error` lines, and construction warned that the logical had no activators. The herald column is not passed into the matcher. The docs treat leakage simulation as a server call.

## Commands and output

Interpreter: `/tmp/dkvenv/bin/python`. Circuit: `deltakit_stim.Circuit.generated("surface_code:rotated_memory_z", distance=3, rounds=3)` with `after_clifford_depolarization`, `before_round_data_depolarization`, `before_measure_flip_probability` and `after_reset_flip_probability` all `2e-3`. Then `Circuit.from_stim_circuit`, `PyMatchingDecoder.construct_decoder_and_stim_circuit`, and `compile_detector_sampler(seed=1).sample(shots=2, separate_observables=True)`. The helper trims the circuit (raw Stim circuit: 24 detectors; decoded shot: 16).

```text
constructed PyMatchingDecoder PyMatching2
decode_batch_to_logical_flip (self, syndrome_batch: np.ndarray) -> np.ndarray
dets.shape (2, 16) obs.shape (2, 1)
preds.shape (2, 1) preds.dtype uint8
preds [[0], [0]]
TypeError: PyMatchingDecoder.decode_batch_to_logical_flip() takes 2 positional arguments but 3 were given
TypeError: PyMatchingDecoder.decode_batch_to_logical_flip() got an unexpected keyword argument 'leakage_batch'
```

`LCDecoder(sc)` with no client produced the `NotImplementedError` quoted above. The same call on the other five cloud classes produced the same sentence.

## What this project uses

This study does not import Deltakit. The open install can run ordinary MWPM on a Stim memory circuit that contains only gates the circuit parser knows. It cannot take the per-shot `HERALDED_ERASE` heralds this study decodes, and `LCDecoder`, the leakage-aware decoder in the docs, will not construct without a `Client` and an account token. The baseline matcher, the per-shot reweighting to probability 1/2, and the C++ union-find decoder are therefore our own code (`src/leakdecode/decoders.py` and `hw/`).
