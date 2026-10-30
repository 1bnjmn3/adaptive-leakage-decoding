"""Build rotated surface-code memory circuits with heralded-erasure "leakage".

We start from Stim's generated rotated surface-code memory circuit (circuit-level
depolarising noise of strength ``p`` on every gate, reset, measurement and idle
data qubit) and insert ``HERALDED_ERASE(p_leak)`` at two kinds of location:

* ``data``    : on every data qubit at the start of every round (the same place
                Stim applies ``before_round_data_depolarization``).
* ``ancilla`` : on every measure (ancilla) qubit just before it is measured.

``HERALDED_ERASE(q)`` means: with probability ``q`` the qubit suffers a uniformly
random Pauli (I, X, Y or Z, each 1/4) AND a herald bit is written to the
measurement record saying "this qubit was erased here". Without the herald, an
erasure is just extra depolarising noise. With it, the decoder knows *where*
the error might be. This is our (simplified) proxy for leakage that is detected
(e.g. by a leakage-detection readout) and immediately reset/returned.

Inserting new measurement-producing instructions shifts every ``rec[-k]``
lookback in DETECTOR / OBSERVABLE_INCLUDE lines, so we rebuild the circuit on a
flattened copy and remap every record target to its new position.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import stim

LEAK_LOCATIONS = ("data", "ancilla")


@dataclass
class LeakyMemory:
    """A memory experiment plus bookkeeping about where the heralds are."""

    distance: int
    rounds: int
    p: float
    p_leak: float
    leak_on: tuple[str, ...]
    circuit: stim.Circuit
    # one entry per herald bit, in measurement-record order
    herald_meas_index: np.ndarray  # absolute index in the measurement record
    herald_qubit: np.ndarray
    herald_round: np.ndarray
    herald_kind: np.ndarray  # 0 = data, 1 = ancilla
    data_qubits: list[int] = field(default_factory=list)
    ancilla_qubits: list[int] = field(default_factory=list)

    @property
    def num_heralds(self) -> int:
        return len(self.herald_meas_index)

    @property
    def num_physical_qubits(self) -> int:
        return len(self.data_qubits) + len(self.ancilla_qubits)

    def with_leak_rate(self, p_leak: float) -> "LeakyMemory":
        """Same circuit layout and heralds, different erasure probability."""
        return build_memory(self.distance, self.rounds, self.p, p_leak, self.leak_on)


def _base_circuit(distance: int, rounds: int, p: float) -> stim.Circuit:
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


def build_memory(
    distance: int,
    rounds: int,
    p: float,
    p_leak: float,
    leak_on: tuple[str, ...] | list[str] = LEAK_LOCATIONS,
) -> LeakyMemory:
    """Return a Z-basis rotated-surface-code memory with heralded erasures.

    The HERALDED_ERASE instructions are always present (even when
    ``p_leak == 0``) so that the circuit layout, measurement record and
    herald indices are identical for every leakage rate.
    """
    if distance < 3 or distance % 2 == 0:
        raise ValueError("distance must be an odd integer >= 3")
    if not 0 <= p_leak <= 1 or not 0 <= p <= 1:
        raise ValueError("probabilities must be in [0, 1]")
    leak_on = tuple(leak_on)
    for loc in leak_on:
        if loc not in LEAK_LOCATIONS:
            raise ValueError(f"unknown leak location {loc!r}")

    base = _base_circuit(distance, rounds, p).flattened()

    # Identify qubit roles from the circuit itself.
    mr_targets: list[int] = []
    data_targets: list[int] = []
    for inst in base:
        if inst.name == "MR":
            mr_targets = [t.value for t in inst.targets_copy()]
        elif inst.name == "M":
            data_targets = [t.value for t in inst.targets_copy()]
    assert mr_targets and data_targets

    out = stim.Circuit()
    old_to_new: list[int] = []  # old measurement index -> new measurement index
    n_old = 0
    n_new = 0
    herald_idx, herald_q, herald_r, herald_k = [], [], [], []
    round_index = 0
    awaiting_round_start = True

    def add_erasure(qubits: list[int], kind: int) -> None:
        nonlocal n_new
        out.append("HERALDED_ERASE", qubits, p_leak)
        for q in qubits:
            herald_idx.append(n_new)
            herald_q.append(q)
            herald_r.append(round_index)
            herald_k.append(kind)
            n_new += 1

    for inst in base:
        name = inst.name
        # Round start: the first H on the ancillas after the previous MR / resets.
        if name == "H" and awaiting_round_start:
            if "data" in leak_on:
                add_erasure(data_targets, 0)
            awaiting_round_start = False
        if name == "MR" and "ancilla" in leak_on:
            add_erasure(mr_targets, 1)

        if name in ("DETECTOR", "OBSERVABLE_INCLUDE"):
            new_targets = []
            for t in inst.targets_copy():
                if t.is_measurement_record_target:
                    old_abs = n_old + t.value  # t.value is negative lookback
                    new_abs = old_to_new[old_abs]
                    new_targets.append(stim.target_rec(new_abs - n_new))
                else:
                    new_targets.append(t)
            out.append(name, new_targets, inst.gate_args_copy())
            continue

        out.append(inst)
        nm = inst.num_measurements
        for _ in range(nm):
            old_to_new.append(n_new)
            n_old += 1
            n_new += 1
        if name == "MR":
            round_index += 1
            awaiting_round_start = True

    assert out.num_measurements == n_new
    assert out.num_detectors == base.num_detectors

    return LeakyMemory(
        distance=distance,
        rounds=rounds,
        p=p,
        p_leak=p_leak,
        leak_on=leak_on,
        circuit=out,
        herald_meas_index=np.asarray(herald_idx, dtype=np.int64),
        herald_qubit=np.asarray(herald_q, dtype=np.int64),
        herald_round=np.asarray(herald_r, dtype=np.int64),
        herald_kind=np.asarray(herald_k, dtype=np.int8),
        data_qubits=sorted(data_targets),
        ancilla_qubits=sorted(mr_targets),
    )


def num_qubits_rotated(distance: int) -> int:
    """Physical qubits in one rotated surface-code patch: d^2 data + d^2-1 ancillas."""
    return 2 * distance * distance - 1
