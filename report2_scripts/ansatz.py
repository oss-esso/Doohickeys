from __future__ import annotations

from typing import Callable

import numpy as np
import pennylane as qml


def build_uccsd_ansatz(n_qubits: int, n_electrons: int,
                       hf_state: np.ndarray | None = None
                       ) -> tuple[Callable[[np.ndarray], None], int, list]:
    """Build a UCCSD (Unitary Coupled-Cluster Singles & Doubles) ansatz circuit.

    If hf_state is None, constructs it as a binary vector of length n_qubits
    with the first n_electrons entries set to 1 (occupied spin-orbitals).

    Uses qml.qchem.excitations(n_electrons, n_qubits) to generate:
      - singles: list of [occupied, virtual] index pairs
      - doubles: list of [occ1, occ2, virt1, virt2] index quadruples
    Concatenate: excitations = singles + doubles, n_params = len(excitations).

    The returned ansatz_fn(params) applies:
      1. qml.BasisState(hf_state, wires=range(n_qubits))
      2. qml.SingleExcitation(params[i], wires=exc) for each single
      3. qml.DoubleExcitation(params[i], wires=exc) for each double

    Note: For H2 (4 qubits, 2 electrons):
      singles = [[0,2], [1,3]], doubles = [[0,1,2,3]], n_params = 3.
      By symmetry, single excitation amplitudes are zero at the optimum.

    Args:
        n_qubits: Number of qubits (= number of spin-orbitals).
        n_electrons: Number of electrons.
        hf_state: Optional explicit HF reference state vector.

    Returns:
        Tuple of (ansatz_fn, n_params, excitations) where excitations
        is the concatenated list of singles + doubles.
    """
    if hf_state is None:
        # Create HF reference: first n_electrons spin-orbitals occupied
        hf_state = np.zeros(n_qubits, dtype=int)
        hf_state[:n_electrons] = 1

    singles, doubles = qml.qchem.excitations(n_electrons, n_qubits)
    excitations = singles + doubles
    n_params = len(excitations)

    def ansatz_fn(params):
        qml.BasisState(hf_state, wires=range(n_qubits))
        for i, excitation in enumerate(excitations):
            if len(excitation) == 2:
                qml.SingleExcitation(params[i], wires=excitation)
            else:
                qml.DoubleExcitation(params[i], wires=excitation)

    return ansatz_fn, n_params, excitations


def build_hardware_efficient_ansatz(n_qubits: int, n_layers: int = 2,
                                     entangler: str = "CNOT"
                                     ) -> tuple[Callable[[np.ndarray], None], int]:
    """Build a hardware-efficient ansatz with alternating rotation and entangling layers.

    Parameter count: 2 * n_qubits * n_layers (one RY + one RZ per qubit per layer).
    Reshape params to (n_layers, n_qubits, 2).

    Each layer applies:
      1. Single-qubit rotations: RY(p[l,q,0]) and RZ(p[l,q,1]) on each qubit q.
      2. Entangling gates: nearest-neighbour CNOT or CZ on wires [q, q+1]
         for q in range(n_qubits - 1).

    Args:
        n_qubits: Number of qubits.
        n_layers: Number of variational layers (default 2).
        entangler: Type of entangling gate: "CNOT" or "CZ".

    Returns:
        Tuple of (ansatz_fn, n_params).
    """
    n_params = 2 * n_qubits * n_layers
    
    def ansatz_fn(params):
        p = params.reshape(n_layers, n_qubits, 2)
        for l in range(n_layers):
            for q in range(n_qubits):
                qml.RY(p[l, q, 0], wires=q)
                qml.RZ(p[l, q, 1], wires=q)
            for q in range(n_qubits - 1):
                if entangler == "CNOT":
                    qml.CNOT(wires=[q, q + 1])
                elif entangler == "CZ":
                    qml.CZ(wires=[q, q + 1])
    
    return ansatz_fn, n_params
    
