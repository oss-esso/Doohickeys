from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pennylane as qml
import qutip


@dataclass
class QubitHamiltonian:
    """Container for a qubit Hamiltonian and metadata."""
    hamiltonian: object  # qml.Hamiltonian
    n_qubits: int
    n_paulis: int
    mapping: str
    fci_energy: float | None = None


def build_hamiltonian(one_electron_integrals: np.ndarray,
                      two_electron_integrals: np.ndarray,
                      n_electrons: int,
                      nuclear_repulsion: float,
                      mapping: str = "jordan_wigner",
                      active_electrons: int | None = None,
                      active_orbitals: int | None = None) -> QubitHamiltonian:
    """Build a qubit Hamiltonian from molecular integrals.

    Constructs the second-quantised fermionic Hamiltonian from one- and
    two-electron integrals in the MO basis, then maps to qubit operators.

    Steps:
      1. Compute n_spin_orbs = 2 * n_mo (each spatial orbital -> 2 spin-orbitals).
      2. Build one-electron terms using qml.fermi.FermiWord:
         For each spatial pair (i,j) with |h_ij| > 1e-12, for sigma in {0,1}:
           p = 2*i + sigma, q = 2*j + sigma  (same spin, spin-diagonal)
           Add h_ij * FermiWord({(0,p): "+", (1,q): "-"})
      3. Build two-electron terms:
         For each (i,j,k,l) with |h_ijkl| > 1e-12, for s1,s2 in {0,1}:
           p = 2*i+s1, q = 2*k+s2, r = 2*j+s1, s = 2*l+s2
           Skip if p == q.
           Add 0.5 * h_ijkl * FermiWord({(0,p):"+", (1,q):"+", (2,s):"-", (3,r):"-"})
      4. Map to qubits:
         "jordan_wigner" -> qml.jordan_wigner(fermi_H)
         "parity" -> qml.parity_transform(fermi_H, n_spin_orbs)
      5. Add nuclear_repulsion * qml.Identity(0) and qml.simplify.

    Args:
        one_electron_integrals: h_pq in MO basis, shape (n_mo, n_mo).
        two_electron_integrals: h_pqrs in chemists' notation, shape (n_mo,)*4.
        n_electrons: Number of electrons.
        nuclear_repulsion: Nuclear repulsion energy in Hartrees.
        mapping: Fermion-to-qubit mapping ("jordan_wigner" or "parity").
        active_electrons: Placeholder for active-space electron count.
        active_orbitals: Placeholder for active-space orbital count.

    Returns:
        QubitHamiltonian with PennyLane qubit operator and metadata.
    """
    ...


def validate_with_qutip(qubit_hamiltonian: QubitHamiltonian) -> float:
    """Validate the qubit Hamiltonian by exact diagonalisation in QuTiP.

    Converts the PennyLane Hamiltonian to a QuTiP Qobj via pennylane_to_qutip,
    computes all eigenvalues with H_qt.eigenenergies(), and stores the
    ground-state energy as qubit_hamiltonian.fci_energy (mutates input).

    Args:
        qubit_hamiltonian: QubitHamiltonian to validate (fci_energy is updated).

    Returns:
        Ground-state energy (smallest eigenvalue) as float.
    """
    ...


def pennylane_to_qutip(hamiltonian: qml.ops.op_math.Sum,
                       n_qubits: int) -> qutip.Qobj:
    """Convert a PennyLane Hamiltonian to a QuTiP Qobj.

    Iterates over (coeff, obs) pairs from hamiltonian.coeffs and hamiltonian.ops.
    For each term, builds a Pauli string of length n_qubits (default all "I"):
      - qml.Identity -> all "I"
      - Single Pauli (PauliX/Y/Z) -> set wire index to "X"/"Y"/"Z"
      - Tensor product (obs.obs) -> set each sub-observable's wire

    Maps labels to QuTiP: "I"->qeye(2), "X"->sigmax(), "Y"->sigmay(), "Z"->sigmaz().
    Each term: float(coeff) * qutip.tensor([op for each qubit]).

    Args:
        hamiltonian: PennyLane Hamiltonian with .coeffs and .ops.
        n_qubits: Total number of qubits.

    Returns:
        QuTiP Qobj of dimension (2^n_qubits, 2^n_qubits).
    """
    ...
