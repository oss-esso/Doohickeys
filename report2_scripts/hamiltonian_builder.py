from __future__ import annotations

from dataclasses import dataclass
import itertools

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
    n_mo = one_electron_integrals.shape[0]
    n_spin_orbs = 2 * n_mo
    
    # Build fermionic Hamiltonian using FermiSentence for proper addition
    fermi_terms = qml.fermi.FermiSentence({})
    
    # One-electron terms: h_ij * a†_pσ a_qσ
    for i in range(n_mo):
        for j in range(n_mo):
            h_ij = one_electron_integrals[i, j]
            if np.abs(h_ij) > 1e-12:
                for sigma in range(2):
                    p = 2 * i + sigma
                    q = 2 * j + sigma
                    fw = qml.fermi.FermiWord({(0, p): "+", (1, q): "-"})
                    fermi_terms += qml.fermi.FermiSentence({fw: h_ij})

    # Two-electron terms: 0.5 * g_ijkl * a†_pσ a†_qτ a_sτ a_rσ
    # Using chemists' notation (ij|kl) with physics ordering
    for i, j, k, l in itertools.product(range(n_mo), repeat=4):
        g_ijkl = two_electron_integrals[i, j, k, l]
        if np.abs(g_ijkl) > 1e-12:
            for s1 in range(2):
                for s2 in range(2):
                    p = 2 * i + s1
                    q = 2 * k + s2
                    r = 2 * j + s1
                    s = 2 * l + s2
                    if p == q or r == s:
                        continue
                    fw = qml.fermi.FermiWord({(0, p): "+", (1, q): "+", (2, s): "-", (3, r): "-"})
                    fw = qml.fermi.FermiWord({(0, p): "+", (1, q): "+", (2, s): "-", (3, r): "-"})
                    fermi_terms += qml.fermi.FermiSentence({fw: 0.5 * g_ijkl})

    if mapping == "jordan_wigner":
        qubit_hamiltonian = qml.jordan_wigner(fermi_terms)
    elif mapping == "parity":
        qubit_hamiltonian = qml.parity_transform(fermi_terms, n_spin_orbs)
    else:
        raise ValueError(f"Unsupported mapping: {mapping}")
    
    # Add nuclear repulsion as constant term
    qubit_hamiltonian = qubit_hamiltonian + nuclear_repulsion * qml.Identity(0)
    qubit_hamiltonian = qml.simplify(qubit_hamiltonian)

    # Count Pauli terms
    if hasattr(qubit_hamiltonian, 'coeffs'):
        n_paulis = len(qubit_hamiltonian.coeffs)
    else:
        n_paulis = 1

    return QubitHamiltonian(
        hamiltonian=qubit_hamiltonian,
        n_qubits=n_spin_orbs,
        n_paulis=n_paulis,
        mapping=mapping,
        fci_energy=None)
    


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
    H_qt = pennylane_to_qutip(qubit_hamiltonian.hamiltonian, qubit_hamiltonian.n_qubits)
    eigenvalues = H_qt.eigenenergies()
    print(f'Energies from QuTiP diagonalisation: {eigenvalues}')
    print(f'Ground-state energy from QuTiP: {eigenvalues[0]}')
    print(f'FCI energy from PySCF: {qubit_hamiltonian.fci_energy}')
    qubit_hamiltonian.fci_energy = eigenvalues[0]
    return qubit_hamiltonian.fci_energy


def pennylane_to_qutip(hamiltonian, n_qubits: int) -> qutip.Qobj:
    """Convert a PennyLane Hamiltonian to a QuTiP Qobj.

    Iterates over (coeff, obs) pairs from hamiltonian.coeffs and hamiltonian.ops.
    For each term, builds a Pauli string of length n_qubits (default all "I"):
      - qml.Identity -> all "I"
      - Single Pauli (PauliX/Y/Z) -> set wire index to "X"/"Y"/"Z"
      - Tensor product (obs.obs or Prod) -> set each sub-observable's wire

    Maps labels to QuTiP: "I"->qeye(2), "X"->sigmax(), "Y"->sigmay(), "Z"->sigmaz().
    Each term: float(coeff) * qutip.tensor([op for each qubit]).

    Args:
        hamiltonian: PennyLane Hamiltonian with .coeffs and .ops (Sum, LinearCombination, etc.).
        n_qubits: Total number of qubits.

    Returns:
        QuTiP Qobj of dimension (2^n_qubits, 2^n_qubits).
    """
    qt_hamiltonian = qutip.tensor([qutip.qeye(2) for _ in range(n_qubits)]) * 0  # Zero matrix
    
    # Handle different PennyLane Hamiltonian types
    if hasattr(hamiltonian, 'coeffs') and hasattr(hamiltonian, 'ops'):
        coeffs = hamiltonian.coeffs
        ops = hamiltonian.ops
    elif hasattr(hamiltonian, 'terms'):
        # For LinearCombination or similar
        coeffs, ops = hamiltonian.terms()
    else:
        raise ValueError(f"Unsupported Hamiltonian type: {type(hamiltonian)}")
    
    for coeff, obs in zip(coeffs, ops):
        pauli_string = ["I"] * n_qubits
        
        def _extract_pauli(op, pauli_str):
            """Recursively extract Pauli operators from an observable."""
            if isinstance(op, qml.Identity):
                pass  # Identity leaves the pauli_str as "I"
            elif isinstance(op, (qml.X, qml.PauliX)):
                pauli_str[op.wires[0]] = "X"
            elif isinstance(op, (qml.Y, qml.PauliY)):
                pauli_str[op.wires[0]] = "Y"
            elif isinstance(op, (qml.Z, qml.PauliZ)):
                pauli_str[op.wires[0]] = "Z"
            elif hasattr(op, 'obs'):  # Tensor product (old style)
                for sub_obs in op.obs:
                    _extract_pauli(sub_obs, pauli_str)
            elif hasattr(op, 'operands'):  # Prod (new style)
                for sub_obs in op.operands:
                    _extract_pauli(sub_obs, pauli_str)
            elif hasattr(op, 'base'):  # SProd or similar
                _extract_pauli(op.base, pauli_str)
            else:
                raise ValueError(f"Unsupported observable type: {type(op)}")
        
        _extract_pauli(obs, pauli_string)
        
        qt_ops = []
        for label in pauli_string:
            if label == "I":
                qt_ops.append(qutip.qeye(2))
            elif label == "X":
                qt_ops.append(qutip.sigmax())
            elif label == "Y":
                qt_ops.append(qutip.sigmay())
            elif label == "Z":
                qt_ops.append(qutip.sigmaz())
        qt_hamiltonian += float(coeff) * qutip.tensor(qt_ops)
    return qt_hamiltonian
