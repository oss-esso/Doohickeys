from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class PESData:
    """Container for PES scan results."""
    bond_lengths: np.ndarray
    energies_vqe: np.ndarray
    energies_hf: np.ndarray
    energies_fci: np.ndarray
    energies_cisd: np.ndarray | None = None
    energies_ccsd_t: np.ndarray | None = None
    vqe_params: list[np.ndarray] | None = None


def scan_pes(molecule_factory: Callable[[float], object],
             bond_lengths: np.ndarray,
             basis: str = "sto-3g",
             mapping: str = "jordan_wigner",
             ansatz_type: str = "uccsd",
             optimizer: str = "Adam",
             maxiter: int = 200,
             use_prev_params: bool = True,
             compute_classical: bool = True,
             callback: Callable[[float, float], None] | None = None) -> PESData:
    """Scan the potential energy surface by running VQE at each bond length.

    For each bond length R in bond_lengths:
      1. Build molecule: mol_spec = molecule_factory(R); set mol_spec.basis = basis.
      2. Run HF via classical_hf.run_hartree_fock -> store HF energy.
      3. Build qubit Hamiltonian via hamiltonian_builder.build_hamiltonian
         using HF one/two-electron integrals, n_electrons, energy_nuc, and mapping.
      4. Build ansatz:
         - "uccsd": ansatz.build_uccsd_ansatz(n_qubits, n_electrons)
         - otherwise: ansatz.build_hardware_efficient_ansatz(n_qubits, n_layers=2)
      5. Run VQE via vqe_solver.run_vqe. If use_prev_params and previous
         optimal params exist with matching dimension, pass as init_params.
      6. Compute FCI energy via classical_hf.get_fci_energy.
         If compute_classical, also compute CISD and CCSD(T).
      7. Call callback(R, vqe_energy) if provided.

    Pre-allocate np.zeros arrays of shape (n_points,). CISD/CCSD(T) arrays
    are only created when compute_classical is True.

    Gotchas:
      - At large bond lengths, RHF may not converge. Consider UHF fallback.
      - VQE param count can change if active space changes; only warm-start
        when dimensions match.

    Args:
        molecule_factory: Callable(bond_length) -> MoleculeSpec.
        bond_lengths: Array of bond lengths to scan (Angstroms).
        basis: Basis set name (default "sto-3g").
        mapping: Fermion-to-qubit mapping (default "jordan_wigner").
        ansatz_type: "uccsd" or "hardware_efficient".
        optimizer: Optimizer for VQE (default "Adam").
        maxiter: Max VQE iterations per point.
        use_prev_params: Warm-start from previous geometry's optimal params.
        compute_classical: Also compute CISD and CCSD(T) energies.
        callback: Optional callable(R, vqe_energy).

    Returns:
        PESData with energies from all methods at each bond length.
    """
    ...
