from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from .classical_hf import HFResult, MoleculeSpec, run_hartree_fock, get_fci_energy, get_cisd_energy, get_ccsd_t_energy
from .hamiltonian_builder import QubitHamiltonian, build_hamiltonian
from .ansatz import build_uccsd_ansatz, build_hardware_efficient_ansatz
from .vqe_solver import run_vqe


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


def scan_pes(molecule_factory: Callable[[float], MoleculeSpec],
             bond_lengths: np.ndarray,
             basis: str = "sto-3g",
             mapping: str = "jordan_wigner",
             ansatz_type: str = "uccsd",
             optimizer: str = "Adam",
             maxiter: int = 200,
             use_prev_params: bool = True,
             compute_classical: bool = True,
             callback: Callable[[float, float], None] | None = None,
             **kwargs) -> PESData:
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
    n_points = len(bond_lengths)
    
    # Pre-allocate result arrays
    energies_vqe = np.zeros(n_points)
    energies_hf = np.zeros(n_points)
    energies_fci = np.zeros(n_points)
    energies_cisd = np.zeros(n_points) if compute_classical else None
    energies_ccsd_t = np.zeros(n_points) if compute_classical else None
    vqe_params_list = []
    
    prev_optimal_params = None
    
    for idx, R in enumerate(bond_lengths):
        print(f'[{idx+1}/{n_points}] Scanning bond length: {R:.3f} Å')
        
        # 1. Build molecule
        mol_spec = molecule_factory(R)
        mol_spec.basis = basis
        
        # 2. Run Hartree-Fock
        hf_result = run_hartree_fock(mol_spec)
        energies_hf[idx] = hf_result.energy_hf
        
        # 3. Build qubit Hamiltonian
        qubit_hamiltonian = build_hamiltonian(
            one_electron_integrals=hf_result.one_electron_integrals,
            two_electron_integrals=hf_result.two_electron_integrals,
            n_electrons=hf_result.n_electrons,
            nuclear_repulsion=hf_result.energy_nuc,
            mapping=mapping
        )
        
        # 4. Build ansatz
        if ansatz_type == "uccsd":
            ansatz_fn, n_params, _ = build_uccsd_ansatz(
                n_qubits=qubit_hamiltonian.n_qubits,
                n_electrons=hf_result.n_electrons
            )
        else:
            ansatz_fn, n_params = build_hardware_efficient_ansatz(
                n_qubits=qubit_hamiltonian.n_qubits,
                n_layers=2
            )
        
        # 5. Initialize parameters (warm-start if possible)
        init_params = None
        if use_prev_params and prev_optimal_params is not None:
            if len(prev_optimal_params) == n_params:
                print('  Using warm-start from previous geometry.')
                init_params = prev_optimal_params
        
        # 6. Run VQE
        vqe_result = run_vqe(
            hamiltonian=qubit_hamiltonian.hamiltonian,
            ansatz_fn=ansatz_fn,
            n_params=n_params,
            n_qubits=qubit_hamiltonian.n_qubits,
            init_params=init_params,
            optimizer=optimizer,
            learning_rate=0.1,
            maxiter=maxiter,
            tol=1e-6,
            device_name="default.qubit",
            diff_method="backprop"
        )
        
        energies_vqe[idx] = vqe_result.optimal_energy
        vqe_params_list.append(vqe_result.optimal_params)
        prev_optimal_params = vqe_result.optimal_params
        
        # 7. Compute classical reference energies
        energies_fci[idx] = get_fci_energy(mol_spec)
        
        if compute_classical:
            energies_cisd[idx] = get_cisd_energy(mol_spec)
            energies_ccsd_t[idx] = get_ccsd_t_energy(mol_spec)
        
        # 8. Callback
        if callback:
            callback(R, vqe_result.optimal_energy)
        
        print(f'  HF: {energies_hf[idx]:.6f}, VQE: {energies_vqe[idx]:.6f}, FCI: {energies_fci[idx]:.6f} Ha')
        print(f'  |VQE - FCI| = {abs(energies_vqe[idx] - energies_fci[idx]):.2e} Ha')
    
    return PESData(
        bond_lengths=bond_lengths,
        energies_vqe=energies_vqe,
        energies_hf=energies_hf,
        energies_fci=energies_fci,
        energies_cisd=energies_cisd,
        energies_ccsd_t=energies_ccsd_t,
        vqe_params=vqe_params_list
    )