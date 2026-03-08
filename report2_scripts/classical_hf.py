from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import ao2mo, gto, scf, fci, ci, cc


@dataclass
class HFResult:
    """Results from a Hartree-Fock calculation."""
    energy_hf: float
    energy_nuc: float
    mo_coefficients: np.ndarray
    mo_energies: np.ndarray
    one_electron_integrals: np.ndarray
    two_electron_integrals: np.ndarray
    overlap_matrix: np.ndarray
    n_electrons: int
    n_orbitals: int
    ao_labels: list[str]


def run_hartree_fock(mol_spec: MoleculeSpec) -> HFResult:
    """Run restricted Hartree-Fock via PySCF and extract molecular integrals.

    Steps:
      1. Build PySCF Mole: atom string "symbol x y z; ...",
         basis, charge, spin = multiplicity - 1, unit = "Angstrom".
      2. Run RHF via scf.RHF(mol).kernel().
      3. Extract MO coefficients (n_ao, n_mo) and orbital energies (n_mo,).
      4. Transform one-electron integrals to MO basis:
         h1_mo = C^T @ get_hcore() @ C.
      5. Transform two-electron integrals to MO basis via ao2mo.full(mol, C),
         then restore to 4-index tensor with ao2mo.restore(1, h2_mo, n_mo).
      6. Extract AO overlap matrix via mol.intor("int1e_ovlp").

    Args:
        mol_spec: Molecular system specification.

    Returns:
        HFResult with HF energy, nuclear repulsion, MO coefficients,
        orbital energies, transformed integrals, overlap, electron/orbital counts,
        and AO labels.
    """
    ...


def get_fci_energy(mol_spec: MoleculeSpec) -> float:
    """Compute the Full Configuration Interaction (FCI) energy.

    Builds the PySCF molecule via _build_pyscf_mol, runs RHF as a
    starting point, then runs FCI solver (fci.FCI(mf).kernel()).

    Args:
        mol_spec: Molecular system specification.

    Returns:
        FCI total energy in Hartrees.
    """
    ...


def get_cisd_energy(mol_spec: MoleculeSpec) -> float:
    """Compute the CISD (Configuration Interaction Singles & Doubles) energy.

    Builds PySCF molecule, runs RHF, then CISD via ci.CISD(mf).run().

    Args:
        mol_spec: Molecular system specification.

    Returns:
        CISD total energy in Hartrees.
    """
    ...


def get_ccsd_t_energy(mol_spec: MoleculeSpec) -> float:
    """Compute the CCSD(T) energy.

    Builds PySCF molecule, runs RHF, then CCSD via cc.CCSD(mf).run().
    The perturbative triples correction is added: e_tot + mycc.ccsd_t().

    Args:
        mol_spec: Molecular system specification.

    Returns:
        CCSD(T) total energy in Hartrees.
    """
    ...


def _build_pyscf_mol(mol_spec: MoleculeSpec) -> gto.Mole:
    """Build and return a PySCF Mole object from a MoleculeSpec.

    Constructs atom string in PySCF format ("symbol x y z; ..."),
    sets basis, charge, spin (multiplicity - 1), unit = "Angstrom",
    and calls mol.build().

    Args:
        mol_spec: Molecular system specification.

    Returns:
        Built PySCF Mole object ready for calculations.
    """
    ...
