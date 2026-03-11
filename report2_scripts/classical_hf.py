from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import ao2mo, gto, scf, fci, ci, cc

from .molecule import MoleculeSpec


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
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol)
    energy_hf = mf.kernel()
    energy_nuc = mol.energy_nuc()

    mo_coefficients = mf.mo_coeff
    mo_energies = mf.mo_energy

    h1_ao = mf.get_hcore()
    h1_mo = mo_coefficients.T @ h1_ao @ mo_coefficients

    # Full 4-index MO integrals (i j | k l)
    eri_mo = ao2mo.kernel(mol, mf.mo_coeff)  # returns compressed by default

    # Restore to full 4D array
    nmo = mf.mo_coeff.shape[1]
    eri_mo_full = ao2mo.restore(1, eri_mo, nmo)  # shape: (nmo, nmo, nmo, nmo)

    overlap_matrix=mol.intor("int1e_ovlp")

    return HFResult(
        energy_hf=energy_hf,
        energy_nuc=energy_nuc,
        mo_coefficients=mo_coefficients,
        mo_energies=mo_energies,
        one_electron_integrals=h1_mo,
        two_electron_integrals=eri_mo_full,
        overlap_matrix=overlap_matrix,
        n_electrons=mol.nelectron,
        n_orbitals=nmo,
        ao_labels=mol.ao_labels())



def get_fci_energy(mol_spec: MoleculeSpec) -> float:
    """Compute the Full Configuration Interaction (FCI) energy.

    Builds the PySCF molecule via _build_pyscf_mol, runs RHF as a
    starting point, then runs FCI solver (fci.FCI(mf).kernel()).

    Args:
        mol_spec: Molecular system specification.

    Returns:
        FCI total energy in Hartrees.
    """
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol)
    mf.kernel()
    fcisolver = fci.FCI(mf)
    energy_fci, _ = fcisolver.kernel()
    return energy_fci


def get_cisd_energy(mol_spec: MoleculeSpec) -> float:
    """Compute the CISD (Configuration Interaction Singles & Doubles) energy.

    Builds PySCF molecule, runs RHF, then CISD via ci.CISD(mf).run().

    Args:
        mol_spec: Molecular system specification.

    Returns:
        CISD total energy in Hartrees.
    """
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol)
    mf.kernel()
    cisolver = ci.CISD(mf)
    cisolver.run()
    return cisolver.e_tot


def get_ccsd_t_energy(mol_spec: MoleculeSpec) -> float:
    """Compute the CCSD(T) energy.

    Builds PySCF molecule, runs RHF, then CCSD via cc.CCSD(mf).run().
    The perturbative triples correction is added: e_tot + mycc.ccsd_t().

    Args:
        mol_spec: Molecular system specification.

    Returns:
        CCSD(T) total energy in Hartrees.
    """
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol)
    mf.kernel()
    cc_solver = cc.CCSD(mf)
    cc_solver.run()
    energy_ccsd_t = cc_solver.e_tot + cc_solver.ccsd_t()
    return energy_ccsd_t


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
    atom_str = "; ".join(
        f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}"
        for atom in mol_spec.atoms)
    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = mol_spec.basis
    mol.charge = mol_spec.charge
    mol.spin = mol_spec.multiplicity - 1
    mol.unit = "Angstrom"
    mol.build()
    return mol
