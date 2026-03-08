from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np


@dataclass
class Atom:
    """A single atom in a molecule."""
    symbol: str
    position: np.ndarray


@dataclass
class MoleculeSpec:
    """Full specification of a molecular system."""
    atoms: list[Atom]
    basis: str = "sto-3g"
    charge: int = 0
    multiplicity: int = 1
    name: str = ""


def create_h2(bond_length: float = 0.735) -> MoleculeSpec:
    """Create an H2 molecule specification along the z-axis.

    Places two hydrogen atoms symmetric about the origin along the z-axis,
    separated by the given bond length. Uses STO-3G basis, neutral charge,
    and singlet multiplicity.

    Args:
        bond_length: H-H distance in Angstroms (default 0.735, near equilibrium).

    Returns:
        MoleculeSpec with two H atoms at (0, 0, ±bond_length/2).
    """
    ...


def create_lih(bond_length: float = 1.595) -> MoleculeSpec:
    """Create a LiH molecule specification along the z-axis.

    Places Li at the origin and H at (0, 0, bond_length).
    Uses STO-3G basis, neutral charge, and singlet multiplicity.

    Args:
        bond_length: Li-H distance in Angstroms (default 1.595, near equilibrium).

    Returns:
        MoleculeSpec with Li at origin and H along z-axis.
    """
    ...


def create_h2o(oh_length: float = 0.957, angle_deg: float = 104.5) -> MoleculeSpec:
    """Create an H2O molecule specification.

    Places oxygen at the origin with two hydrogens at the given O-H distance
    and H-O-H bond angle. The molecule lies in the xz-plane with the oxygen
    at the origin and hydrogens symmetric about the z-axis.

    Convert angle_deg to radians. H positions:
      H1: (+oh_length * sin(angle/2), 0, -oh_length * cos(angle/2))
      H2: (-oh_length * sin(angle/2), 0, -oh_length * cos(angle/2))

    Args:
        oh_length: O-H bond distance in Angstroms (default 0.957).
        angle_deg: H-O-H bond angle in degrees (default 104.5).

    Returns:
        MoleculeSpec for water with STO-3G basis.
    """
    ...


def set_bond_length(mol: MoleculeSpec, atom_idx_a: int,
                    atom_idx_b: int, new_length: float) -> MoleculeSpec:
    """Return a deep copy of mol with the bond between two atoms set to new_length.

    Adjusts both atom positions symmetrically about their midpoint along
    the existing bond axis. The direction vector from atom_a to atom_b is
    preserved (normalised to a unit vector); both atoms are moved to
    ±new_length/2 from the midpoint.

    Args:
        mol: Original molecule specification (not modified).
        atom_idx_a: Index of the first atom in mol.atoms.
        atom_idx_b: Index of the second atom in mol.atoms.
        new_length: Desired bond length in Angstroms.

    Returns:
        New MoleculeSpec with updated atom positions.
    """
    ...
