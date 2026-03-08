from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import gto


@dataclass
class OrbitalVolumeData:
    """Volumetric data for a molecular orbital on a 3D grid."""
    values: np.ndarray
    grid_origin: tuple[float, float, float]
    grid_spacing: tuple[float, float, float]
    grid_shape: tuple[int, int, int]
    orbital_index: int


def evaluate_orbital_on_grid(mo_coefficients: np.ndarray,
                              orbital_index: int,
                              atom_coords: np.ndarray,
                              atom_symbols: list[str],
                              basis_name: str,
                              grid_extent: float = 5.0,
                              grid_spacing: float = 0.1) -> OrbitalVolumeData:
    """Evaluate a single molecular orbital on a uniform 3D grid.

    Steps:
      1. Build a PySCF Mole (unit="Bohr") from atom_coords and atom_symbols.
         Atom string format: "symbol x y z; ...".
      2. Compute grid centre as np.mean(atom_coords, axis=0).
      3. Create 1D arrays x, y, z via np.arange from (centre - extent) to
         (centre + extent) with the given spacing.
      4. Build full 3D grid: np.meshgrid(x, y, z, indexing='ij'),
         then flatten to (N, 3) coordinate array via np.stack + ravel.
      5. Evaluate all AOs at grid points: mol.eval_gto("GTOval_cart", coords)
         -> shape (N, n_ao).
      6. Contract with MO coefficients:
         mo_values = ao_values @ mo_coefficients[:, orbital_index] -> shape (N,).
      7. Reshape to 3D: volume = mo_values.reshape(nx, ny, nz).

    Optimisation: PySCF's eval_gto exploits Gaussian sparsity. For large
    grids, process in batches of ~10000 points to manage memory.

    Gotcha: PySCF uses Bohr internally. If atom_coords are in Angstroms,
    convert: coords_bohr = coords_angstrom * 1.8897259886.

    Args:
        mo_coefficients: MO coefficient matrix, shape (n_ao, n_mo).
        orbital_index: Which MO to evaluate (0-indexed).
        atom_coords: Atomic coordinates, shape (n_atoms, 3), in Bohr.
        atom_symbols: List of element symbols (e.g., ["H", "H"]).
        basis_name: Basis set name for PySCF (e.g., "sto-3g").
        grid_extent: Half-width of cubic grid in Bohr (default 5.0).
        grid_spacing: Distance between grid points in Bohr (default 0.1).

    Returns:
        OrbitalVolumeData with the 3D volume and grid metadata.
    """
    ...


def evaluate_electron_density(mo_coefficients: np.ndarray,
                               n_occupied: int,
                               atom_coords: np.ndarray,
                               atom_symbols: list[str],
                               basis_name: str,
                               grid_extent: float = 5.0,
                               grid_spacing: float = 0.1) -> OrbitalVolumeData:
    """Compute the total electron density on a 3D grid.

    Electron density for RHF: rho(r) = sum_{i=0}^{n_occ-1} 2 * |phi_i(r)|^2
    (factor of 2 for doubly occupied spatial orbitals).

    Iterates over occupied orbitals [0, n_occupied), calling
    evaluate_orbital_on_grid for each, and accumulates 2 * values**2.

    The returned OrbitalVolumeData uses orbital_index = -1 to indicate
    this is a density (not a single orbital).

    Args:
        mo_coefficients: MO coefficient matrix, shape (n_ao, n_mo).
        n_occupied: Number of occupied spatial orbitals.
        atom_coords: Atomic coordinates, shape (n_atoms, 3), in Bohr.
        atom_symbols: List of element symbols.
        basis_name: Basis set name for PySCF.
        grid_extent: Half-width of cubic grid in Bohr.
        grid_spacing: Distance between grid points in Bohr.

    Returns:
        OrbitalVolumeData with electron density volume (orbital_index=-1).
    """
    ...
