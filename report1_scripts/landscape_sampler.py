"""Energy landscape sampling over the (γ, β) parameter grid."""

from __future__ import annotations

import numpy as np
import pennylane as qml
from dataclasses import dataclass

from qaoa_circuit import qaoa_ansatz


@dataclass
class LandscapeData:
    """Container for a sampled 2-D energy landscape.

    Attributes:
        gamma_values: 1-D array of γ grid points.
        beta_values: 1-D array of β grid points.
        energies: 2-D array of shape ``(n_gamma, n_beta)`` with ⟨H_C⟩ values.
        approximation_ratios: ``energies / C_max`` (same shape).
    """

    gamma_values: np.ndarray
    beta_values: np.ndarray
    energies: np.ndarray
    approximation_ratios: np.ndarray


def sample_landscape(
    qnode: qml.QNode,
    p: int,
    gamma_range: tuple[float, float] = (0, np.pi),
    beta_range: tuple[float, float] = (0, np.pi / 2),
    n_gamma: int = 50,
    n_beta: int = 50,
    fixed_params: np.ndarray | None = None,
    c_max: float = 1.0,
) -> LandscapeData:
    """Evaluate ⟨H_C⟩ on an evenly-spaced (γ, β) grid.

    For depth ``p = 1`` the full parameter vector is simply ``[γ, β]``.
    For ``p > 1``, *fixed_params* supplies the remaining angles and only
    the first γ (index 0) and first β (index p) are swept.

    Args:
        qnode: QAOA QNode (expects a flat params array of length ``2·p``).
        p: QAOA depth.
        gamma_range: ``(γ_min, γ_max)`` sweep limits.
        beta_range: ``(β_min, β_max)`` sweep limits.
        n_gamma: Number of γ grid points.
        n_beta: Number of β grid points.
        fixed_params: Base parameter vector for ``p > 1`` (shape ``(2·p,)``).
        c_max: Exact MaxCut value, used to compute approximation ratios.

    Returns:
        ``LandscapeData`` holding the grid coordinates and energy values.
    """
    gamma_values = np.linspace(gamma_range[0], gamma_range[1], n_gamma)
    beta_values = np.linspace(beta_range[0], beta_range[1], n_beta)

    energies = np.zeros((n_gamma, n_beta))
    for i, gamma in enumerate(gamma_values):
        for j, beta in enumerate(beta_values):
            params = np.zeros(2 * p)
            params[0] = gamma
            params[p] = beta
            if fixed_params is not None:
                params += fixed_params
            energies[i, j] = qnode(params)
    approximation_ratios = energies / c_max
    return LandscapeData(
        gamma_values=gamma_values,
        beta_values=beta_values,
        energies=energies,
        approximation_ratios=approximation_ratios
    )



def sample_landscape_parallel(
    qnode: qml.QNode,
    p: int,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
    n_gamma: int,
    n_beta: int,
    c_max: float = 1.0,
) -> LandscapeData:
    """Vectorised landscape sampling (delegates to sample_landscape for now)."""
    return sample_landscape(
        qnode, p, gamma_range, beta_range, n_gamma, n_beta, c_max=c_max,
    )
