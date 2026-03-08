from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pennylane as qml
import scipy.optimize as sp_opt


@dataclass
class VQEResult:
    """Container for VQE optimisation results."""
    optimal_energy: float
    optimal_params: np.ndarray
    trajectory: list[tuple[np.ndarray, float]]
    n_iterations: int
    converged: bool
    gradient_norms: list[float]


def run_vqe(hamiltonian: object, ansatz_fn: Callable[[np.ndarray], None],
            n_params: int, n_qubits: int,
            init_params: np.ndarray | None = None,
            optimizer: str = "Adam", learning_rate: float = 0.05,
            maxiter: int = 200, tol: float = 1e-6,
            device_name: str = "default.qubit",
            diff_method: str = "backprop",
            callback: Callable[[int, np.ndarray, float], None] | None = None
            ) -> VQEResult:
    """Run the Variational Quantum Eigensolver to find the ground-state energy.

    Creates a PennyLane device and a QNode cost function that applies ansatz_fn
    then measures qml.expval(hamiltonian).

    If init_params is None, initialise to zeros + small random perturbation
    (np.random.normal(0, 0.01)) to break symmetry.

    Supported optimizers:
      - "Adam": qml.AdamOptimizer(stepsize=learning_rate).
        Uses opt.step_and_cost for param updates, qml.grad for gradient norms.
      - "GradientDescent": qml.GradientDescentOptimizer(stepsize=learning_rate).
        Same interface as Adam.
      - "COBYLA": scipy.optimize.minimize(method="COBYLA", rhobeg=0.5).
        Gradient-free; gradient_norms recorded as 0.0.
      - "SPSA": Simultaneous Perturbation Stochastic Approximation.
        Hyperparams: a0=0.1, c0=0.1, alpha=0.602, gamma=0.101, A=maxiter*0.1.
        Bernoulli ±1 perturbation vectors. Gradient estimate:
          g_hat = (f(p + c_k*delta) - f(p - c_k*delta)) / (2*c_k*delta)

    Convergence: check |ΔE| < tol between consecutive iterations.
    Trajectory records (params_copy, energy) at each step.

    Returns VQEResult with the globally best energy across all iterations
    (argmin over trajectory energies).

    Args:
        hamiltonian: PennyLane Hamiltonian operator.
        ansatz_fn: Callable that applies the ansatz circuit given params.
        n_params: Number of variational parameters.
        n_qubits: Number of qubits for the device.
        init_params: Optional initial parameter vector.
        optimizer: One of "Adam", "GradientDescent", "COBYLA", "SPSA".
        learning_rate: Step size for gradient-based optimizers.
        maxiter: Maximum number of iterations.
        tol: Convergence tolerance on energy change.
        device_name: PennyLane device name (default "default.qubit").
        diff_method: Differentiation method for QNode (default "backprop").
        callback: Optional function called as callback(iteration, params, energy).

    Returns:
        VQEResult with optimal energy, parameters, and optimisation history.
    """
    ...
