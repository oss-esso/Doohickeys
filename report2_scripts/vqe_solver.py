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
    device = qml.device(device_name, wires=n_qubits)
    @qml.qnode(device, diff_method=diff_method)

    def cost_fn(params):
        ansatz_fn(params)
        return qml.expval(hamiltonian)
    
    if init_params is None:
        init_params = np.random.normal(0, 0.01, size=n_params)

    trajectory = []
    gradient_norms = []

    if optimizer == 'Adam':
        opt = qml.AdamOptimizer(stepsize=learning_rate)
        params = init_params.copy()
        energy = cost_fn(params)
        for n in range(maxiter):
            params, energy = opt.step_and_cost(cost_fn, params)
            grad = qml.grad(cost_fn)(params)
            grad_norm = np.linalg.norm(grad)
            gradient_norms.append(grad_norm)
            trajectory.append((params.copy(), energy))
            if callback:
                callback(n, params, energy)
            if n > 0 and abs(trajectory[-2][1] - energy) < tol:
                break
    elif optimizer == 'GradientDescent':
        opt = qml.GradientDescentOptimizer(stepsize=learning_rate)
        params = init_params.copy()
        energy = cost_fn(params)
        for n in range(maxiter):
            params, energy = opt.step_and_cost(cost_fn, params)
            grad = qml.grad(cost_fn)(params)
            grad_norm = np.linalg.norm(grad)
            gradient_norms.append(grad_norm)
            trajectory.append((params.copy(), energy))
            if callback:
                callback(n, params, energy)
            if n > 0 and abs(trajectory[-2][1] - energy) < tol:
                break
    elif optimizer == 'COBYLA':
        # COBYLA trajectory tracking via callback
        cobyla_trajectory = []
        
        def scipy_cost_fn(params):
            energy = float(cost_fn(params))
            cobyla_trajectory.append((params.copy(), energy))
            if callback:
                callback(len(cobyla_trajectory) - 1, params, energy)
            return energy
        
        res = sp_opt.minimize(scipy_cost_fn, init_params, method='COBYLA', 
                              options={'maxiter': maxiter, 'rhobeg': 0.5})
        params = res.x
        energy = res.fun
        trajectory = cobyla_trajectory if cobyla_trajectory else [(params.copy(), energy)]
        gradient_norms = [0.0] * len(trajectory)  # COBYLA is gradient-free
    elif optimizer == 'SPSA':
        a0 = 0.1
        c0 = 0.1
        alpha = 0.602
        gamma = 0.101
        A = maxiter * 0.1
        params = init_params.copy()
        energy = cost_fn(params)
        for n in range(maxiter):
            ak = a0 / (n + 1 + A) ** alpha
            ck = c0 / (n + 1) ** gamma
            delta = np.random.choice([-1, 1], size=n_params)
            energy_plus = cost_fn(params + ck * delta)
            energy_minus = cost_fn(params - ck * delta)
            grad_estimate = (energy_plus - energy_minus) / (2 * ck * delta)
            params -= ak * grad_estimate
            energy = cost_fn(params)
            gradient_norms.append(np.linalg.norm(grad_estimate))
            trajectory.append((params.copy(), energy))
            if callback:
                callback(n, params, energy)
            if n > 0 and abs(trajectory[-2][1] - energy) < tol:
                break
            
    else:
        raise ValueError(f"Unsupported optimizer: {optimizer}") 
    
    optimal_energy = min(trajectory, key=lambda x: x[1])[1]
    optimal_params = min(trajectory, key=lambda x: x[1])[0]
    n_iterations = len(trajectory)
    converged = n_iterations < maxiter

    return VQEResult(
        optimal_energy=optimal_energy,
        optimal_params=optimal_params,
        trajectory=trajectory,
        n_iterations=n_iterations,
        converged=converged,
        gradient_norms=gradient_norms)


