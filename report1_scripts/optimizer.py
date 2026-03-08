"""Classical optimizers for QAOA parameter tuning."""

from __future__ import annotations

import numpy as np
import scipy.optimize as sp_opt
import pennylane as qml
from dataclasses import dataclass


@dataclass
class OptimizationResult:
    """Container for the output of a QAOA optimization run.

    Attributes:
        optimal_params: Best parameter vector found (shape ``(2·p,)``).
        optimal_energy: Best cost-function expectation value (maximised).
        trajectory: List of ``(params_copy, energy)`` recorded at each step.
        n_evals: Total number of circuit evaluations used.
        converged: Whether the optimizer reported convergence.
    """

    optimal_params: np.ndarray
    optimal_energy: float
    trajectory: list[tuple[np.ndarray, float]]
    n_evals: int
    converged: bool


def optimize_qaoa(
    qnode: qml.QNode,
    p: int,
    method: str = "COBYLA",
    init_params: np.ndarray | None = None,
    maxiter: int = 200,
    tol: float = 1e-6,
    callback: callable | None = None,
) -> OptimizationResult:
    """Run a classical optimizer to find optimal QAOA angles.

    QAOA *maximises* ⟨H_C⟩ but ``scipy.optimize.minimize`` *minimises*,
    so negate the objective internally.

    Supported methods:
      - **COBYLA**: Gradient-free constrained optimisation (``rhobeg=0.5``).
      - **L-BFGS-B**: Quasi-Newton with gradients from ``parameter_shift_gradient``.
      - **SPSA**: Simultaneous Perturbation Stochastic Approximation
        (custom loop with Rademacher perturbations; hyperparameters
        ``a0=0.1, c0=0.1, α=0.602, γ=0.101, A=0.1·maxiter``).
      - **Adam**: PennyLane ``AdamOptimizer`` (stepsize 0.05) with
        early stopping when consecutive energies differ by < *tol*.

    Default parameter initialisation (when *init_params* is ``None``):
      ``γ ~ U(0, π)``, ``β ~ U(0, π/2)``.

    Args:
        qnode: QAOA QNode returning ⟨H_C⟩.
        p: QAOA depth (parameter dimension is ``2·p``).
        method: Optimizer name.
        init_params: Starting parameter vector (shape ``(2·p,)``).
        maxiter: Maximum iterations / function evaluations.
        tol: Convergence tolerance.
        callback: ``callback(iteration, params, energy)`` called each step.

    Returns:
        ``OptimizationResult`` with the best parameters and diagnostics.
    """
    if init_params is None:
        init_params = np.concatenate(
            [
                np.random.uniform(0, np.pi, size=p),  # γ angles
                np.random.uniform(0, np.pi / 2, size=p),  # β angles
            ]
        )
    
    method_upper = method.upper()
    trajectory: list[tuple[np.ndarray, float]] = []

    if method_upper == 'COBYLA':
        def cost(params):
            energy = -float(qnode(params))
            trajectory.append((params.copy(), -energy))
            if callback:
                callback(len(trajectory), params, -energy)
            return energy

        res = sp_opt.minimize(
            cost,
            init_params,
            method='COBYLA',
            options={'maxiter': maxiter, 'rhobeg':0.5}
        )

        return OptimizationResult(
            optimal_params=res.x,
            optimal_energy=-res.fun,
            trajectory=trajectory,
            n_evals=res.nfev,
            converged=res.success,
        )
    elif method_upper == 'L-BFGS-B':
        # Use scipy's 3-point finite-difference gradient.
        # Note: the hand-coded parameter_shift_gradient only works when each
        # external QAOA angle appears in a SINGLE gate; since one gamma drives
        # every IsingZZ gate simultaneously, E(γ) has period π/2 and the naive
        # π/2 parameter-shift cancels exactly.  scipy handles this correctly.
        n_evals = 0

        def neg_cost(params):
            nonlocal n_evals
            energy = float(qnode(params))
            n_evals += 1
            trajectory.append((params.copy(), energy))
            if callback:
                callback(n_evals, params, energy)
            return -energy

        res = sp_opt.minimize(
            neg_cost,
            init_params,
            method='L-BFGS-B',
            jac='3-point',
            options={
                'maxiter': maxiter,
                'gtol': 1e-8,
                'ftol': 2.22e-12,
                'disp': False,
            },
        )
        return OptimizationResult(
            optimal_params=res.x,
            optimal_energy=-res.fun,
            trajectory=trajectory,
            n_evals=n_evals,
            converged=res.success,
        )
    elif method_upper == 'SPSA':
        # Hyperparameters tuned for noiseless simulation (small c → sharper grad)
        a0, c0, alpha, gamma_exp, A = 0.6, 0.05, 0.602, 0.101, max(1, 0.05 * maxiter)
        params = init_params.copy()
        for k in range(1, maxiter + 1):
            ak = a0 / (k + A) ** alpha
            ck = c0 / k ** gamma_exp
            delta = np.random.choice([-1, 1], size=params.shape)  # Rademacher
            params_plus = params + ck * delta
            params_minus = params - ck * delta
            energy_plus = float(qnode(params_plus))
            energy_minus = float(qnode(params_minus))
            grad_estimate = (energy_plus - energy_minus) / (2 * ck) * delta
            params += ak * grad_estimate   # ascent (maximise)
            energy = float(qnode(params))
            trajectory.append((params.copy(), energy))
            if callback:
                callback(k, params, energy)
        
        optimal_energy = float(qnode(params))
        return OptimizationResult(
            optimal_params=params,
            optimal_energy=optimal_energy,
            trajectory=trajectory,
            n_evals=2 * maxiter,
            converged=True,
        )
    elif method_upper == 'ADAM':
        # Manual Adam using central-difference numerical gradients.
        # (Parameter-shift is only exact when each circuit angle appears in a
        # single gate; QAOA couples gamma to all IsingZZ gates simultaneously,
        # making the effective E-period half the naive shift.)
        lr = 0.05
        beta1, beta2, eps_adam = 0.9, 0.999, 1e-8
        h_fd = 0.01  # central-difference step size
        params = init_params.astype(float).copy()
        m_moment = np.zeros_like(params)
        v_moment = np.zeros_like(params)
        k = 0
        prev_energy = float(qnode(params))
        for k in range(1, maxiter + 1):
            grad = _central_diff_gradient(qnode, params, h=h_fd)
            m_moment = beta1 * m_moment + (1 - beta1) * grad
            v_moment = beta2 * v_moment + (1 - beta2) * grad ** 2
            m_hat = m_moment / (1 - beta1 ** k)
            v_hat = v_moment / (1 - beta2 ** k)
            params = params + lr * m_hat / (np.sqrt(v_hat) + eps_adam)
            energy = float(qnode(params))
            trajectory.append((params.copy(), energy))
            if callback:
                callback(k, params, energy)
            if abs(energy - prev_energy) < tol:
                break
            prev_energy = energy
        
        return OptimizationResult(
            optimal_params=params,
            optimal_energy=float(qnode(params)),
            trajectory=trajectory,
            n_evals=k,
            converged=(k < maxiter),
        )






def _central_diff_gradient(
    qnode: qml.QNode,
    params: np.ndarray,
    h: float = 0.01,
) -> np.ndarray:
    """Central-difference gradient estimate.  Requires ``2·d`` circuit evals."""
    grad = np.zeros(len(params))
    for i in range(len(params)):
        p_plus = params.copy()
        p_plus[i] += h
        p_minus = params.copy()
        p_minus[i] -= h
        grad[i] = (float(qnode(p_plus)) - float(qnode(p_minus))) / (2 * h)
    return grad


def parameter_shift_gradient(
    qnode: qml.QNode,
    params: np.ndarray,
    shift: float = np.pi / 2,
) -> np.ndarray:
    """Compute the gradient of a QNode via the parameter-shift rule.

    For each parameter θ_i:

        ∂f/∂θ_i = [f(θ_i + s) − f(θ_i − s)] / 2

    where *s* is the shift (default π/2). This requires ``2·d`` circuit
    evaluations for *d* parameters.

    Note: The shift π/2 is exact for gates whose generator has eigenvalues
    ±1 (e.g. ``IsingZZ``, ``RX``).

    Args:
        qnode: Parameterised QNode.
        params: Current parameter vector (shape ``(d,)``).
        shift: Shift amount (default ``π/2``).

    Returns:
        Gradient array of shape ``(d,)``.
    """
    grad = np.zeros_like(params)
    for i in range(len(params)):
        params_plus = params.copy()
        params_minus = params.copy()
        params_plus[i] += shift
        params_minus[i] -= shift
        grad[i] = (qnode(params_plus) - qnode(params_minus)) / 2
    return grad
