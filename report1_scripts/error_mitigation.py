"""Zero-noise extrapolation (ZNE) via unitary folding."""

from __future__ import annotations

from typing import Callable

import numpy as np
import pennylane as qml


def global_fold(
    circuit_fn: Callable,
    scale_factor: float,
) -> Callable:
    """Return a new circuit function with globally-folded unitary.

    Global unitary folding replaces U with U (U† U)^n U_partial where:

        scale_factor = 1 + 2n + 2·(partial_fraction)

    For integer odd scale factors (1, 3, 5, …): ``n = (s − 1) / 2``,
    giving U followed by *n* copies of U†U.

    For non-integer scale factors, additionally fold the last
    ``frac · total_gates`` gates (requires tape introspection).

    Implementation approach:
      1. Apply the original ``circuit_fn``.
      2. Append *n* repetitions of ``qml.adjoint(circuit_fn)`` then
         ``circuit_fn``.
      3. For fractional parts, fold the last *k* gates where
         ``k = round(frac × total_gates)``.

    Args:
        circuit_fn: A callable that applies PennyLane quantum operations
            (no return value).
        scale_factor: Noise amplification factor (≥ 1.0).

    Returns:
        A new callable with the same signature that applies the folded
        circuit.
    """
    n = np.floor((scale_factor - 1) / 2)
    frac = (scale_factor - 1) / 2 - n
    n = int(n)
    def folded_circuit(*args, **kwargs):
        circuit_fn(*args, **kwargs)
        for _ in range(n):
            qml.adjoint(circuit_fn)(*args, **kwargs)
            circuit_fn(*args, **kwargs)
        if frac > 0:
            with qml.tape.QuantumTape() as tape:
                circuit_fn(*args, **kwargs)
            total_gates = len(tape.operations)
            k = int(round(frac * total_gates))
            if k > 0:
                last_ops = tape.operations[-k:]
                for op in reversed(last_ops):
                    qml.adjoint(op)
                for op in last_ops:
                    qml.apply(op)

    return folded_circuit


def zne_extrapolate(
    qnode: qml.QNode,
    params: np.ndarray,
    scale_factors: list[float] | None = None,
    extrapolation: str = "linear",
) -> tuple[float, dict]:
    """Perform zero-noise extrapolation on a noisy QNode.

    Evaluates the QNode at each scale factor via global unitary folding
    and extrapolates to the zero-noise limit.

    Args:
        qnode: Noisy QAOA QNode.
        params: Parameter vector.
        scale_factors: Noise scale factors (default ``[1.0, 2.0, 3.0]``).
        extrapolation: ``"linear"`` or ``"richardson"``.

    Returns:
        A 2-tuple ``(mitigated_value, details)`` where *details* is a
        dict with keys ``"noisy_values"``, ``"scale_factors"``,
        ``"extrapolated_value"``.
    """
    if scale_factors is None:
        scale_factors = [1.0, 2.0, 3.0]

    noisy_values = []
    for s in scale_factors:
        if s == 1.0:
            val = float(qnode(params))
        else:
            # Use PennyLane's fold_global transform
            folded_qnode = qml.noise.fold_global(qnode, scale_factor=s)
            val = float(folded_qnode(params))
        noisy_values.append(val)

    noisy_values = np.array(noisy_values)
    sf = np.array(scale_factors)

    if extrapolation == "linear":
        mitigated_value = linear_extrapolation(sf, noisy_values)
    elif extrapolation == "richardson":
        mitigated_value = richardson_extrapolation(sf, noisy_values)
    else:
        raise ValueError(f"Unsupported extrapolation method: {extrapolation}")

    details = {
        "noisy_values": noisy_values,
        "scale_factors": sf,
        "extrapolated_value": mitigated_value,
    }
    return mitigated_value, details



def richardson_extrapolation(
    scale_factors: np.ndarray,
    noisy_values: np.ndarray,
) -> float:
    """Richardson extrapolation to zero noise.

    Fits a polynomial of degree ``M − 1`` through the *M* data points
    ``(s_i, E(s_i))`` and returns ``E(0)``.

    Mathematically, solve the Vandermonde system:

        V · a = noisy_values

    where ``V_{i,j} = s_i^j`` (increasing powers) and the zero-noise
    estimate is ``a[0]``.

    Args:
        scale_factors: 1-D array of scale factors (length M).
        noisy_values: 1-D array of noisy expectation values (length M).

    Returns:
        Extrapolated zero-noise expectation value.
    """
    V = np.vander(scale_factors, increasing=True)
    a = np.linalg.solve(V, noisy_values)
    return a[0]
    


def linear_extrapolation(
    scale_factors: np.ndarray,
    noisy_values: np.ndarray,
) -> float:
    """Linear extrapolation to zero noise.

    Fits ``E(s) = a + b·s`` via ``np.polyfit`` and returns the intercept
    ``a = E(0)``.

    Args:
        scale_factors: 1-D array of scale factors.
        noisy_values: 1-D array of noisy expectation values.

    Returns:
        Extrapolated zero-noise expectation value.
    """
    coeffs = np.polyfit(scale_factors, noisy_values, 1)
    return coeffs[1]
