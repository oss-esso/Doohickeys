from __future__ import annotations

import os

import numpy as np
import pennylane as qml


def create_ibm_device(n_qubits: int, backend_name: str = "ibm_brisbane",
                      shots: int = 8192,
                      ibm_token: str | None = None) -> qml.Device:
    """Create a PennyLane device connected to an IBM Quantum backend.

    Uses ibm_token if provided, otherwise reads from the IBMQ_TOKEN
    environment variable. Constructs a "qiskit.ibmq" device with the
    given number of wires, backend, and shot count.

    Args:
        n_qubits: Number of qubits (wires).
        backend_name: IBM backend name (default "ibm_brisbane").
        shots: Number of measurement shots (default 8192).
        ibm_token: Optional IBM Quantum API token.

    Returns:
        A PennyLane device configured for the IBM backend.
    """
    if ibm_token is None:
        ibm_token = os.getenv("IBMQ_TOKEN")
        if ibm_token is None:
            raise ValueError("IBM Quantum API token not provided. Set the "
                             "IBMQ_TOKEN environment variable or pass ibm_token.")
    provider = qml.qiskit.IBMProvider(token=ibm_token)
    return qml.qiskit.IBMQDevice(wires=n_qubits, backend=backend_name, provider=provider, shots=shots)



def create_ionq_device(n_qubits: int, backend: str = "ionq_harmony",
                       shots: int = 1024,
                       api_key: str | None = None) -> qml.Device:
    """Create a PennyLane device connected to an IonQ backend.

    Uses api_key if provided, otherwise reads from the IONQ_API_KEY
    environment variable. Constructs an "ionq.simulator" device
    (use "ionq.qpu" for real hardware).

    Args:
        n_qubits: Number of qubits (wires).
        backend: IonQ backend name (default "ionq_harmony").
        shots: Number of measurement shots (default 1024).
        api_key: Optional IonQ API key.

    Returns:
        A PennyLane device configured for the IonQ backend.
    """
    if api_key is None:
        api_key = os.getenv("IONQ_API_KEY")
        if api_key is None:
            raise ValueError("IonQ API key not provided. Set the "
                             "IONQ_API_KEY environment variable or pass api_key.")
    return qml.ionq.IonQDevice(wires=n_qubits, backend=backend, api_key=api_key, shots=shots)


def estimate_shot_budget(hamiltonian: object,
                         target_precision: float = 1.6e-3) -> int:
    """Estimate the number of measurement shots needed for a target precision.

    Uses the variance bound for independent Pauli measurements:
        Var[E] = sum_k |c_k|^2 / S

    To achieve sqrt(Var[E]) <= target_precision:
        S >= sum_k |c_k|^2 / target_precision^2

    This is a conservative upper bound; Pauli grouping (commuting cliques)
    reduces the actual requirement significantly.

    Clamp result to [1_000, 10_000_000].

    Example: For H2 (STO-3G), sum|c_k|^2 ~ 0.3, target = 1.6e-3
    => S ~ 117,000 shots.

    Args:
        hamiltonian: PennyLane Hamiltonian with .coeffs attribute
                     (array of Pauli term coefficients).
        target_precision: Desired energy precision in Hartrees.

    Returns:
        Estimated shot count (int), clamped to [1_000, 10_000_000].
    """
    coeffs = np.array(hamiltonian.coeffs)
    sum_squared = np.sum(np.abs(coeffs)**2)
    shots = int(np.ceil(sum_squared / target_precision**2))
    return max(1_000, min(shots, 10_000_000))
