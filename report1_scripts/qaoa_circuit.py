"""QAOA circuit construction: cost layer, mixer layer, and full ansatz."""

from __future__ import annotations

import pennylane as qml


def cost_layer(gamma: float, graph_edges: list[tuple[int, int, float]]) -> None:
    """Apply the QAOA cost (phase-separation) layer.

    For each weighted edge ``(u, v, w)``, apply the ZZ-rotation gate:

        exp(−i γ w/2  Z_u ⊗ Z_v)

    implemented via ``qml.IsingZZ(γ·w, wires=[u, v])``.

    Args:
        gamma: Cost-layer variational angle.
        graph_edges: List of ``(u, v, w)`` tuples — qubit indices and weight.
    """
    for u, v, w in graph_edges:
        qml.IsingZZ(-gamma * w, wires=[u, v])


def mixer_layer(beta: float, n_qubits: int) -> None:
    """Apply the QAOA mixer layer.

    For each qubit *i*, apply:

        exp(−i β X_i)  =  RX(2β)

    Args:
        beta: Mixer-layer variational angle.
        n_qubits: Number of qubits in the register.
    """
    for i in range(n_qubits):
        qml.RX(2 * beta, wires=i)


def qaoa_ansatz(
    params: list[float],
    cost_hamiltonian: qml.Hamiltonian,
    n_qubits: int,
    p: int,
) -> None:
    """Prepare the full depth-*p* QAOA state.

    Steps:
      1. Initialise all qubits in |+⟩ (Hadamard on each wire).
      2. For each layer k = 1 … p, apply ``cost_layer(γ_k, edges)``
         followed by ``mixer_layer(β_k, n_qubits)``.

    The parameter vector layout is ``[γ_1, …, γ_p, β_1, …, β_p]``.

    Hint: use ``_extract_edges_from_hamiltonian`` to obtain the edge list
    from *cost_hamiltonian*.

    Args:
        params: Flat array of length 2·p (gammas then betas).
        cost_hamiltonian: PennyLane Hamiltonian built by
            ``graph_generator.build_cost_hamiltonian``.
        n_qubits: Number of qubits.
        p: QAOA depth (number of layers).
    """
    edges = _extract_edges_from_hamiltonian(cost_hamiltonian)
    gamma = params[:p]
    beta = params[p:]
    for i in range(n_qubits):
        qml.H(i)
    
    for layer in range(p):
        cost_layer(gamma[layer], edges)
        mixer_layer(beta[layer], n_qubits)



def build_qaoa_circuit(
    cost_hamiltonian: qml.Hamiltonian,
    n_qubits: int,
    p: int,
    device: qml.Device,
) -> qml.QNode:
    """Return a QNode that evaluates ⟨H_C⟩ for given variational parameters.

    The returned callable has signature ``circuit(params) -> float`` where
    *params* is a flat array of shape ``(2·p,)``.

    Use ``diff_method="best"`` so PennyLane selects the most efficient
    differentiation strategy for the device.

    Args:
        cost_hamiltonian: QAOA cost Hamiltonian.
        n_qubits: Number of qubits.
        p: QAOA depth.
        device: PennyLane device (simulator or hardware).

    Returns:
        ``qml.QNode`` computing ``⟨ψ(γ,β)|H_C|ψ(γ,β)⟩``.
    """
    @qml.qnode(device, diff_method="best")
    def circuit(params: list[float]) -> float:
        qaoa_ansatz(params, cost_hamiltonian, n_qubits, p)
        return qml.expval(cost_hamiltonian)

    return circuit


def _extract_edges_from_hamiltonian(
    H: qml.Hamiltonian,
) -> list[tuple[int, int, float]]:
    """Extract the weighted edge list from a MaxCut Hamiltonian.

    Iterates over the terms of *H*, skips Identity terms, and recovers
    the original edge weight from the ZZ coefficient via ``w = −2·coeff``
    (since each ZZ term was stored as ``−w/2``).

    Args:
        H: PennyLane Hamiltonian produced by ``build_cost_hamiltonian``.

    Returns:
        List of ``(u, v, weight)`` tuples.
    """
    edges = []
    coeffs, ops = H.terms()
    for coeff, op in zip(coeffs, ops):
        if isinstance(op, qml.Identity):
            continue
        wires = op.wires.tolist()
        if len(wires) == 2:
            w = -2 * coeff
            edges.append((wires[0], wires[1], w))
    return edges
