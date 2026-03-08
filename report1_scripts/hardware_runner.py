"""Quantum hardware execution via IBM and IQM backends."""

from __future__ import annotations
import time

import numpy as np
import pennylane as qml
from dataclasses import dataclass
from qiskit_ibm_runtime import QiskitRuntimeService
import collections
import networkx as nx


@dataclass
class HardwareResult:
    """Container for a single hardware execution result.

    Attributes:
        expectation_value: Measured ⟨H_C⟩ (or sample-estimated energy).
        raw_counts: Bitstring histogram ``{"0110": 42, ...}`` (empty if
            using ``qml.expval`` mode).
        shots: Number of measurement shots.
        backend_name: Identifier string of the backend that was used.
        execution_time: Wall-clock time in seconds.
    """

    expectation_value: float
    raw_counts: dict[str, int]
    shots: int
    backend_name: str
    execution_time: float


def create_ibm_device(
    n_qubits: int,
    backend_name: str = "ibm_brisbane",
    shots: int = 8192,
    ibm_token: str | None = None,
) -> qml.Device:
    """Create a PennyLane device targeting an IBM Quantum backend.

    Uses the ``qiskit.ibmq`` PennyLane-Qiskit plugin device.
    The token is read from *ibm_token* or the ``IBMQ_TOKEN`` environment
    variable.  Transpilation is handled internally by Qiskit.

    Args:
        n_qubits: Number of wires.
        backend_name: IBM backend identifier (e.g. ``"ibm_brisbane"``).
        shots: Number of measurement shots per circuit execution.
        ibm_token: IBM Quantum API token (optional if env var is set).

    Returns:
        PennyLane ``qiskit.ibmq`` device.
    """
    return qml.device(
        "qiskit.ibmq",
        wires=n_qubits,
        backend=backend_name,
        token=ibm_token,
        shots=shots,
    )


#def create_iqm_device(
#    n_qubits: int,
#    backend_url: str,
#    shots: int = 8192,
#) -> qml.Device:
#    """Create a PennyLane device targeting an IQM backend.
#
#    Uses the ``iqm.iqm_device`` PennyLane-IQM plugin.
#
#    Args:
#        n_qubits: Number of wires.
#        backend_url: IQM backend endpoint URL.
#        shots: Number of measurement shots.
#
#    Returns:
#        PennyLane ``iqm.iqm_device`` device.
#    """
#    ...


def execute_on_hardware(
    qnode: qml.QNode,
    params: np.ndarray,
    cost_hamiltonian: qml.Hamiltonian,
    mode: str = "expectation",
) -> HardwareResult:
    """Execute a QAOA circuit on hardware and return the energy estimate.

    Two modes of operation:
      - **Expectation mode** (default): call ``qnode(params)`` directly
        and let PennyLane decompose the Hamiltonian into Pauli measurements.
      - **Sample mode** (alternative): use ``qml.sample()`` to collect
        bitstrings, convert to Ising spins ``s_i = 1 − 2·z_i``, and
        compute the cut energy per sample.

    The function records wall-clock execution time.

    Args:
        qnode: Hardware-backed QAOA QNode.
        params: Variational parameters (shape ``(2·p,)``).
        cost_hamiltonian: Cost Hamiltonian (for sample-mode energy
            computation and metadata).

    Returns:
        ``HardwareResult`` with the measured expectation value and metadata.
    """

    start_time = time.perf_counter()
    if mode == "expectation":
        # Directly evaluate the expectation value via PennyLane's built-in
        # measurement handling.
        energy = qnode(params)
        counts = {}
    elif mode == "sample":
        # Collect raw bitstring samples and compute the energy distribution.
        samples = qnode(params)
        counts = {}
        for sample in samples:
            bitstring = "".join(str(bit) for bit in sample)
            counts[bitstring] = counts.get(bitstring, 0) + 1
        energy = sum(
            (cost_hamiltonian.evaluate(sample) * count for sample, count in counts.items())
        ) / sum(counts.values())
    else:
        raise ValueError(f"Unsupported mode: {mode}")
    execution_time = time.perf_counter() - start_time
    return HardwareResult(
        expectation_value=energy,
        raw_counts=counts,
        shots=sum(counts.values()) if counts else 0,
        backend_name=qnode.device.short_name,
        execution_time=execution_time,
    )

        


def select_optimal_qubits(
    backend_name: str,
    n_qubits: int,
    ibm_token: str | None = None,
) -> list[int]:
    # ── 1. Load backend properties ────────────────────────────────────────
    service  = QiskitRuntimeService(channel="ibm_quantum", token=ibm_token)
    backend  = service.backend(backend_name)
    props    = backend.properties()
    n_phys   = len(props.qubits)

    def get_qubit_param(q: int, name: str) -> float:
        try:
            return props.qubit_property(q, name)[0].value
        except Exception:
            return 0.0

    # ── 2. Build networkx graph ───────────────────────────────────────────
    G = nx.Graph()

    for q in range(n_phys):
        t1 = get_qubit_param(q, "T1")
        t2 = get_qubit_param(q, "T2")
        G.add_node(q, quality=t1 * t2)

    coupling_map = backend.configuration().coupling_map
    for i, j in coupling_map:
        if i >= j:          # undirected — skip duplicates
            continue
        try:
            cx_error = props.gate_error("cx", [i, j])
        except Exception:
            cx_error = 1.0  # worst-case weight if data missing
        if cx_error > 0:
            G.add_edge(i, j, weight=1.0 / cx_error)

    # ── 3. Enumerate connected subgraphs of size n_qubits via BFS ─────────
    def bfs_subgraphs(start: int) -> list[frozenset[int]]:
        """All connected subsets of size n_qubits reachable from start."""
        found = []
        # Each queue entry is a frozenset of nodes in the current subgraph
        queue = collections.deque([frozenset([start])])
        seen  = {frozenset([start])}

        while queue:
            current = queue.popleft()
            if len(current) == n_qubits:
                found.append(current)
                continue
            # Expand: add any neighbour of any node already in the set
            for node in current:
                for neighbour in G.neighbors(node):
                    if neighbour not in current:
                        candidate = current | {neighbour}
                        if candidate not in seen:
                            seen.add(candidate)
                            queue.append(candidate)
        return found

    all_subgraphs: set[frozenset[int]] = set()
    for start in G.nodes:
        all_subgraphs.update(bfs_subgraphs(start))

    if not all_subgraphs:
        raise ValueError(
            f"No connected subgraph of size {n_qubits} found on {backend_name}."
        )

    # ── 4. Score each candidate ───────────────────────────────────────────
    def score(nodes: frozenset[int]) -> float:
        edge_score = sum(
            G[u][v]["weight"]
            for u, v in G.subgraph(nodes).edges()
        )
        node_score = sum(G.nodes[q]["quality"] for q in nodes)
        return edge_score + 0.1 * node_score

    # ── 5. Return best subgraph ───────────────────────────────────────────
    best = max(all_subgraphs, key=score)
    return sorted(best)