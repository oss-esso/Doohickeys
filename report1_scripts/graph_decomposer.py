"""Graph decomposition for QAOA on graphs exceeding the qubit threshold.

Strategy: Kernighan-Lin recursive bisection splits G into parts whose node
count fits the qubit budget.  Each part is solved independently with QAOA.
Cross-partition edges are then optimised with a greedy local-search pass.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import networkx as nx
import numpy as np
from networkx.algorithms.community import kernighan_lin_bisection

log = logging.getLogger(__name__)

QUBIT_THRESHOLD_DEFAULT = 12


@dataclass
class Partition:
    """A single induced subgraph with its QAOA result."""

    nodes: list[int]
    subgraph: nx.Graph
    qaoa_result: dict | None = None  # populated after QAOA run
    assignments: dict[int, int] = field(default_factory=dict)  # node → {0,1}


@dataclass
class DecompositionResult:
    """Aggregated result of decomposed-QAOA on a large graph."""

    partitions: list[Partition]
    cross_edges: list[tuple[int, int, float]]
    assignments: dict[int, int]
    total_cut: float
    intra_cut: float
    cross_cut: float
    c_max_upper: float
    approximation_ratio: float
    n_partitions: int
    largest_partition: int


# ── Bisection ─────────────────────────────────────────────────────────────────


def _recursive_bisect(
    graph: nx.Graph,
    threshold: int,
    seed: int,
) -> list[set[int]]:
    """Recursively bisect *graph* until every part has ≤ *threshold* nodes."""
    if graph.number_of_nodes() <= threshold:
        return [set(graph.nodes())]

    try:
        part_a, part_b = kernighan_lin_bisection(
            graph, weight="weight", seed=seed,
        )
    except Exception:
        # Fallback: arbitrary split if KL fails (e.g. disconnected graph)
        nodes = list(graph.nodes())
        mid = len(nodes) // 2
        part_a, part_b = set(nodes[:mid]), set(nodes[mid:])

    result: list[set[int]] = []
    for part in (part_a, part_b):
        sub = graph.subgraph(part)
        result.extend(_recursive_bisect(sub, threshold, seed))
    return result


def decompose_graph(
    graph: nx.Graph,
    threshold: int = QUBIT_THRESHOLD_DEFAULT,
    seed: int = 42,
) -> tuple[list[Partition], list[tuple[int, int, float]]]:
    """Partition *graph* into subgraphs of ≤ *threshold* nodes.

    Uses recursive Kernighan-Lin bisection.  Returns the list of Partition
    objects and the cross-partition edge list ``(u, v, weight)``.

    Args:
        graph: Input weighted graph.
        threshold: Maximum allowed nodes per partition.
        seed: RNG seed for KL bisection.

    Returns:
        ``(partitions, cross_edges)``
    """
    parts = _recursive_bisect(graph, threshold, seed)

    node_to_part: dict[int, int] = {}
    for idx, part in enumerate(parts):
        for node in part:
            node_to_part[node] = idx

    partitions: list[Partition] = []
    for part in parts:
        nodes = sorted(part)
        sub = graph.subgraph(nodes).copy()
        partitions.append(Partition(nodes=nodes, subgraph=sub))

    cross_edges: list[tuple[int, int, float]] = [
        (u, v, d.get("weight", 1.0))
        for u, v, d in graph.edges(data=True)
        if node_to_part[u] != node_to_part[v]
    ]

    log.info(
        "Decomposed: %d nodes → %d partitions (max=%d nodes), %d cross-edges",
        graph.number_of_nodes(), len(partitions),
        max(len(p.nodes) for p in partitions),
        len(cross_edges),
    )
    return partitions, cross_edges


# ── Cross-edge greedy post-processing ─────────────────────────────────────────


def _incident_cross_cut(
    node: int,
    cross_edges: list[tuple[int, int, float]],
    assignments: dict[int, int],
) -> float:
    """Sum of cross-edge weights cut by *node* with its current assignment."""
    total = 0.0
    for u, v, w in cross_edges:
        if u == node or v == node:
            if assignments.get(u, 0) != assignments.get(v, 0):
                total += w
    return total


def greedy_cross_cut(
    cross_edges: list[tuple[int, int, float]],
    intra_assignments: dict[int, int],
) -> dict[int, int]:
    """Greedily improve cross-partition cut given existing intra assignments.

    For each node touching a cross edge, flips its bit if that strictly
    increases the cut weight.  Iterates until no further improvement.

    Args:
        cross_edges: ``(u, v, weight)`` edges spanning partition boundaries.
        intra_assignments: QAOA-derived ``{node: bit}`` for all intra nodes.

    Returns:
        Updated ``{node: bit}`` assignment dict.
    """
    assignments = dict(intra_assignments)
    cross_nodes = {u for u, _, _ in cross_edges} | {v for _, v, _ in cross_edges}

    changed = True
    while changed:
        changed = False
        for node in cross_nodes:
            before = _incident_cross_cut(node, cross_edges, assignments)
            assignments[node] ^= 1
            after = _incident_cross_cut(node, cross_edges, assignments)
            if after > before:
                changed = True
            else:
                assignments[node] ^= 1  # revert

    return assignments


# ── Bitstring extraction ───────────────────────────────────────────────────────


def _best_bitstring(
    params: np.ndarray,
    cost_hamiltonian: "qml.Hamiltonian",
    n_qubits: int,
    p: int,
) -> str:
    """Return the most-probable computational-basis bitstring at *params*.

    Builds a new QNode returning ``qml.probs`` over all wires and picks
    the maximum-probability index.
    """
    import pennylane as qml
    from qaoa_circuit import qaoa_ansatz

    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def _prob_circuit(p_arr: np.ndarray) -> np.ndarray:
        qaoa_ansatz(p_arr, cost_hamiltonian, n_qubits, p)
        return qml.probs(wires=list(range(n_qubits)))

    probs = np.array(_prob_circuit(params))
    best_idx = int(np.argmax(probs))
    return format(best_idx, f"0{n_qubits}b")


# ── Main decomposed-QAOA entry point ─────────────────────────────────────────


def run_decomposed_qaoa(
    graph: nx.Graph,
    threshold: int = QUBIT_THRESHOLD_DEFAULT,
    p: int = 1,
    optimizer: str = "COBYLA",
    maxiter: int = 100,
    seed: int = 42,
) -> DecompositionResult:
    """Run QAOA on a large graph via Kernighan-Lin decomposition.

    Steps:
      1. Partition *graph* into subgraphs of ≤ *threshold* nodes.
      2. For each partition, relabel nodes to 0…k-1, build cost Hamiltonian,
         optimise with QAOA, and extract the best bitstring.
      3. Collect intra-partition assignments; run greedy cross-edge pass.
      4. Return combined cut and per-partition diagnostics.

    Args:
        graph: Large input graph (must have ``weight`` on every edge).
        threshold: Qubit threshold — max nodes per partition.
        p: QAOA depth.
        optimizer: Classical optimizer (COBYLA | L-BFGS-B | SPSA | Adam).
        maxiter: Max optimizer iterations per partition.
        seed: Random seed.

    Returns:
        :class:`DecompositionResult` with combined cut information.
    """
    import pennylane as qml
    from graph_generator import build_cost_hamiltonian
    from qaoa_circuit import build_qaoa_circuit
    from optimizer import optimize_qaoa

    np.random.seed(seed)

    partitions, cross_edges = decompose_graph(graph, threshold, seed)
    intra_assignments: dict[int, int] = {}
    intra_cut = 0.0

    for idx, part in enumerate(partitions):
        n_sub = len(part.nodes)
        log.info(
            "  → Partition %d/%d: %d nodes, %d intra-edges",
            idx + 1, len(partitions), n_sub, part.subgraph.number_of_edges(),
        )

        # Relabel to contiguous 0…k-1
        relabel = {v: i for i, v in enumerate(part.nodes)}
        inv_relabel = {i: v for v, i in relabel.items()}
        sub_r = nx.relabel_nodes(part.subgraph, relabel)

        cost_H = build_cost_hamiltonian(sub_r)
        dev = qml.device("default.qubit", wires=n_sub)
        qnode = build_qaoa_circuit(cost_H, n_sub, p, dev)

        # Coarse grid warm-start (same logic as benchmark_runner._grid_warmstart)
        gammas_ws = np.linspace(0.05, np.pi - 0.05, 5)
        betas_ws = np.linspace(0.05, np.pi - 0.05, 5)
        best_e, best_init = -np.inf, np.array([np.pi / 4, np.pi / 8])
        for gam_ws in gammas_ws:
            for bet_ws in betas_ws:
                p_arr = np.concatenate([np.full(p, gam_ws), np.full(p, bet_ws)])
                e_ws = float(qnode(p_arr))
                if e_ws > best_e:
                    best_e, best_init = e_ws, p_arr.copy()

        result = optimize_qaoa(
            qnode, p=p, method=optimizer,
            maxiter=maxiter, init_params=best_init,
        )
        part.qaoa_result = {
            "optimal_energy": float(result.optimal_energy),
            "optimal_params": [float(x) for x in result.optimal_params],
            "n_evals": int(result.n_evals),
            "n_nodes": n_sub,
        }
        # Record per-partition QAOA expectation (for diagnostics only)
        intra_cut += float(result.optimal_energy)

        bs = _best_bitstring(result.optimal_params, cost_H, n_sub, p)
        for local_i, bit in enumerate(bs):
            part.assignments[inv_relabel[local_i]] = int(bit)
        intra_assignments.update(part.assignments)

    final_assignments = greedy_cross_cut(cross_edges, intra_assignments)

    # Compute the TRUE integer cut from the combined bitstring assignment.
    # (Summing independent QAOA expectations + greedy cross would be apples-to-
    # oranges and could exceed the global C_max.)
    total_cut = float(sum(
        d.get("weight", 1.0)
        for u, v, d in graph.edges(data=True)
        if final_assignments.get(u, 0) != final_assignments.get(v, 0)
    ))
    # Recompute per-component cuts for diagnostics
    intra_cut = float(sum(
        d.get("weight", 1.0)
        for part in partitions
        for u, v, d in part.subgraph.edges(data=True)
        if final_assignments.get(u, 0) != final_assignments.get(v, 0)
    ))
    cross_cut = float(sum(
        w for u, v, w in cross_edges
        if final_assignments.get(u, 0) != final_assignments.get(v, 0)
    ))

    c_max_upper = float(sum(
        d.get("weight", 1.0) for _, _, d in graph.edges(data=True)
    ))
    approx_ratio = total_cut / c_max_upper if c_max_upper > 0 else 0.0

    log.info(
        "Decomposed QAOA: intra=%.3f  cross=%.3f  total=%.3f  "
        "ratio=%.4f (upper bound C_max=%.2f)",
        intra_cut, cross_cut, total_cut, approx_ratio, c_max_upper,
    )

    return DecompositionResult(
        partitions=partitions,
        cross_edges=cross_edges,
        assignments=final_assignments,
        total_cut=total_cut,
        intra_cut=intra_cut,
        cross_cut=cross_cut,
        c_max_upper=c_max_upper,
        approximation_ratio=approx_ratio,
        n_partitions=len(partitions),
        largest_partition=max(len(p.nodes) for p in partitions),
    )
