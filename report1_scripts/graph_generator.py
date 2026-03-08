"""Graph generation and Ising Hamiltonian construction for MaxCut QAOA."""

from __future__ import annotations

import numpy as np
import networkx as nx
import pennylane as qml


def create_graph(
    graph_type: str,
    n_nodes: int,
    seed: int | None = None,
    **kwargs,
) -> nx.Graph:
    """Create a weighted undirected NetworkX graph of the specified family.

    Supported graph types:
      - 'cycle': Cycle graph on *n_nodes* vertices.
      - 'complete': Complete graph K_n.
      - 'erdos_renyi': Erdős–Rényi G(n, p) random graph.
      - 'random_regular': Random *d*-regular graph.
      - 'grid_2d': 2-D grid graph (rows × cols ≈ n_nodes).
      - 'petersen': Petersen graph (ignores *n_nodes*).
      - 'custom': Build from an explicit edge list.

    After construction, ensure every edge carries a ``"weight"`` attribute
    (default 1.0) so downstream code can assume weighted edges.

    Args:
        graph_type: Graph family name (see list above).
        n_nodes: Number of vertices.
        seed: Random seed for stochastic families (Erdős–Rényi, regular).
        **kwargs: Family-specific parameters:
            - ``p`` (float): Edge probability for Erdős–Rényi (default 0.5).
            - ``d`` (int): Degree for random regular graphs (default 3).
            - ``rows`` (int): Row count for 2-D grid.
            - ``edges`` (list): Edge list for 'custom' — each element is
              ``(i, j)`` or ``(i, j, weight)``.

    Returns:
        Weighted undirected ``nx.Graph``.

    Raises:
        ValueError: If *graph_type* is not recognised.
    """
    if graph_type == "cycle":
        G = nx.cycle_graph(n_nodes)
    elif graph_type == 'complete':
        G = nx.complete_graph(n_nodes)
    elif graph_type == 'erdos_renyi':
        G = nx.erdos_renyi_graph(n_nodes, kwargs.get('p', 0.5))
    elif graph_type == 'random_regular':
        G = nx.random_regular_graph(kwargs.get('d', 3), n_nodes, seed=seed)
    elif graph_type == 'grid_2d':
        G = nx.grid_2d_graph(kwargs.get('rows', int(n_nodes**0.5)), kwargs.get('cols', int(n_nodes**0.5)))
    elif graph_type == 'petersen':
        G = nx.petersen_graph()
    elif graph_type == 'custom':
        G = nx.from_edgelist(kwargs['edges'])
    else:
        raise ValueError(f"Unknown graph type: {graph_type}")

    # Ensure all edges have a 'weight' attribute
    for u, v in G.edges():
        if 'weight' not in G[u][v]:
            G[u][v]['weight'] = 1.0

    return G


def adjacency_to_ising(
    graph: nx.Graph,
) -> tuple[list[float], list[tuple[int, int]], float]:
    """Convert a weighted graph to Ising ZZ-coefficients for MaxCut.

    The MaxCut cost Hamiltonian is:

        H_C = Σ_{(i,j)} w_ij / 2 · (I − Z_i Z_j)

    This function extracts the ZZ coefficients (−w/2) and the constant
    energy offset (Σ w/2) so that the Hamiltonian can be reconstructed as:

        H_C = offset · I  +  Σ_k coeffs[k] · Z_{obs[k][0]} Z_{obs[k][1]}

    Args:
        graph: Weighted undirected graph (every edge must have ``"weight"``).

    Returns:
        A 3-tuple ``(coeffs, obs, offset)`` where:
          - *coeffs*: list of float ZZ coefficients (−w/2 per edge).
          - *obs*: list of ``(u, v)`` qubit-index pairs.
          - *offset*: constant energy shift (float).
    """
    coeffs = []
    obs = []
    offset = 0

    for u, v in graph.edges():
        coeffs.append(-graph[u][v]['weight'] / 2)
        obs.append((u, v))
        offset += graph[u][v]['weight'] / 2

    return coeffs, obs, offset


def build_cost_hamiltonian(graph: nx.Graph) -> qml.Hamiltonian:
    """Build the PennyLane Hamiltonian for the MaxCut cost function.

    Constructs:

        H_C = Σ_{(u,v)} (−w/2) Z_u Z_v  +  (Σ w/2) I

    The Identity term preserves the absolute energy scale so that
    ``⟨H_C⟩`` directly equals the expected cut value.

    Args:
        graph: Weighted undirected graph with ``"weight"`` on every edge.

    Returns:
        ``qml.Hamiltonian`` ready for use as a QAOA cost observable.
    """
    coeffs, obs, offset = adjacency_to_ising(graph)
    observables = [qml.PauliZ(u) @ qml.PauliZ(v) for (u, v) in obs]
    
    # Append the identity offset term explicitly
    coeffs.append(offset)
    observables.append(qml.Identity(wires=list(graph.nodes())))
    
    return qml.Hamiltonian(coeffs, observables)
