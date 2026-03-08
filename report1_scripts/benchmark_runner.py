"""Comparative benchmark: multiple graph families × four classical optimizers.

Runs QAOA (p=1) on a suite of graphs for each optimizer (COBYLA, L-BFGS-B,
SPSA, Adam) and saves comparison plots to output/benchmark/.

Graphs exceeding QUBIT_THRESHOLD nodes are handled transparently via the
Kernighan-Lin decomposition in graph_decomposer.py.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import networkx as nx
import numpy as np

from graph_decomposer import (
    QUBIT_THRESHOLD_DEFAULT,
    decompose_graph,
    run_decomposed_qaoa,
)
from graph_generator import create_graph, build_cost_hamiltonian
from optimizer import optimize_qaoa

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).parent / "output" / "benchmark"
QUBIT_THRESHOLD = QUBIT_THRESHOLD_DEFAULT  # 12 nodes

# ── Benchmark suite ───────────────────────────────────────────────────────────

GRAPH_CASES: list[dict] = [
    {"graph_type": "cycle",          "n_nodes": 4,  "label": "C₄"},
    {"graph_type": "cycle",          "n_nodes": 6,  "label": "C₆"},
    {"graph_type": "cycle",          "n_nodes": 8,  "label": "C₈"},
    {"graph_type": "complete",       "n_nodes": 4,  "label": "K₄"},
    {"graph_type": "complete",       "n_nodes": 5,  "label": "K₅"},
    {"graph_type": "erdos_renyi",    "n_nodes": 6,  "label": "ER₆"},
    {"graph_type": "erdos_renyi",    "n_nodes": 8,  "label": "ER₈"},
    {"graph_type": "random_regular", "n_nodes": 6,  "label": "Reg₆"},
    {"graph_type": "random_regular", "n_nodes": 8,  "label": "Reg₈"},
    # Large graphs — decomposed automatically
    {"graph_type": "erdos_renyi",    "n_nodes": 16, "label": "ER₁₆*"},
    {"graph_type": "cycle",          "n_nodes": 20, "label": "C₂₀*"},
]

OPTIMIZERS: list[str] = ["COBYLA", "L-BFGS-B", "SPSA", "Adam"]

COLORS: dict[str, str] = {
    "COBYLA":   "#2196F3",
    "L-BFGS-B": "#FF9800",
    "SPSA":     "#4CAF50",
    "Adam":     "#E91E63",
}

MAXITER_MAP: dict[str, int] = {
    "COBYLA":   150,
    "L-BFGS-B": 80,
    "SPSA":     200,
    "Adam":     100,
}

SEED = 42


# ── Result container ──────────────────────────────────────────────────────────


@dataclass
class BenchmarkRecord:
    """One (graph, optimizer) measurement."""

    graph_label: str
    graph_type: str
    n_nodes: int
    n_edges: int
    optimizer: str
    p: int = 1
    decomposed: bool = False
    optimal_energy: float = float("nan")
    c_max: float = float("nan")
    approximation_ratio: float = float("nan")
    n_evals: int = 0
    wall_time: float = float("nan")
    converged: bool = False
    trajectory_energies: list[float] = field(default_factory=list)


# ── Brute-force MaxCut (small graphs only) ───────────────────────────────────


def _brute_force_maxcut(graph: nx.Graph) -> float:
    """Exhaustive MaxCut for n ≤ 20."""
    n = graph.number_of_nodes()
    if n > 20:
        # Upper bound: sum of all edge weights
        return float(sum(d.get("weight", 1.0) for _, _, d in graph.edges(data=True)))
    nodes = list(graph.nodes())
    best = 0.0
    for k in range(2 ** n):
        bs = format(k, f"0{n}b")
        cut = sum(
            d.get("weight", 1.0)
            for u, v, d in graph.edges(data=True)
            if bs[nodes.index(u)] != bs[nodes.index(v)]
        )
        if cut > best:
            best = cut
    return best


# ── Warm-start helper ─────────────────────────────────────────────────────────


def _grid_warmstart(
    qnode: "qml.QNode",
    p: int = 1,
    n_gamma: int = 6,
    n_beta: int = 6,
) -> np.ndarray:
    """Return the best parameter vector from a coarse grid scan.

    Evaluates ``n_gamma × n_beta`` evenly-spaced (γ, β) pairs over
    [0, π] × [0, π/2] and returns the one with the highest energy.
    For p > 1 the scan covers only the first (γ₀, β₀) dimension;
    remaining parameters are set to π/4 and π/8 respectively.

    Args:
        qnode: QAOA qnode.
        p: QAOA depth.
        n_gamma: Grid resolution along γ axis.
        n_beta: Grid resolution along β axis.

    Returns:
        Best ``params`` array found (shape ``(2·p,)``).
    """
    gammas = np.linspace(0.05, np.pi - 0.05, n_gamma)
    betas = np.linspace(0.05, np.pi - 0.05, n_beta)  # full [0,π] to capture all optima

    best_energy = -np.inf
    best_params = np.concatenate([
        np.full(p, np.pi / 4),
        np.full(p, np.pi / 8),
    ])

    for g in gammas:
        for b in betas:
            params = np.concatenate([
                np.full(p, g),
                np.full(p, b),
            ])
            e = float(qnode(params))
            if e > best_energy:
                best_energy = e
                best_params = params.copy()

    return best_params


# ── Single benchmark run ──────────────────────────────────────────────────────


def run_single(
    graph: nx.Graph,
    graph_label: str,
    graph_type: str,
    optimizer: str,
    p: int = 1,
    seed: int = SEED,
) -> BenchmarkRecord:
    """Run one (graph, optimizer) combination and return a BenchmarkRecord."""
    import pennylane as qml
    from qaoa_circuit import build_qaoa_circuit

    n = graph.number_of_nodes()
    c_max = _brute_force_maxcut(graph)
    needs_decomp = n > QUBIT_THRESHOLD
    maxiter = MAXITER_MAP[optimizer]

    rec = BenchmarkRecord(
        graph_label=graph_label,
        graph_type=graph_type,
        n_nodes=n,
        n_edges=graph.number_of_edges(),
        optimizer=optimizer,
        p=p,
        decomposed=needs_decomp,
        c_max=c_max,
    )

    t0 = time.perf_counter()

    if needs_decomp:
        # ── Decomposed path ────────────────────────────────────────────────
        log.info("  [DECOMPOSED] %s | %s  (n=%d > threshold=%d)",
                 graph_label, optimizer, n, QUBIT_THRESHOLD)
        dr = run_decomposed_qaoa(
            graph, threshold=QUBIT_THRESHOLD, p=p,
            optimizer=optimizer, maxiter=maxiter, seed=seed,
        )
        rec.optimal_energy = dr.total_cut
        rec.approximation_ratio = (
            dr.total_cut / c_max if c_max > 0 else dr.approximation_ratio
        )
        rec.n_evals = sum(
            pt.qaoa_result["n_evals"]
            for pt in dr.partitions
            if pt.qaoa_result is not None
        )
        rec.converged = True
        # Build a synthetic trajectory from per-partition energies
        rec.trajectory_energies = [
            sum(
                pt.qaoa_result["optimal_energy"]
                for pt in dr.partitions
                if pt.qaoa_result is not None
            )
        ]
    else:
        # ── Direct QAOA path ──────────────────────────────────────────────
        log.info("  %s | %s  (n=%d)", graph_label, optimizer, n)
        cost_H = build_cost_hamiltonian(graph)
        dev = qml.device("default.qubit", wires=n)
        qnode = build_qaoa_circuit(cost_H, n, p, dev)

        np.random.seed(seed)
        # Warm-start: coarse 6×6 grid scan (36 evaluations) to find a good
        # initial point, avoiding barren-plateau starts for gradient-based opts.
        init_params = _grid_warmstart(qnode, p=p, n_gamma=6, n_beta=6)
        result = optimize_qaoa(
            qnode, p=p, method=optimizer,
            maxiter=maxiter, init_params=init_params,
        )

        rec.optimal_energy = result.optimal_energy
        rec.approximation_ratio = (
            result.optimal_energy / c_max if c_max > 0 else 0.0
        )
        rec.n_evals = result.n_evals
        rec.converged = result.converged
        rec.trajectory_energies = [float(e) for _, e in result.trajectory]

    rec.wall_time = time.perf_counter() - t0
    log.info(
        "    → E=%.4f  r=%.4f  evals=%d  t=%.1fs",
        rec.optimal_energy, rec.approximation_ratio, rec.n_evals, rec.wall_time,
    )
    return rec


# ── Full benchmark sweep ──────────────────────────────────────────────────────


def run_all(seed: int = SEED) -> list[BenchmarkRecord]:
    """Run the full benchmark matrix and return all records."""
    records: list[BenchmarkRecord] = []
    total = len(GRAPH_CASES) * len(OPTIMIZERS)
    done = 0
    for case in GRAPH_CASES:
        graph = create_graph(case["graph_type"], case["n_nodes"], seed=seed)
        for opt in OPTIMIZERS:
            done += 1
            log.info("[%d/%d] %s — %s", done, total, case["label"], opt)
            rec = run_single(
                graph,
                graph_label=case["label"],
                graph_type=case["graph_type"],
                optimizer=opt,
                seed=seed,
            )
            records.append(rec)
    return records


# ── Plotting ──────────────────────────────────────────────────────────────────


def _by(records: list[BenchmarkRecord], **filters) -> list[BenchmarkRecord]:
    """Filter records by keyword equality."""
    out = records
    for k, v in filters.items():
        out = [r for r in out if getattr(r, k) == v]
    return out


def plot_convergence_comparison(
    records: list[BenchmarkRecord],
    representative_labels: list[str],
    save_path: Path,
) -> None:
    """Multi-panel convergence curves: one panel per graph, 4 optimizer lines."""
    n_panels = len(representative_labels)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 4.5),
                              sharey=False)
    if n_panels == 1:
        axes = [axes]

    for ax, label in zip(axes, representative_labels):
        for opt in OPTIMIZERS:
            matches = _by(records, graph_label=label, optimizer=opt)
            if not matches:
                continue
            rec = matches[0]
            traj = rec.trajectory_energies
            if not traj:
                ax.scatter([rec.n_evals], [rec.optimal_energy],
                           color=COLORS[opt], label=opt, zorder=5)
                continue
            ax.plot(range(1, len(traj) + 1), traj,
                    color=COLORS[opt], label=opt, linewidth=1.8, alpha=0.9)

        # Mark best possible
        c_max_vals = [r.c_max for r in _by(records, graph_label=label)]
        if c_max_vals:
            ax.axhline(c_max_vals[0], color="red", linestyle="--",
                       linewidth=1.2, label=f"$C_{{\\max}}={c_max_vals[0]:.1f}$")
        ax.set_title(label)
        ax.set_xlabel("Iteration")
        ax.set_ylabel(r"$\langle H_C \rangle$")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.25)

    fig.suptitle("Optimizer Convergence Comparison (p = 1)", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved convergence comparison → %s", save_path)


def plot_approx_ratio_heatmap(
    records: list[BenchmarkRecord],
    save_path: Path,
) -> None:
    """Heatmap of approximation ratios: graph labels (rows) × optimizers (cols).

    Decomposed graphs are marked with '*' in the row label.
    """
    # Use all graph labels that appear in records
    graph_labels = list(dict.fromkeys(r.graph_label for r in records))
    data = np.full((len(graph_labels), len(OPTIMIZERS)), float("nan"))

    for i, gl in enumerate(graph_labels):
        for j, opt in enumerate(OPTIMIZERS):
            matches = _by(records, graph_label=gl, optimizer=opt)
            if matches:
                data[i, j] = matches[0].approximation_ratio

    fig, ax = plt.subplots(figsize=(len(OPTIMIZERS) * 1.6 + 1.5,
                                     len(graph_labels) * 0.55 + 1.5))
    im = ax.imshow(data, aspect="auto", cmap="RdYlGn",
                   vmin=max(0.0, np.nanmin(data) - 0.05), vmax=1.0)
    fig.colorbar(im, ax=ax, label="Approximation ratio $r$", shrink=0.8)

    ax.set_xticks(range(len(OPTIMIZERS)))
    ax.set_xticklabels(OPTIMIZERS, rotation=30, ha="right")
    ax.set_yticks(range(len(graph_labels)))
    ax.set_yticklabels(graph_labels)

    for i in range(len(graph_labels)):
        for j in range(len(OPTIMIZERS)):
            if not np.isnan(data[i, j]):
                ax.text(j, i, f"{data[i, j]:.3f}", ha="center", va="center",
                        fontsize=8, color="black" if data[i, j] > 0.5 else "white")

    ax.set_title("Approximation Ratio  $r = E_{\\mathrm{opt}} / C_{\\max}$\n"
                 "(* = Kernighan-Lin decomposed)", fontsize=11)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved approximation-ratio heatmap → %s", save_path)


def plot_nevals(
    records: list[BenchmarkRecord],
    save_path: Path,
) -> None:
    """Grouped bar chart: circuit evaluations per (graph, optimizer)."""
    graph_labels = list(dict.fromkeys(r.graph_label for r in records))
    n_graphs = len(graph_labels)
    n_opts = len(OPTIMIZERS)
    x = np.arange(n_graphs)
    width = 0.8 / n_opts

    fig, ax = plt.subplots(figsize=(max(10, n_graphs * 1.2), 5))
    for k, opt in enumerate(OPTIMIZERS):
        heights = []
        for gl in graph_labels:
            m = _by(records, graph_label=gl, optimizer=opt)
            heights.append(m[0].n_evals if m else 0)
        ax.bar(x + k * width - 0.4 + width / 2, heights, width,
               label=opt, color=COLORS[opt], edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(graph_labels, rotation=35, ha="right")
    ax.set_ylabel("Circuit evaluations")
    ax.set_title("Number of Circuit Evaluations per Optimizer\n"
                 "(* = reported per-partition total for decomposed graphs)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved n-evals chart → %s", save_path)


def plot_wall_time(
    records: list[BenchmarkRecord],
    save_path: Path,
) -> None:
    """Wall-time comparison bar chart."""
    graph_labels = list(dict.fromkeys(r.graph_label for r in records))
    n_graphs = len(graph_labels)
    n_opts = len(OPTIMIZERS)
    x = np.arange(n_graphs)
    width = 0.8 / n_opts

    fig, ax = plt.subplots(figsize=(max(10, n_graphs * 1.2), 5))
    for k, opt in enumerate(OPTIMIZERS):
        heights = []
        for gl in graph_labels:
            m = _by(records, graph_label=gl, optimizer=opt)
            heights.append(m[0].wall_time if m else 0.0)
        ax.bar(x + k * width - 0.4 + width / 2, heights, width,
               label=opt, color=COLORS[opt], edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(graph_labels, rotation=35, ha="right")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Wall-Clock Time per Optimizer")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved wall-time chart → %s", save_path)


def plot_decomposition_graph(
    graph: nx.Graph,
    graph_label: str,
    save_path: Path,
    seed: int = SEED,
) -> None:
    """Draw a large graph with Kernighan-Lin partition colouring."""
    from graph_decomposer import decompose_graph

    partitions, cross_edges = decompose_graph(graph, QUBIT_THRESHOLD, seed)

    # Assign a colour to each partition
    cmap = plt.get_cmap("tab10")
    node_colors: dict[int, str] = {}
    for idx, part in enumerate(partitions):
        color = cmap(idx % 10)
        hex_color = "#{:02x}{:02x}{:02x}".format(
            int(color[0] * 255), int(color[1] * 255), int(color[2] * 255),
        )
        for node in part.nodes:
            node_colors[node] = hex_color

    fig, ax = plt.subplots(figsize=(8, 7))
    pos = nx.spring_layout(graph, seed=seed, k=1.5)

    # Draw intra-partition edges
    intra_edge_set = {
        (u, v) for part in partitions
        for u, v in part.subgraph.edges()
    }
    cross_edge_set = [(u, v) for u, v, _ in cross_edges]

    nx.draw_networkx_edges(graph, pos, edgelist=list(intra_edge_set),
                           edge_color="gray", width=1.5, ax=ax, alpha=0.7)
    nx.draw_networkx_edges(graph, pos, edgelist=cross_edge_set,
                           edge_color="red", width=2.0, ax=ax, alpha=0.9,
                           style="dashed")
    nx.draw_networkx_nodes(graph, pos, node_size=400,
                           node_color=[node_colors[n] for n in graph.nodes()],
                           ax=ax)
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=9)

    # Legend
    for idx, part in enumerate(partitions):
        c = cmap(idx % 10)
        ax.scatter([], [], color=c,
                   label=f"Partition {idx + 1}  ({len(part.nodes)} nodes)")
    ax.scatter([], [], marker="None", label="— intra-edge  -- cross-edge")
    ax.legend(fontsize=8, loc="upper right")

    title = (f"{graph_label}:  {graph.number_of_nodes()} nodes  "
             f"→  {len(partitions)} partitions  "
             f"({len(cross_edges)} cross-edges)")
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved decomposition graph → %s", save_path)


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    """Run the full benchmark suite and save all plots."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    log.info("=== Starting benchmark: %d graphs × %d optimizers ===",
             len(GRAPH_CASES), len(OPTIMIZERS))
    t_start = time.perf_counter()

    records = run_all(seed=SEED)

    elapsed = time.perf_counter() - t_start
    log.info("=== Benchmark complete in %.1f s ===", elapsed)

    # ── Save raw results ───────────────────────────────────────────────────
    summary = [
        {
            "graph_label": r.graph_label,
            "graph_type": r.graph_type,
            "n_nodes": r.n_nodes,
            "n_edges": r.n_edges,
            "optimizer": r.optimizer,
            "decomposed": bool(r.decomposed),
            "optimal_energy": float(r.optimal_energy),
            "c_max": float(r.c_max),
            "approximation_ratio": float(r.approximation_ratio),
            "n_evals": int(r.n_evals),
            "wall_time": float(r.wall_time),
            "converged": bool(r.converged),
        }
        for r in records
    ]
    with open(OUTPUT_DIR / "benchmark_results.json", "w") as f:
        json.dump(summary, f, indent=2)
    log.info("Saved raw results → %s", OUTPUT_DIR / "benchmark_results.json")

    # ── Plots ──────────────────────────────────────────────────────────────

    # 1. Convergence comparison (2 representative graphs)
    rep_labels = [c["label"] for c in GRAPH_CASES
                  if not c["n_nodes"] > QUBIT_THRESHOLD]
    # Pick two that are non-trivially interesting
    for candidate in ["C₆", "ER₈", "ER₆", "C₄"]:
        if candidate in rep_labels and rep_labels.index(candidate) < len(rep_labels):
            first_choice = candidate
            break
    else:
        first_choice = rep_labels[0]

    for candidate in ["Reg₈", "K₅", "ER₈"]:
        if candidate in rep_labels and candidate != first_choice:
            second_choice = candidate
            break
    else:
        second_choice = rep_labels[-1] if rep_labels[-1] != first_choice else rep_labels[0]

    plot_convergence_comparison(
        records, [first_choice, second_choice],
        OUTPUT_DIR / "convergence_comparison.pdf",
    )

    # 2. Approx-ratio heatmap
    plot_approx_ratio_heatmap(records, OUTPUT_DIR / "approx_ratio_heatmap.pdf")

    # 3. N_evals bar chart
    plot_nevals(records, OUTPUT_DIR / "nevals_comparison.pdf")

    # 4. Wall time
    plot_wall_time(records, OUTPUT_DIR / "wall_time.pdf")

    # 5. Decomposition graph(s)
    large_cases = [c for c in GRAPH_CASES if c["n_nodes"] > QUBIT_THRESHOLD]
    for case in large_cases[:2]:
        g = create_graph(case["graph_type"], case["n_nodes"], seed=SEED)
        fname = f"decomp_{case['label'].replace('*', '').replace('₁', '1').replace('₂', '2').replace('₀', '0')}.pdf"
        plot_decomposition_graph(g, case["label"],
                                 OUTPUT_DIR / fname, seed=SEED)

    log.info("All benchmark plots saved to %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
