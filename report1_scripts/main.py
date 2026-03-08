"""QAOA Energy Landscape Explorer — CLI entry point and pipeline orchestration."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
import numpy as np
import pennylane as qml

from graph_decomposer import run_decomposed_qaoa
from graph_generator import create_graph, build_cost_hamiltonian
from qaoa_circuit import build_qaoa_circuit, qaoa_ansatz
from landscape_sampler import sample_landscape, LandscapeData
from optimizer import optimize_qaoa, OptimizationResult

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).parent / "output"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the QAOA explorer."""
    ap = argparse.ArgumentParser(description="QAOA Energy Landscape Explorer")
    ap.add_argument("--graph-type", default="cycle",
                    choices=["cycle", "complete", "erdos_renyi",
                             "random_regular", "grid_2d", "petersen"])
    ap.add_argument("--n-nodes", type=int, default=4)
    ap.add_argument("--p", type=int, default=1, help="QAOA depth")
    ap.add_argument("--optimizer", default="COBYLA",
                    choices=["COBYLA", "L-BFGS-B", "SPSA", "Adam"])
    ap.add_argument("--grid-size", type=int, default=50)
    ap.add_argument("--noise", action="store_true",
                    help="Enable noisy simulation")
    ap.add_argument("--zne", action="store_true",
                    help="Enable zero-noise extrapolation")
    ap.add_argument("--render", action="store_true",
                    help="Launch OpenGL 3-D visualisation")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument(
        "--qubit-threshold", type=int, default=12,
        help="Max nodes before Kernighan-Lin decomposition is used (default 12).",
    )
    return ap.parse_args()


def brute_force_maxcut(graph: nx.Graph) -> tuple[float, str]:
    """Exhaustive search for the maximum cut value.

    Feasible for n <= ~20 nodes.

    Returns:
        (max_cut_value, best_bitstring)
    """
    n = graph.number_of_nodes()
    nodes = list(graph.nodes())
    best_cut = 0.0
    best_bs = "0" * n
    for k in range(2 ** n):
        bs = format(k, f"0{n}b")
        cut = 0.0
        for u, v, d in graph.edges(data=True):
            ui = nodes.index(u)
            vi = nodes.index(v)
            if bs[ui] != bs[vi]:
                cut += d.get("weight", 1.0)
        if cut > best_cut:
            best_cut = cut
            best_bs = bs
    return best_cut, best_bs


def plot_graph(graph: nx.Graph, save_path: Path) -> None:
    """Draw the graph and save to file."""
    fig, ax = plt.subplots(figsize=(5, 5))
    pos = nx.spring_layout(graph, seed=42)
    nx.draw(graph, pos, ax=ax, with_labels=True, node_color="skyblue",
            edge_color="gray", node_size=500, font_size=12, width=2)
    edge_labels = {(u, v): f"{d['weight']:.1f}"
                   for u, v, d in graph.edges(data=True)}
    nx.draw_networkx_edge_labels(graph, pos, edge_labels, ax=ax)
    ax.set_title(f"Graph ({graph.number_of_nodes()} nodes, "
                 f"{graph.number_of_edges()} edges)")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved graph plot → %s", save_path)


def plot_landscape_2d(data: LandscapeData, save_path: Path) -> None:
    """Heatmap of the energy landscape."""
    fig, ax = plt.subplots(figsize=(7, 5))
    extent = [data.beta_values[0], data.beta_values[-1],
              data.gamma_values[0], data.gamma_values[-1]]
    im = ax.imshow(data.energies, origin="lower", aspect="auto",
                   extent=extent, cmap="RdYlBu_r")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\gamma$")
    ax.set_title(r"$\langle H_C \rangle(\gamma, \beta)$")
    fig.colorbar(im, ax=ax, label=r"$\langle H_C \rangle$")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved 2-D landscape → %s", save_path)


def plot_landscape_3d(data: LandscapeData, save_path: Path,
                      trajectory: list[tuple[np.ndarray, float]] | None = None,
                      ) -> None:
    """3-D surface plot of the energy landscape with optional trajectory."""
    G, B = np.meshgrid(data.gamma_values, data.beta_values, indexing="ij")
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(G, B, data.energies, cmap="RdYlBu_r",
                           alpha=0.85, edgecolor="none")
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$\beta$")
    ax.set_zlabel(r"$\langle H_C \rangle$")
    ax.set_title("QAOA Energy Landscape")
    fig.colorbar(surf, ax=ax, shrink=0.5, label=r"$\langle H_C \rangle$")

    if trajectory:
        gammas = [t[0][0] for t in trajectory]
        betas = [t[0][1] for t in trajectory]  # index p for p=1 → index 1
        energies = [t[1] for t in trajectory]
        ax.plot(gammas, betas, energies, "o-", color="black",
                markersize=3, linewidth=1.5, zorder=10,
                label="Optimizer trajectory")
        ax.scatter([gammas[-1]], [betas[-1]], [energies[-1]],
                   color="red", s=80, zorder=11, label="Optimum")
        ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved 3-D landscape → %s", save_path)


def plot_convergence(trajectory: list[tuple[np.ndarray, float]],
                     save_path: Path) -> None:
    """Plot optimization convergence curve."""
    energies = [t[1] for t in trajectory]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(energies, "o-", markersize=3)
    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"$\langle H_C \rangle$")
    ax.set_title("Optimizer Convergence")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved convergence plot → %s", save_path)


def plot_approximation_ratio(data: LandscapeData, save_path: Path) -> None:
    """Heatmap of approximation ratio r = E/C_max."""
    fig, ax = plt.subplots(figsize=(7, 5))
    extent = [data.beta_values[0], data.beta_values[-1],
              data.gamma_values[0], data.gamma_values[-1]]
    im = ax.imshow(data.approximation_ratios, origin="lower", aspect="auto",
                   extent=extent, cmap="viridis", vmin=0, vmax=1)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\gamma$")
    ax.set_title(r"Approximation ratio $r = \langle H_C \rangle / C_{\max}$")
    fig.colorbar(im, ax=ax, label="$r$")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved approximation-ratio heatmap → %s", save_path)


def plot_noisy_comparison(ideal_energy: float, noisy_energy: float,
                          mitigated_energy: float | None,
                          c_max: float, save_path: Path) -> None:
    """Bar chart comparing ideal, noisy, and mitigated energies."""
    labels = ["Ideal", "Noisy"]
    values = [ideal_energy, noisy_energy]
    colors = ["#2196F3", "#FF9800"]
    if mitigated_energy is not None:
        labels.append("ZNE-Mitigated")
        values.append(mitigated_energy)
        colors.append("#4CAF50")

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, values, color=colors, edgecolor="black", width=0.5)
    ax.axhline(c_max, color="red", linestyle="--", label=f"$C_{{\\max}} = {c_max:.2f}$")
    ax.set_ylabel(r"$\langle H_C \rangle$")
    ax.set_title("Energy Comparison: Ideal vs Noisy vs Mitigated")
    ax.legend()
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f"{val:.3f}", ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved noise comparison → %s", save_path)


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full QAOA exploration pipeline."""
    np.random.seed(args.seed)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Graph generation ───────────────────────────────────────────────
    log.info("=== Stage 1: Graph Generation ===")
    graph = create_graph(args.graph_type, args.n_nodes, seed=args.seed)
    cost_H = build_cost_hamiltonian(graph)
    c_max, best_bs = brute_force_maxcut(graph)
    log.info("Graph: %s with %d nodes, %d edges",
             args.graph_type, graph.number_of_nodes(), graph.number_of_edges())
    log.info("C_max = %.2f  (bitstring %s)", c_max, best_bs)
    plot_graph(graph, OUTPUT_DIR / "graph.pdf")

    # ── 2. Circuit construction (ideal simulator) ─────────────────────────
    log.info("=== Stage 2: Circuit Construction ===")
    n_qubits = graph.number_of_nodes()
    threshold = args.qubit_threshold

    if n_qubits > threshold:
        log.info(
            "Graph has %d nodes > threshold %d — using Kernighan-Lin decomposition.",
            n_qubits, threshold,
        )
        dr = run_decomposed_qaoa(
            graph, threshold=threshold, p=args.p,
            optimizer=args.optimizer, maxiter=args.maxiter, seed=args.seed,
        )
        results = {
            "graph_type": args.graph_type,
            "n_nodes": int(args.n_nodes),
            "n_edges": int(graph.number_of_edges()),
            "p": int(args.p),
            "c_max_upper": float(dr.c_max_upper),
            "optimizer": args.optimizer,
            "total_cut": float(dr.total_cut),
            "intra_cut": float(dr.intra_cut),
            "cross_cut": float(dr.cross_cut),
            "approximation_ratio": float(dr.approximation_ratio),
            "n_partitions": int(dr.n_partitions),
            "largest_partition": int(dr.largest_partition),
        }
        results_path = OUTPUT_DIR / "results.json"
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)
        log.info("Saved results → %s", results_path)
        log.info("=== Pipeline complete (decomposed) ===")
        return

    dev = qml.device("default.qubit", wires=n_qubits)
    qnode = build_qaoa_circuit(cost_H, n_qubits, args.p, dev)
    log.info("Device: default.qubit (%d wires), depth p=%d", n_qubits, args.p)

    # ── 3. Landscape sampling ─────────────────────────────────────────────
    log.info("=== Stage 3: Landscape Sampling ===")
    landscape = sample_landscape(
        qnode, p=args.p,
        gamma_range=(0, np.pi),
        beta_range=(0, np.pi / 2),
        n_gamma=args.grid_size,
        n_beta=args.grid_size,
        c_max=c_max,
    )
    log.info("Sampled %d × %d grid", args.grid_size, args.grid_size)
    plot_landscape_2d(landscape, OUTPUT_DIR / "landscape_2d.pdf")
    plot_approximation_ratio(landscape, OUTPUT_DIR / "approx_ratio.pdf")

    # ── 4. Optimization ─────────────────────────────────────────────────
    log.info("=== Stage 4: Optimization ===")
    opt_result = optimize_qaoa(
        qnode, p=args.p, method=args.optimizer, maxiter=args.maxiter,
    )
    log.info("Optimizer: %s  converged=%s  n_evals=%d",
             args.optimizer, opt_result.converged, opt_result.n_evals)
    log.info("Optimal energy: %.6f  (approx ratio: %.4f)",
             opt_result.optimal_energy, opt_result.optimal_energy / c_max)
    log.info("Optimal params: %s", opt_result.optimal_params)
    plot_convergence(opt_result.trajectory, OUTPUT_DIR / "convergence.pdf")

    # For p=1, overlay trajectory on landscape
    if args.p == 1:
        traj_for_plot = [
            (t[0], t[1]) for t in opt_result.trajectory
        ]
        plot_landscape_3d(landscape, OUTPUT_DIR / "landscape_3d.pdf",
                          trajectory=opt_result.trajectory)
    else:
        plot_landscape_3d(landscape, OUTPUT_DIR / "landscape_3d.pdf")

    ideal_energy = opt_result.optimal_energy

    # ── 5. Noisy simulation (optional) ────────────────────────────────────
    noisy_energy = None
    mitigated_energy = None
    if args.noise:
        log.info("=== Stage 5: Noisy Simulation ===")
        from noise_model import create_noisy_device, noisy_qaoa_qnode
        noise_params = {
            "depol_rate": 0.01,
            "t1": 50e-6,
            "t2": 70e-6,
            "gate_time_1q": 50e-9,
            "gate_time_2q": 300e-9,
        }
        noisy_qn = noisy_qaoa_qnode(cost_H, n_qubits, args.p, noise_params)
        gammas = opt_result.optimal_params[:args.p]
        betas = opt_result.optimal_params[args.p:]
        noisy_energy = float(noisy_qn(gammas, betas))
        log.info("Noisy energy: %.6f", noisy_energy)

        # ── 5b. ZNE (optional) ─────────────────────────────────────────────
        if args.zne:
            log.info("=== Stage 5b: Zero-Noise Extrapolation ===")
            from error_mitigation import zne_extrapolate
            # Build a qnode that takes a single params array for ZNE
            noisy_dev = qml.device("default.mixed", wires=n_qubits)
            noisy_qnode_for_zne = build_qaoa_circuit(
                cost_H, n_qubits, args.p, noisy_dev,
            )
            mitigated_energy, zne_details = zne_extrapolate(
                noisy_qnode_for_zne,
                opt_result.optimal_params,
                scale_factors=[1.0, 2.0, 3.0],
                extrapolation="richardson",
            )
            log.info("ZNE mitigated energy: %.6f", mitigated_energy)
            log.info("ZNE details: %s", {
                k: (v.tolist() if isinstance(v, np.ndarray) else v)
                for k, v in zne_details.items()
            })

        plot_noisy_comparison(
            ideal_energy, noisy_energy, mitigated_energy, c_max,
            OUTPUT_DIR / "noise_comparison.pdf",
        )

    # ── 6. Save results summary ───────────────────────────────────────────
    results = {
        "graph_type": args.graph_type,
        "n_nodes": int(args.n_nodes),
        "n_edges": int(graph.number_of_edges()),
        "p": int(args.p),
        "c_max": float(c_max),
        "best_bitstring": best_bs,
        "optimizer": args.optimizer,
        "optimal_energy": float(opt_result.optimal_energy),
        "approximation_ratio": float(opt_result.optimal_energy / c_max),
        "optimal_params": [float(x) for x in opt_result.optimal_params],
        "n_evals": int(opt_result.n_evals),
        "converged": bool(opt_result.converged),
    }
    if noisy_energy is not None:
        results["noisy_energy"] = float(noisy_energy)
    if mitigated_energy is not None:
        results["mitigated_energy"] = float(mitigated_energy)

    results_path = OUTPUT_DIR / "results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    log.info("Saved results → %s", results_path)

    # ── 7. OpenGL renderer (optional) ─────────────────────────────────────
    if args.render:
        log.info("=== Stage 7: OpenGL Rendering ===")
        try:
            _launch_opengl(landscape, opt_result)
        except Exception as exc:
            log.warning("OpenGL render failed: %s", exc)

    log.info("=== Pipeline complete ===")


def _launch_opengl(
    landscape: LandscapeData,
    opt_result: OptimizationResult,
) -> None:
    """Launch the interactive OpenGL 3-D surface viewer."""
    from renderer import (
        init_opengl_context, compile_shader_program,
        create_surface_vbo, draw_optimizer_trajectory, ArcballCamera,
    )
    import glfw
    from OpenGL.GL import (
        glClear, glClearColor, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT,
        glUseProgram, glGetUniformLocation, glUniformMatrix4fv, glUniform3f,
        glBindVertexArray, glDrawElements, GL_TRIANGLES, GL_UNSIGNED_INT,
        GL_FALSE,
    )
    import ctypes

    shader_dir = Path(__file__).parent / "shaders"
    window = init_opengl_context()
    surf_prog = compile_shader_program(
        (shader_dir / "surface.vert").read_text(),
        (shader_dir / "surface.frag").read_text(),
    )
    traj_prog = compile_shader_program(
        (shader_dir / "trajectory.vert").read_text(),
        (shader_dir / "trajectory.frag").read_text(),
    )

    vao, vbo, ebo = create_surface_vbo(
        landscape.gamma_values, landscape.beta_values, landscape.energies,
    )
    n_indices = 6 * (len(landscape.gamma_values) - 1) * (len(landscape.beta_values) - 1)

    cam = ArcballCamera(
        position=np.array([np.pi / 2, np.pi / 4, 5.0]),
        target=np.array([np.pi / 2, np.pi / 4, landscape.energies.mean()]),
    )

    # Mouse state
    mouse_state = {"last_x": 0.0, "last_y": 0.0, "pressed": False}

    def cursor_cb(win, xpos, ypos):
        if mouse_state["pressed"]:
            cam.on_mouse_drag(xpos - mouse_state["last_x"],
                              mouse_state["last_y"] - ypos)
        mouse_state["last_x"] = xpos
        mouse_state["last_y"] = ypos

    def button_cb(win, button, action, mods):
        if button == glfw.MOUSE_BUTTON_LEFT:
            mouse_state["pressed"] = (action == glfw.PRESS)

    def scroll_cb(win, xoff, yoff):
        cam.on_scroll(yoff)

    glfw.set_cursor_pos_callback(window, cursor_cb)
    glfw.set_mouse_button_callback(window, button_cb)
    glfw.set_scroll_callback(window, scroll_cb)

    trajectory_3d = [
        (t[0][0], t[0][1] if len(t[0]) > 1 else 0.0, t[1])
        for t in opt_result.trajectory
    ]

    glClearColor(0.1, 0.1, 0.15, 1.0)

    while not glfw.window_should_close(window):
        glfw.poll_events()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        w, h = glfw.get_framebuffer_size(window)
        aspect = w / max(h, 1)
        view = cam.get_view_matrix()
        proj = cam.get_projection_matrix(aspect)
        model = np.eye(4, dtype=np.float32)

        # Draw surface
        glUseProgram(surf_prog)
        glUniformMatrix4fv(glGetUniformLocation(surf_prog, "uModel"),
                           1, GL_FALSE, model)
        glUniformMatrix4fv(glGetUniformLocation(surf_prog, "uView"),
                           1, GL_FALSE, view)
        glUniformMatrix4fv(glGetUniformLocation(surf_prog, "uProjection"),
                           1, GL_FALSE, proj)
        glUniform3f(glGetUniformLocation(surf_prog, "uLightPos"),
                    5.0, 5.0, 5.0)
        glUniform3f(glGetUniformLocation(surf_prog, "uViewPos"),
                    *cam.position)

        glBindVertexArray(vao)
        glDrawElements(GL_TRIANGLES, n_indices, GL_UNSIGNED_INT,
                       ctypes.c_void_p(0))
        glBindVertexArray(0)

        # Draw trajectory
        mvp = proj @ view @ model
        glUseProgram(traj_prog)
        glUniformMatrix4fv(glGetUniformLocation(traj_prog, "uMVP"),
                           1, GL_FALSE, mvp)
        draw_optimizer_trajectory(trajectory_3d, traj_prog)

        glfw.swap_buffers(window)

    glfw.terminate()


def main() -> None:
    """Entry point: parse arguments and run the pipeline."""
    args = parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
