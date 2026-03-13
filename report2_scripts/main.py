from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the VQE PES Explorer.

    Arguments:
      --molecule: "h2", "lih", or "h2o" (default "h2")
      --basis: Basis set name (default "sto-3g")
      --bond-range: "start,stop,n_points" for PES scan (default "0.3,3.0,20")
      --ansatz: "uccsd" or "hardware_efficient" (default "uccsd")
      --optimizer: "Adam", "GradientDescent", "COBYLA", or "SPSA" (default "Adam")
      --maxiter: Max VQE iterations, int (default 200)
      --hardware: store_true flag to enable hardware execution
      --backend: Hardware backend name (default "ibm_brisbane")
      --render: store_true flag to enable 3D orbital visualisation
      --benchmark: store_true flag to enable multi-method benchmark
      --orbital: MO index to visualise, int (default 0)
      --output: Output directory for results (default "output")

    Returns:
        Parsed argparse.Namespace.
    """
    parser = argparse.ArgumentParser(
        description="VQE Potential Energy Surface Explorer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Molecule specification
    parser.add_argument("--molecule", type=str, default="h2",
                        choices=["h2", "lih", "h2o"],
                        help="Molecule to simulate")
    parser.add_argument("--basis", type=str, default="sto-3g",
                        help="Basis set name")
    parser.add_argument("--bond-range", type=str, default="0.3,3.0,20",
                        help="Bond length range: start,stop,n_points")
    
    # VQE configuration
    parser.add_argument("--ansatz", type=str, default="uccsd",
                        choices=["uccsd", "hardware_efficient"],
                        help="Ansatz type")
    parser.add_argument("--optimizer", type=str, default="Adam",
                        choices=["Adam", "GradientDescent", "COBYLA", "SPSA"],
                        help="Optimizer for VQE")
    parser.add_argument("--maxiter", type=int, default=200,
                        help="Maximum VQE iterations")
    parser.add_argument("--learning-rate", type=float, default=0.1,
                        help="Learning rate for gradient-based optimizers")
    
    # Hardware execution
    parser.add_argument("--hardware", action="store_true",
                        help="Enable hardware execution")
    parser.add_argument("--backend", type=str, default="ibm_brisbane",
                        help="Hardware backend name")
    
    # Visualization
    parser.add_argument("--render", action="store_true",
                        help="Enable 3D orbital visualisation")
    parser.add_argument("--orbital", type=int, default=0,
                        help="MO index to visualise")
    
    # Benchmarking
    parser.add_argument("--benchmark", action="store_true",
                        help="Enable multi-method benchmark (HF, CISD, CCSD(T), FCI)")
    
    # Output
    parser.add_argument("--output", type=str, default="output",
                        help="Output directory for results")
    
    return parser.parse_args()


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full VQE pipeline based on parsed arguments.

    Steps:
      1. Parse args.bond_range string ("start,stop,n_points") into
         np.linspace(start, stop, n_points).
      2. Select molecule factory from molecule.create_h2 / create_lih / create_h2o
         using a dict mapping args.molecule -> factory function.
      3. Run PES scan via pes_scanner.scan_pes with selected factory, bond_lengths,
         basis, ansatz_type, optimizer, maxiter, compute_classical=args.benchmark.
      4. If args.benchmark: run benchmarker.run_benchmark(pes_data),
         print benchmarker.print_benchmark_table(bench).
      5. If args.render:
         a. Build equilibrium molecule via factory() (default bond length).
         b. Run classical_hf.run_hartree_fock on it.
         c. Convert atom coords to Bohr: coords * 1.8897259886.
         d. Evaluate MO on grid via orbital_grid.evaluate_orbital_on_grid.
         e. Extract dual isosurface via isosurface.extract_dual_isosurface
            (isovalue=0.05).
         f. Init OpenGL via renderer.init_opengl_context, create VAO,
            set up ArcballCamera. (Render loop: see report Part 4.)
      6. Print summary: equilibrium VQE energy, bond length, FCI energy, |VQE-FCI|.
         Equilibrium index = np.argmin(pes_data.energies_vqe).

    Args:
        args: Parsed command-line arguments from parse_args().
    """
    # Import modules here to allow the script to show help without dependencies
    from .molecule import create_h2, create_lih, create_h2o
    from .pes_scanner import scan_pes
    from .benchmarker import run_benchmark, print_benchmark_table
    
    # 1. Parse bond range
    parts = args.bond_range.split(",")
    start, stop, n_points = float(parts[0]), float(parts[1]), int(parts[2])
    bond_lengths = np.linspace(start, stop, n_points)
    
    # 2. Select molecule factory
    molecule_factories = {
        "h2": create_h2,
        "lih": create_lih,
        "h2o": create_h2o
    }
    factory = molecule_factories[args.molecule]
    
    print("=" * 60)
    print(f"VQE Potential Energy Surface Explorer")
    print("=" * 60)
    print(f"Molecule:     {args.molecule.upper()}")
    print(f"Basis:        {args.basis}")
    print(f"Bond range:   {start:.2f} - {stop:.2f} Å ({n_points} points)")
    print(f"Ansatz:       {args.ansatz}")
    print(f"Optimizer:    {args.optimizer}")
    print(f"Max iter:     {args.maxiter}")
    print(f"Benchmark:    {args.benchmark}")
    print("=" * 60)
    
    # 3. Run PES scan
    print("\n[Phase 1] Running PES scan...")
    pes_data = scan_pes(
        molecule_factory=factory,
        bond_lengths=bond_lengths,
        basis=args.basis,
        mapping="jordan_wigner",
        ansatz_type=args.ansatz,
        optimizer=args.optimizer,
        maxiter=args.maxiter,
        use_prev_params=True,
        compute_classical=args.benchmark
    )
    
    # 4. Benchmarking
    if args.benchmark:
        print("\n[Phase 2] Running benchmark analysis...")
        bench_result = run_benchmark(pes_data)
        table = print_benchmark_table(bench_result)
        print("\nBenchmark Results:")
        print(table)
    
    # 5. Visualization (optional)
    if args.render:
        print("\n[Phase 3] Preparing 3D visualization...")
        _run_visualization(args, factory)
    
    # 6. Summary
    eq_idx = np.argmin(pes_data.energies_vqe)
    eq_bond_length = pes_data.bond_lengths[eq_idx]
    eq_vqe_energy = pes_data.energies_vqe[eq_idx]
    eq_fci_energy = pes_data.energies_fci[eq_idx]
    error = abs(eq_vqe_energy - eq_fci_energy)
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Equilibrium bond length: {eq_bond_length:.4f} Å")
    print(f"VQE energy:              {eq_vqe_energy:.8f} Ha")
    print(f"FCI energy:              {eq_fci_energy:.8f} Ha")
    print(f"|VQE - FCI|:             {error:.2e} Ha")
    print(f"Chemical accuracy:       {'YES' if error < 1.6e-3 else 'NO'} (< 1.6 mHa)")
    print("=" * 60)
    
    # Plot PES curves
    _plot_pes(args, pes_data)

    # Save results
    _save_results(args, pes_data)


def _plot_pes(args: argparse.Namespace, pes_data) -> None:
    """Generate and save PES curve plot (E vs R for all methods)."""
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    R = pes_data.bond_lengths

    ax.plot(R, pes_data.energies_fci, "k-", linewidth=2, label="FCI (exact)")
    ax.plot(R, pes_data.energies_vqe, "o-", color="#2196F3", markersize=4,
            linewidth=1.5, label="VQE (UCCSD)")
    ax.plot(R, pes_data.energies_hf, "--", color="#F44336", linewidth=1.5,
            label="Hartree\u2013Fock")

    if pes_data.energies_cisd is not None:
        ax.plot(R, pes_data.energies_cisd, "s-", color="#4CAF50", markersize=3,
                linewidth=1, label="CISD")
    if pes_data.energies_ccsd_t is not None:
        ax.plot(R, pes_data.energies_ccsd_t, "^-", color="#FF9800", markersize=3,
                linewidth=1, label="CCSD(T)")

    ax.set_xlabel("Bond Length (\u00c5)")
    ax.set_ylabel("Energy (Hartree)")
    ax.set_title(f"Potential Energy Surface \u2014 {args.molecule.upper()} / {args.basis}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    save_path = output_dir / f"pes_{args.molecule}_{args.basis}.pdf"
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPES plot saved to: {save_path}")


def _run_visualization(args: argparse.Namespace, factory) -> None:
    """Run the 3D orbital visualization."""
    from .classical_hf import run_hartree_fock
    from .orbital_grid import evaluate_orbital_on_grid
    from .isosurface import extract_dual_isosurface
    from .renderer import (init_opengl_context, create_orbital_vao,
                           ArcballCamera, load_shader_from_file)
    import glfw
    from OpenGL.GL import (
        glViewport, glClearColor, glClear, glUseProgram, 
        glUniformMatrix4fv, glGetUniformLocation, glUniform3f,
        glBindVertexArray, glDrawElements,
        GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, GL_FALSE,
        GL_TRIANGLES, GL_UNSIGNED_INT
    )
    
    ANGSTROM_TO_BOHR = 1.8897259886
    
    # Build equilibrium molecule
    mol_spec = factory()
    hf_result = run_hartree_fock(mol_spec)
    
    # Get atom coordinates in Bohr
    atom_coords = np.array([atom.position for atom in mol_spec.atoms]) * ANGSTROM_TO_BOHR
    atom_symbols = [atom.symbol for atom in mol_spec.atoms]
    
    # Evaluate MO on grid
    print(f"  Evaluating MO {args.orbital} on grid...")
    orbital_data = evaluate_orbital_on_grid(
        mo_coefficients=hf_result.mo_coefficients,
        orbital_index=args.orbital,
        atom_coords=atom_coords,
        atom_symbols=atom_symbols,
        basis_name=mol_spec.basis,
        grid_extent=5.0,
        grid_spacing=0.15
    )
    
    # Extract isosurface
    print("  Extracting isosurface...")
    isovalue = 0.05
    vol = orbital_data.values
    # Clamp isovalue to be within the actual volume range so both lobes are attempted
    max_abs = float(np.abs(vol).max())
    if max_abs < isovalue:
        isovalue = max_abs * 0.3
        print(f"  Isovalue clamped to {isovalue:.4f} (volume max |val| = {max_abs:.4f})")
    mesh = extract_dual_isosurface(
        volume=orbital_data.values,
        isovalue=isovalue,
        grid_origin=orbital_data.grid_origin,
        grid_spacing=orbital_data.grid_spacing
    )
    if len(mesh.faces) == 0:
        print("  Warning: isosurface is empty — no surface to render.")
        return
    
    print(f"  Mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} triangles")
    
    # Initialize OpenGL
    window = init_opengl_context(1400, 900, f"VQE Orbital Viewer - MO {args.orbital}")
    
    # Load shaders
    shader_dir = Path(__file__).parent / "shaders"
    shader = load_shader_from_file(str(shader_dir), "orbital.vert", "orbital.frag")
    
    # Create VAO
    vao, vbo, ebo, n_indices = create_orbital_vao(mesh)
    
    # Camera setup — view from the side so the bond axis (z) is horizontal
    camera = ArcballCamera(
        position=np.array([12.0, 4.0, 0.0]),
        target=np.array([0.0, 0.0, 0.0])
    )
    
    # Mouse state
    mouse_state = {"last_x": 0, "last_y": 0, "dragging": False}
    
    def cursor_callback(window, xpos, ypos):
        if mouse_state["dragging"]:
            dx = xpos - mouse_state["last_x"]
            dy = ypos - mouse_state["last_y"]
            camera.on_mouse_drag(dx, dy)
        mouse_state["last_x"] = xpos
        mouse_state["last_y"] = ypos
    
    def mouse_button_callback(window, button, action, mods):
        if button == glfw.MOUSE_BUTTON_LEFT:
            mouse_state["dragging"] = (action == glfw.PRESS)
    
    def scroll_callback(window, xoffset, yoffset):
        camera.on_scroll(yoffset)
    
    glfw.set_cursor_pos_callback(window, cursor_callback)
    glfw.set_mouse_button_callback(window, mouse_button_callback)
    glfw.set_scroll_callback(window, scroll_callback)
    
    print("  Starting render loop (press ESC to exit)...")
    
    while not glfw.window_should_close(window):
        glfw.poll_events()
        
        if glfw.get_key(window, glfw.KEY_ESCAPE) == glfw.PRESS:
            break
        
        # Get window size
        width, height = glfw.get_framebuffer_size(window)
        if height == 0:
            height = 1
        aspect = width / height
        
        # Clear
        glViewport(0, 0, width, height)
        glClearColor(0.1, 0.1, 0.15, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        # Matrices
        view = camera.get_view_matrix()
        proj = camera.get_projection_matrix(aspect)
        model = np.eye(4, dtype=np.float32)
        
        # Draw orbital
        glUseProgram(shader)
        glUniformMatrix4fv(glGetUniformLocation(shader, "uModel"), 1, GL_FALSE, model.flatten())
        glUniformMatrix4fv(glGetUniformLocation(shader, "uView"), 1, GL_FALSE, view.flatten())
        glUniformMatrix4fv(glGetUniformLocation(shader, "uProjection"), 1, GL_FALSE, proj.flatten())
        glUniform3f(glGetUniformLocation(shader, "uLightPos"), 10.0, 10.0, 10.0)
        glUniform3f(glGetUniformLocation(shader, "uViewPos"), *camera.position)
        
        glBindVertexArray(vao)
        glDrawElements(GL_TRIANGLES, n_indices, GL_UNSIGNED_INT, None)
        glBindVertexArray(0)
        
        glfw.swap_buffers(window)
    
    glfw.terminate()


def _save_results(args: argparse.Namespace, pes_data) -> None:
    """Save PES scan results to JSON file."""
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {
        "molecule": args.molecule,
        "basis": args.basis,
        "ansatz": args.ansatz,
        "optimizer": args.optimizer,
        "maxiter": args.maxiter,
        "bond_lengths": pes_data.bond_lengths.tolist(),
        "energies_vqe": pes_data.energies_vqe.tolist(),
        "energies_hf": pes_data.energies_hf.tolist(),
        "energies_fci": pes_data.energies_fci.tolist(),
    }
    
    if pes_data.energies_cisd is not None:
        results["energies_cisd"] = pes_data.energies_cisd.tolist()
    if pes_data.energies_ccsd_t is not None:
        results["energies_ccsd_t"] = pes_data.energies_ccsd_t.tolist()
    
    output_file = output_dir / f"pes_{args.molecule}_{args.basis}.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")


def main() -> None:
    """Entry point: parse args and run the pipeline."""
    args = parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()
