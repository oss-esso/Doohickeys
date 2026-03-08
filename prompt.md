# Prompt for Claude Opus 4.6 — Report 1
## QAOA Energy Landscape Explorer: Full Research & Implementation Report

---

You are an expert in quantum computing, variational algorithms, and scientific software engineering. I am a quantum applications researcher preparing a one-day coding project that is both research-worthy and interview-demonstrable. I need a **thorough, publication-quality research report** on building a QAOA Energy Landscape Explorer from first principles to real hardware execution, with an interactive OpenGL visualization layer.

**Do not simplify or omit technical detail. Write for a quantum computing researcher, not a beginner.**

---

## What I Need

Produce a single, exhaustive research report structured exactly as follows:

---

### Part 1 — Theoretical Foundations

1. **Combinatorial Optimization & the Ising Model**
   - Derive the Ising Hamiltonian from first principles
   - Show the explicit QUBO ↔ Ising equivalence with worked example on MaxCut for a 4-node graph
   - Define the energy landscape: what it is geometrically, why it matters for optimization

2. **QAOA: Full Derivation**
   - Derive the QAOA circuit from the Quantum Approximate Optimization framework (Farhi et al.)
   - Define cost Hamiltonian H_C and mixer Hamiltonian H_B
   - Write out the full circuit unitary U(β,γ) = ∏ U_B(β_i) U_C(γ_i) with explicit gate decompositions
   - Explain the p→∞ limit and connection to adiabatic quantum computing
   - Derive the expectation value ⟨β,γ|H_C|β,γ⟩ analytically for p=1 on a simple graph

3. **The (β, γ) Parameter Landscape**
   - Explain why the landscape is periodic (give the exact periods)
   - Describe the known landscape structure: concentration of measure, symmetry properties, good initial parameter strategies (INTERP, FOURIER heuristics)
   - Define and explain barren plateaus formally: variance of gradient ∂⟨H_C⟩/∂θ and its scaling
   - Describe the phase transition behaviour with local field strength (reference the 2025 *Quantum* journal result on anti-ferromagnetic, critical, and trivial regimes)

4. **Classical Optimization of QAOA Parameters**
   - Parameter-shift rule: full derivation for ∂⟨H⟩/∂θ
   - Compare COBYLA, L-BFGS-B, SPSA, Adam for this landscape — when to use each
   - Explain the gradient-free vs gradient-based tradeoff on noisy hardware

---

### Part 2 — Software Architecture

Provide a **full module-by-module code architecture** with:
- File/folder structure tree
- For every module: its purpose, its inputs/outputs, the key functions with their signatures and docstrings (no implementation, just skeletons with full type hints and docstrings explaining the logic)
- Data flow diagram in ASCII art showing how modules connect

Modules must include:
1. `graph_generator.py` — graph construction, adjacency to Ising conversion
2. `qaoa_circuit.py` — circuit construction using PennyLane
3. `optimizer.py` — parameter landscape sampling + classical optimization loop
4. `landscape_sampler.py` — grid sampling of (β, γ) space for fixed p
5. `noise_model.py` — depolarizing + thermal relaxation noise model via PennyLane's `default.mixed` device; use QuTiP Lindblad master equation for deeper noise characterisation and channel comparison
6. `hardware_runner.py` — IBM Quantum / IQM submission via PennyLane's `qml.device("qiskit.ibmq", ...)` plugin or IQM Cloud PennyLane plugin; no raw Qiskit Runtime
7. `error_mitigation.py` — ZNE (zero-noise extrapolation) implementation
8. `renderer.py` — OpenGL 3D landscape surface renderer
9. `main.py` — orchestration and CLI

---

### Part 3 — Implementation Guidelines (Instructional Pseudocode)

For each module listed above, provide **instructional pseudocode** — not full working code, but detailed enough that a skilled developer can implement it without ambiguity. This means:
- All major loops written out
- All matrix/tensor shapes annotated in comments
- All non-obvious algorithmic choices explained inline
- Key PennyLane and QuTiP calls named explicitly (assume full fluency with NumPy, SciPy, JAX, and all standard math libraries — do not explain these)
- Edge cases flagged

Pay special attention to:
- How to construct the cost Hamiltonian as a `qml.Hamiltonian` object from an adjacency matrix
- How to implement the parameter-shift rule manually (not via autograd) for hardware compatibility
- How to do grid sampling of the landscape efficiently (vectorised)
- How to pass landscape data as a VBO to OpenGL for real-time surface rendering with a colour map (energy → hue)
- How to use PennyLane's `qml.ExpvalCost` / `qml.execute` with the IBM backend plugin for hardware expectation values

---

### Part 4 — OpenGL Visualization Design

1. Describe the full rendering pipeline for a 3D (β, γ, ⟨H_C⟩) surface:
   - Vertex buffer layout
   - Shader code (GLSL) for the surface with energy-based colouring
   - Camera controls (arcball rotation)
   - How to update the surface live as optimization progresses (streaming VBO updates)

2. Explain how to overlay the optimizer trajectory as a 3D polyline on the surface

3. Describe a side panel showing: current circuit depth p, current ⟨H_C⟩, approximation ratio, gradient norm

---

### Part 5 — Hardware Execution Plan

1. Walk through the full IBM Quantum submission pipeline using PennyLane's IBM plugin:
   - Transpilation for a specific backend (e.g. ibm_brisbane) via PennyLane-Qiskit device
   - How to choose the right qubits using coupling map analysis
   - Setting shots, using `qml.sample` vs `qml.expval` on hardware devices
   - How to extract expectation values from shot distributions and handle readout error

2. Error mitigation strategy:
   - ZNE with linear + Richardson extrapolation: full method
   - When to use Pauli twirling as a precursor
   - Realistic expectation: what fidelity to expect at p=1, p=2 on 10 qubits

3. A benchmarking checklist: what metrics to record, how to plot noiseless vs noisy vs error-mitigated results

---

### Part 6 — Research Contribution Angle

1. Describe 3 concrete ways this project could be extended into a publishable paper
2. Point to the most relevant 2024–2025 papers to cite and relate to
3. Suggest what graph families or problem types would give the most interesting landscape structure to visualize

---

**Formatting instructions:**
- Use LaTeX math notation throughout (inline $ $ and block $$ $$)
- Use Python code blocks for all pseudocode
- Use headers and subheaders exactly as structured above
- Be exhaustive — do not summarise or skip steps
- Length is not a concern; completeness is


# Prompt for Claude Opus 4.6 — Report 2
## VQE Potential Energy Surface & 3D Molecular Orbital Renderer: Full Research & Implementation Report

---

You are an expert in quantum chemistry, variational quantum algorithms, and scientific visualization. I am a quantum applications researcher preparing a one-day coding project that is both research-worthy and visually compelling for interviews. I need a **thorough, publication-quality research report** on building a VQE-based Potential Energy Surface (PES) explorer with real-time 3D molecular orbital isosurface rendering using OpenGL, from quantum chemistry theory all the way to real hardware execution.

**Do not simplify or omit technical detail. Write for a quantum computing researcher with some quantum chemistry background but not a specialist.**

---

## What I Need

Produce a single, exhaustive research report structured exactly as follows:

---

### Part 1 — Theoretical Foundations

1. **Quantum Chemistry from First Principles**
   - Write out the molecular Hamiltonian (Born-Oppenheimer approximation): nuclear and electronic terms
   - Define the electronic structure problem: why it's exponentially hard classically
   - Explain the second-quantized Hamiltonian: creation/annihilation operators, one- and two-electron integrals
   - Define molecular orbitals, basis sets (STO-3G, 6-31G*), and the Fock matrix
   - Explain Hartree-Fock (HF) theory and its failure near bond dissociation — why correlation energy matters

2. **Fermion-to-Qubit Mappings**
   - Derive the Jordan-Wigner transformation explicitly: how fermionic operators become Pauli strings
   - Explain the Parity mapping and its qubit-reduction advantages
   - Show a worked example: map the H₂ STO-3G Hamiltonian to a qubit operator explicitly (list all Pauli terms with coefficients)
   - Compare qubit counts and circuit depth implications for H₂ vs LiH vs H₂O

3. **Variational Quantum Eigensolver (VQE)**
   - Derive the variational principle: why ⟨ψ(θ)|H|ψ(θ)⟩ ≥ E_ground
   - Define the ansatz: explain UCCSD (Unitary Coupled Cluster Singles and Doubles) from the coupled-cluster wavefunction
   - Derive the UCCSD operator T = T₁ + T₂, the UCC ansatz e^(T-T†)|HF⟩, and how it's Trotterised into gates
   - Explain hardware-efficient ansatz as an alternative and its tradeoffs
   - Describe the full VQE loop: state preparation → measurement → classical optimization → convergence

4. **Potential Energy Surface (PES)**
   - Define the PES: E(R) as a function of nuclear coordinates
   - Explain what the PES reveals: equilibrium geometry, dissociation energy, transition states
   - Describe the specific challenge of PES accuracy: why FCI is the gold standard and where VQE sits relative to CCSD(T), CISD, MP2
   - For H₂: give the known exact dissociation curve and what chemical accuracy (1.6 mHartree) means

5. **Molecular Orbitals and Wavefunction Visualization**
   - Define molecular orbital: LCAO expansion, coefficients from HF diagonalisation
   - Explain how to evaluate ψ(r) on a 3D grid from basis function coefficients
   - Define electron density ρ(r) = Σ_i |ψ_i(r)|² and the one-particle reduced density matrix
   - Explain isosurface extraction: what isovalue to choose, why ±0.05 a.u. is conventional

---

### Part 2 — Software Architecture

Provide a **full module-by-module code architecture** with:
- File/folder structure tree
- For every module: its purpose, its inputs/outputs, the key functions with their signatures, full type hints, and detailed docstrings explaining the logic
- Data flow diagram in ASCII art showing how all modules connect end-to-end

Modules must include:
1. `molecule.py` — molecule definition, geometry specification, basis set selection
2. `classical_hf.py` — PySCF interface: run HF, extract integrals, MO coefficients, orbital energies
3. `hamiltonian_builder.py` — build second-quantized Hamiltonian from PySCF integrals, apply Jordan-Wigner / Parity mapping via `qml.qchem`; use QuTiP's `Qobj` for operator algebra validation and exact diagonalisation cross-checks
4. `ansatz.py` — UCCSD and hardware-efficient ansatz construction in PennyLane
5. `vqe_solver.py` — VQE optimization loop with gradient and gradient-free options
6. `pes_scanner.py` — PES sweep: loop over bond lengths, call VQE at each geometry, store results
7. `orbital_grid.py` — evaluate MO wavefunctions on a 3D Cartesian grid using basis functions
8. `isosurface.py` — marching cubes algorithm to extract isosurface mesh from volumetric data
9. `renderer.py` — OpenGL renderer: orbital isosurface + PES curve + bond geometry
10. `hardware_runner.py` — IBM Quantum / IonQ hardware submission via PennyLane's IBM or IonQ device plugins; no raw Qiskit Runtime
11. `benchmarker.py` — compare VQE, HF, CISD, CCSD(T), FCI across the PES
12. `main.py` — orchestration, CLI, config file parsing

---

### Part 3 — Implementation Guidelines (Instructional Pseudocode)

For each module listed above, provide **instructional pseudocode** — detailed enough that a skilled developer can implement without ambiguity:
- All major loops written out explicitly
- All array/tensor shapes annotated in comments (e.g., `# shape: (n_qubits, n_params)`)
- All non-obvious algorithmic choices explained inline
- Key PySCF, PennyLane (`qml.qchem`, `qml.grad`, `qml.device`), and QuTiP calls named explicitly (assume full fluency with NumPy, SciPy, JAX, and all standard math libraries — do not explain these)
- Edge cases and gotchas flagged

Pay special attention to:
- How to extract one- and two-electron integrals from PySCF and pass them into `qml.qchem.molecular_hamiltonian`
- How to construct the UCCSD excitation list via `qml.qchem.excitations` and map to parametrised `qml.DoubleExcitation` / `qml.SingleExcitation` gates
- How to use QuTiP `mesolve` for Lindblad noise simulation and compare against PennyLane `default.mixed`
- How to implement the marching cubes isosurface extraction (explain the lookup table approach)
- How to stream the isosurface mesh as a VAO/VBO into OpenGL and update it as bond length changes
- How to use the parameter-shift rule with PennyLane's `qml.grad` for UCCSD circuits

---

### Part 4 — OpenGL Visualization Design

Describe the complete visualization, which must include three components rendered simultaneously:

1. **Molecular Orbital Isosurface**
   - Full VAO/VBO layout for the isosurface mesh (vertices, normals, isovalue sign for +/- lobe colouring)
   - GLSL vertex + fragment shader: Phong lighting model, +lobe blue / −lobe red colouring with transparency
   - How to recompute and re-upload the mesh when the bond length slider changes (double-buffered VBO update)

2. **Potential Energy Surface Curve**
   - Render a 2D PES curve (bond length vs energy) as a 3D ribbon in world space
   - Colour-encode: HF in grey, CISD in green, VQE in blue, FCI/exact in gold
   - Show a moving marker at the current bond length

3. **Molecular Geometry**
   - Render atom spheres (van der Waals radii scaled) and bond cylinders
   - Animate bond stretching as the PES scan progresses

4. Camera, lighting, and UI layout description

---

### Part 5 — Hardware Execution Plan

1. Walk through submitting a VQE circuit for H₂ at equilibrium to IBM Quantum via PennyLane's IBM plugin:
   - Circuit transpilation and qubit layout selection through the PennyLane-Qiskit device
   - Using `qml.expval` on hardware with Pauli grouping to minimise measurement circuits (explain commuting clique partitioning)
   - Shot budget: how to estimate required shots for chemical accuracy

2. Error mitigation strategy for chemistry:
   - Why ZNE is well-suited for VQE energy estimation
   - Symmetry verification as an additional filter (particle number, spin)
   - What to expect: typical VQE energy error on hardware for H₂ at p=1 vs p=2 layers

3. Comparison checklist: noiseless simulation vs noisy simulation vs hardware — what metrics to record and plot

---

### Part 6 — Benchmarking & Validation

1. How to validate your VQE PES against PySCF's FCI solver (give the exact PySCF calls)
2. How to measure and plot:
   - Correlation energy recovered: (E_VQE − E_HF) / (E_FCI − E_HF)
   - Fidelity of VQE state with exact ground state
   - Circuit depth vs energy accuracy tradeoff
3. A table template for recording results at 5 bond lengths across 5 methods

---

### Part 7 — Research Contribution Angle

1. Describe 3 concrete ways this project extends into a publishable result (e.g., SHARC-VQE partitioned Hamiltonian benchmarking, adaptive ansatz PES, noise-aware VQE)
2. Point to the most relevant 2024–2025 papers in VQE and quantum chemistry simulation to cite
3. Explain what molecules beyond H₂ are tractable on today's hardware and why

---

**Formatting instructions:**
- Use LaTeX math notation throughout (inline $ $ and block $$ $$)
- Use Python code blocks for all pseudocode
- Use headers and subheaders exactly as structured above
- Be exhaustive — do not summarise or skip steps
- Length is not a concern; completeness is