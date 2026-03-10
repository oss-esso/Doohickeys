# AI Prompt for Generating VQE PES Report

## Context

You are writing a comprehensive technical report documenting the VQE (Variational Quantum Eigensolver) Potential Energy Surface Explorer, following the style and structure of `report1_qaoa_worked_example.tex`. The report should be self-contained with full mathematical derivations.

---

## Prompt

Generate a LaTeX report titled **"VQE Molecular Energy Surface Explorer: A Complete Mathematical Walkthrough for H₂ using the UCCSD Ansatz"** with the following structure:

### 1. Introduction (1 page)
- Introduce VQE as a hybrid quantum-classical algorithm for finding ground-state energies
- Motivation: near-term quantum advantage for chemistry problems
- Overview of the H₂ molecule as the simplest test case
- Outline what the report will cover

### 2. Problem Formulation (2 pages)

#### 2.1 The Electronic Structure Problem
- Write the molecular Hamiltonian in second quantization:
  $$\hat{H} = \sum_{pq} h_{pq} \hat{a}_p^\dagger \hat{a}_q + \frac{1}{2}\sum_{pqrs} g_{pqrs} \hat{a}_p^\dagger \hat{a}_r^\dagger \hat{a}_s \hat{a}_q + E_{\text{nuc}}$$
- Define one-electron integrals $h_{pq}$ and two-electron integrals $g_{pqrs}$
- Explain the Born-Oppenheimer approximation

#### 2.2 The H₂ Molecule
- Geometry: two H atoms at distance $R$ along the z-axis
- STO-3G basis: 2 spatial orbitals → 4 spin-orbitals
- Number of electrons: 2
- Hilbert space dimension: $\binom{4}{2} = 6$ configurations

#### 2.3 Classical Reference Methods
- Explain Hartree-Fock (HF) as mean-field approximation
- Define electron correlation energy: $E_{\text{corr}} = E_{\text{exact}} - E_{\text{HF}}$
- Mention CISD, CCSD(T), and FCI as benchmarks

### 3. Qubit Hamiltonian Construction (2 pages)

#### 3.1 Jordan-Wigner Transformation
- Map fermionic operators to Pauli strings:
  $$\hat{a}_j^\dagger \rightarrow \frac{1}{2}(X_j - iY_j) \otimes Z_{j-1} \otimes \cdots \otimes Z_0$$
- Show that the H₂ Hamiltonian in STO-3G becomes a sum of ~15 Pauli terms
- Present the explicit Hamiltonian with coefficients

#### 3.2 Parity Transformation (Alternative)
- Briefly describe the parity encoding
- Note: can reduce qubit count via Z₂ symmetries

#### 3.3 Worked Example: H₂ at R = 0.74 Å
- Show full derivation of Hamiltonian terms
- Table of Pauli strings and coefficients
- Verification via exact diagonalization (QuTiP)

### 4. VQE Circuit Design (2 pages)

#### 4.1 The VQE Algorithm
- State preparation: $|\psi(\vec{\theta})\rangle = U(\vec{\theta})|\phi_{\text{ref}}\rangle$
- Energy evaluation: $E(\vec{\theta}) = \langle\psi(\vec{\theta})|\hat{H}|\psi(\vec{\theta})\rangle$
- Classical optimization loop
- Measurement: decompose $\hat{H}$ into measurable Pauli groups

#### 4.2 UCCSD Ansatz
- Unitary Coupled Cluster with Singles and Doubles:
  $$U_{\text{UCCSD}} = e^{\hat{T} - \hat{T}^\dagger}$$
- Define excitation operators:
  - Singles: $\hat{T}_1 = \sum_{i,a} \theta_i^a \hat{a}_a^\dagger \hat{a}_i$
  - Doubles: $\hat{T}_2 = \sum_{i<j,a<b} \theta_{ij}^{ab} \hat{a}_a^\dagger \hat{a}_b^\dagger \hat{a}_j \hat{a}_i$
- For H₂: 1 double excitation only (HOMO→LUMO)
  $$U = e^{\theta(\hat{a}_2^\dagger\hat{a}_3^\dagger\hat{a}_1\hat{a}_0 - h.c.)}$$

#### 4.3 Circuit Implementation
- Show PennyLane `qml.DoubleExcitation` gate
- Circuit diagram for H₂ UCCSD (4 qubits, 1 parameter)
- Approximate gate count and depth

### 5. Classical Optimization (1.5 pages)

#### 5.1 Gradient-Based Optimizers
- Parameter-shift rule for gradient computation
- Adam optimizer: adaptive learning rates
- Gradient descent with momentum

#### 5.2 Gradient-Free Optimizers
- COBYLA: constrained optimization by linear approximation
- SPSA: simultaneous perturbation stochastic approximation
- Comparison for noisy evaluations

#### 5.3 Convergence Analysis
- Plot: VQE energy vs. iteration for each optimizer
- Table: Final energies and iteration counts
- Discussion of local minima and barren plateaus

### 6. Potential Energy Surface Scan (2 pages)

#### 6.1 Methodology
- Scan bond length: $R \in [0.3, 3.0]$ Å with 20 points
- Warm-start: use previous optimal parameters as initial guess
- Parallel computation of HF, CISD, CCSD(T), FCI at each $R$

#### 6.2 Results
- Plot: $E(R)$ for HF, VQE, CISD, CCSD(T), FCI
- Table: Selected energies at equilibrium, stretched, and compressed geometries
- VQE achieves chemical accuracy (< 1.6 mHa) everywhere

#### 6.3 Correlation Energy Recovery
- Define: $\eta = \frac{E_{\text{VQE}} - E_{\text{HF}}}{E_{\text{FCI}} - E_{\text{HF}}}$
- Plot $\eta(R)$ showing VQE recovers >99% correlation at all geometries
- Discuss behavior at dissociation ($R \rightarrow \infty$): RHF fails, VQE handles it

### 7. Molecular Orbital Visualization (1 page)

#### 7.1 Orbital Grid Evaluation
- Expand MO as linear combination of AOs
- Evaluate on 3D grid using PySCF `mol.eval_gto`
- Grid parameters: spacing = 0.1 Bohr, extent = 5 Bohr

#### 7.2 Isosurface Rendering
- Marching cubes algorithm for isosurface extraction
- Dual surfaces: $+\psi$ and $-\psi$ lobes
- Phong shading for 3D appearance
- Show bonding (σ) and antibonding (σ*) orbitals for H₂

### 8. Discussion (1 page)

#### 8.1 Sources of Error
- Basis set incompleteness (STO-3G is minimal)
- Ansatz expressibility (UCCSD may miss higher excitations for larger molecules)
- Optimization convergence (local minima)
- Measurement shot noise (not considered in noiseless simulation)

#### 8.2 Comparison with Classical Methods
- VQE matches FCI for H₂ (UCCSD is exact for 2-electron systems)
- Scaling: VQE $O(N^4)$ parameters vs FCI $O(N!)$ states
- When VQE wins: larger systems where FCI is intractable

#### 8.3 Extensions
- Larger molecules (LiH, H₂O)
- Active space reduction
- Hardware execution with error mitigation
- Adaptive VQE for automatic ansatz construction

### 9. Conclusion (0.5 page)
- Summary of key results
- VQE accurately reproduces H₂ PES
- Code available for exploration

### Appendices

#### A. Code Listings
- Key Python code for Hamiltonian construction
- UCCSD ansatz implementation
- VQE optimization loop

#### B. Mathematical Details
- Full Jordan-Wigner derivation
- UCCSD parameter count formula
- Gradient via parameter-shift rule

---

## Data to Include

Run the following command to generate raw data:

```bash
# In WSL
source .venv/bin/activate
python -m report2_scripts.main \
    --molecule h2 \
    --basis sto-3g \
    --bond-range "0.3,3.0,20" \
    --ansatz uccsd \
    --optimizer Adam \
    --maxiter 200 \
    --benchmark \
    --output output/report_data
```

Extract from `output/report_data/pes_h2_sto-3g.json`:
- `bond_lengths`: array of R values
- `energies_vqe`: VQE results
- `energies_hf`: Hartree-Fock results
- `energies_fci`: Exact (FCI) results
- `energies_cisd`, `energies_ccsd_t`: Classical post-HF methods

---

## Style Guidelines

1. **Mathematical rigor**: Show all derivation steps explicitly
2. **Self-contained**: Reader should understand without external references
3. **Visual**: Include figures for PES curves, correlation recovery, orbitals
4. **Code snippets**: Show key PennyLane/PySCF code in `lstlisting`
5. **Tables**: Use `booktabs` for professional formatting
6. **Length**: ~15-20 pages total

---

## Example Sections to Emulate

From `report1_qaoa_worked_example.tex`:
- Section 3: Ising Hamiltonian (explicit term derivation)
- Section 4: QAOA Circuit Construction (gate-by-gate)
- Table 1: Eigenvalue verification
- Section 5: Energy Landscape (parameter space visualization)

---

## LaTeX Template Headers

```latex
\documentclass[11pt,a4paper]{article}

% Packages
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{physics}
\usepackage{graphicx}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue]{hyperref}
\usepackage{listings}
\usepackage{booktabs}
\usepackage{float}
\usepackage{geometry}
\geometry{margin=1in}

% Code listing style
\lstset{
    language=Python,
    basicstyle=\ttfamily\footnotesize,
    keywordstyle=\color{blue},
    commentstyle=\color{gray},
    breaklines=true
}

\title{VQE Molecular Energy Surface Explorer:\\
    A Complete Mathematical Walkthrough\\
    for H$_2$ using the UCCSD Ansatz}
\author{Generated Report}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\newpage
% ... content ...
\end{document}
```
