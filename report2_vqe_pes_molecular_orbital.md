# VQE Potential Energy Surface & 3D Molecular Orbital Renderer: Full Research & Implementation Report

---

## Part 1 — Theoretical Foundations

### 1. Quantum Chemistry from First Principles

#### 1.1 The Molecular Hamiltonian Under the Born-Oppenheimer Approximation

The full non-relativistic Hamiltonian for a molecule with $N_e$ electrons and $N_n$ nuclei is (in atomic units, $\hbar = m_e = e = 4\pi\epsilon_0 = 1$):

$$\hat{H}_{\text{mol}} = \hat{T}_n + \hat{T}_e + \hat{V}_{nn} + \hat{V}_{en} + \hat{V}_{ee}$$

where:

$$\hat{T}_n = -\sum_{A=1}^{N_n} \frac{1}{2M_A}\nabla_A^2 \qquad \text{(nuclear kinetic energy)}$$

$$\hat{T}_e = -\sum_{i=1}^{N_e} \frac{1}{2}\nabla_i^2 \qquad \text{(electronic kinetic energy)}$$

$$\hat{V}_{nn} = \sum_{A<B} \frac{Z_A Z_B}{|\mathbf{R}_A - \mathbf{R}_B|} \qquad \text{(nuclear-nuclear repulsion)}$$

$$\hat{V}_{en} = -\sum_{i=1}^{N_e}\sum_{A=1}^{N_n} \frac{Z_A}{|\mathbf{r}_i - \mathbf{R}_A|} \qquad \text{(electron-nuclear attraction)}$$

$$\hat{V}_{ee} = \sum_{i<j} \frac{1}{|\mathbf{r}_i - \mathbf{r}_j|} \qquad \text{(electron-electron repulsion)}$$

Here $M_A$ and $Z_A$ are the mass and atomic number of nucleus $A$, $\mathbf{R}_A$ are nuclear coordinates, and $\mathbf{r}_i$ are electron coordinates.

The **Born-Oppenheimer (BO) approximation** exploits the mass ratio $m_e / M_A \sim 10^{-3}$–$10^{-5}$: nuclei are vastly heavier than electrons, so electrons respond instantaneously to nuclear motion. Under this separation, the total wavefunction factorises:

$$\Psi(\mathbf{r}, \mathbf{R}) \approx \psi_{\text{el}}(\mathbf{r}; \mathbf{R}) \cdot \chi_{\text{nuc}}(\mathbf{R})$$

where $\psi_{\text{el}}$ depends parametrically on nuclear positions $\mathbf{R}$. The **electronic Hamiltonian** at fixed nuclear geometry is:

$$\hat{H}_{\text{el}}(\mathbf{R}) = \hat{T}_e + \hat{V}_{en} + \hat{V}_{ee} + \hat{V}_{nn}$$

(the nuclear repulsion $\hat{V}_{nn}$ is a constant at fixed $\mathbf{R}$). Solving the electronic Schrödinger equation:

$$\hat{H}_{\text{el}}|\psi_k(\mathbf{R})\rangle = E_k(\mathbf{R})|\psi_k(\mathbf{R})\rangle$$

yields the electronic energy $E_k(\mathbf{R})$ as a function of nuclear coordinates — this is the **potential energy surface (PES)**.

#### 1.2 The Electronic Structure Problem: Exponential Complexity

The electronic wavefunction $\psi_{\text{el}}(\mathbf{r}_1, \sigma_1, \ldots, \mathbf{r}_{N_e}, \sigma_{N_e})$ is a function of $3N_e$ spatial coordinates and $N_e$ spin variables. For fermions, it must be antisymmetric under particle exchange.

Given a basis of $M$ spin-orbitals $\{\phi_p\}$, the Hilbert space of $N_e$ electrons occupying $M$ orbitals has dimension $\binom{M}{N_e}$. For $M = 100$ orbitals and $N_e = 20$ electrons, $\binom{100}{20} \approx 5.4 \times 10^{20}$ — this is the Full Configuration Interaction (FCI) dimension. Storing a single FCI vector at double precision would require $\sim 4 \times 10^9$ TB. This exponential scaling is why classical exact diagonalisation is intractable beyond $\sim 20$ orbitals and is the fundamental motivation for quantum computing approaches.

#### 1.3 Second-Quantized Hamiltonian

In second quantisation, we introduce creation $\hat{a}_p^\dagger$ and annihilation $\hat{a}_p$ operators for each spin-orbital $\phi_p$, satisfying the fermionic anti-commutation relations:

$$\{\hat{a}_p, \hat{a}_q^\dagger\} = \delta_{pq}, \qquad \{\hat{a}_p, \hat{a}_q\} = 0, \qquad \{\hat{a}_p^\dagger, \hat{a}_q^\dagger\} = 0$$

The electronic Hamiltonian becomes:

$$\hat{H}_{\text{el}} = \sum_{p,q} h_{pq}\, \hat{a}_p^\dagger \hat{a}_q + \frac{1}{2}\sum_{p,q,r,s} h_{pqrs}\, \hat{a}_p^\dagger \hat{a}_q^\dagger \hat{a}_s \hat{a}_r + E_{\text{nuc}}$$

where:

- **One-electron integrals** (kinetic energy + electron-nuclear attraction):
$$h_{pq} = \int \phi_p^*(\mathbf{x})\left(-\frac{1}{2}\nabla^2 - \sum_A \frac{Z_A}{|\mathbf{r} - \mathbf{R}_A|}\right)\phi_q(\mathbf{x})\,d\mathbf{x}$$

- **Two-electron integrals** (electron-electron repulsion, in chemists' notation):
$$h_{pqrs} = \iint \frac{\phi_p^*(\mathbf{x}_1)\phi_r(\mathbf{x}_1)\phi_q^*(\mathbf{x}_2)\phi_s(\mathbf{x}_2)}{|\mathbf{r}_1 - \mathbf{r}_2|}\,d\mathbf{x}_1\,d\mathbf{x}_2$$

Here $\mathbf{x} = (\mathbf{r}, \sigma)$ denotes combined spatial and spin coordinates. The two-electron integrals in **physicists' notation** (antisymmetrised) are:

$$\langle pq || rs \rangle = h_{pqrs} - h_{pqsr}$$

$E_{\text{nuc}} = \hat{V}_{nn}$ is the constant nuclear repulsion energy at the given geometry.

The number of one-electron integrals scales as $O(M^2)$ and two-electron integrals as $O(M^4)$. These integrals are computed classically by quantum chemistry packages (PySCF, Psi4, etc.) and serve as input to the quantum algorithm.

#### 1.4 Molecular Orbitals, Basis Sets, and the Fock Matrix

**Molecular orbitals (MOs)** are expressed as linear combinations of atomic orbitals (LCAO):

$$\phi_p(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu p}\, \chi_\mu(\mathbf{r})$$

where $\{\chi_\mu\}$ is the **atomic orbital (AO) basis set** and $C_{\mu p}$ are the MO coefficients obtained by solving the Hartree-Fock equations.

**Basis sets** approximate the true atomic orbitals with Gaussian-type orbitals (GTOs):

- **STO-3G (Minimal basis):** Each Slater-type orbital (STO) is approximated by 3 Gaussian primitives. For H₂: 2 basis functions (one per atom). For H₂O: 7 basis functions (1s, 2s, 2p×3 on O; 1s on each H). This is the smallest possible basis — qualitatively correct but quantitatively poor.

- **6-31G\* (Split-valence + polarisation):** Valence orbitals are split into two functions (3+1 Gaussians), core orbitals use 6 Gaussians, and polarisation d-functions are added on heavy atoms. For H₂O: 19 basis functions. Significantly more accurate, especially for geometries and relative energies.

The **Fock matrix** $\mathbf{F}$ in the AO basis is:

$$F_{\mu\nu} = h_{\mu\nu}^{\text{core}} + \sum_{\lambda\sigma} P_{\lambda\sigma}\left[(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)\right]$$

where $h^{\text{core}}$ is the one-electron (core) Hamiltonian matrix, $P_{\lambda\sigma}$ is the density matrix, and $(\mu\nu|\lambda\sigma)$ are two-electron integrals in the AO basis. The Hartree-Fock equations (Roothaan-Hall equations) are the generalised eigenvalue problem:

$$\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\epsilon}$$

where $\mathbf{S}$ is the overlap matrix $S_{\mu\nu} = \langle\chi_\mu|\chi_\nu\rangle$ and $\boldsymbol{\epsilon}$ contains the orbital energies. This is solved self-consistently (SCF): the Fock matrix depends on the density matrix $\mathbf{P} = \mathbf{C}_{\text{occ}}\mathbf{C}_{\text{occ}}^T$, which depends on $\mathbf{C}$, which comes from diagonalising $\mathbf{F}$.

#### 1.5 Hartree-Fock Theory and Its Failure at Bond Dissociation

**Restricted Hartree-Fock (RHF)** constrains the wavefunction to be a single Slater determinant:

$$|\Phi_{\text{HF}}\rangle = \hat{a}_{1}^\dagger \hat{a}_{2}^\dagger \cdots \hat{a}_{N_e}^\dagger |0\rangle$$

This captures the mean-field electron-electron interaction (each electron moves in the average potential of all others) but misses **correlation energy** — the energy difference between the exact solution and the HF solution:

$$E_{\text{corr}} = E_{\text{exact}} - E_{\text{HF}} < 0$$

HF works well near equilibrium geometries where the single-determinant picture is qualitatively correct. However, it fails catastrophically at bond dissociation because:

1. **Static (strong) correlation:** At large internuclear separation, the correct wavefunction becomes a superposition of multiple determinants. For H₂ dissociation, the exact wavefunction at $R \to \infty$ is:

$$|\Psi\rangle \to \frac{1}{\sqrt{2}}(|\sigma_g\bar{\sigma}_g\rangle - |\sigma_u\bar{\sigma}_u\rangle)$$

(a superposition of doubly-occupied bonding and anti-bonding configurations). RHF forces a single determinant $|\sigma_g\bar{\sigma}_g\rangle$, which incorrectly mixes in ionic terms ($\text{H}^+\text{H}^-$) at dissociation, predicting an energy $\sim 0.25$ Hartree too high.

2. **Dynamic correlation:** Short-range electron-electron cusp conditions are never captured by a single determinant. This contributes a correlation energy typically $\sim 1\%$ of total energy but $\sim 100\%$ of chemical bond energies.

Post-HF methods (MP2, CISD, CCSD(T), FCI) recover correlation energy at increasing computational cost. VQE on a quantum computer can in principle reach FCI accuracy with polynomial quantum resources.

---

### 2. Fermion-to-Qubit Mappings

#### 2.1 Jordan-Wigner Transformation

The Jordan-Wigner (JW) transformation maps fermionic operators to qubit (Pauli) operators by encoding the occupation of each spin-orbital in a qubit and using $Z$-strings to enforce anti-commutation:

$$\hat{a}_p^\dagger \mapsto \frac{1}{2}(\hat{X}_p - i\hat{Y}_p) \otimes \prod_{q=0}^{p-1} \hat{Z}_q = \frac{1}{2}\hat{\sigma}_p^- \bigotimes_{q<p} \hat{Z}_q$$

$$\hat{a}_p \mapsto \frac{1}{2}(\hat{X}_p + i\hat{Y}_p) \otimes \prod_{q=0}^{p-1} \hat{Z}_q = \frac{1}{2}\hat{\sigma}_p^+ \bigotimes_{q<p} \hat{Z}_q$$

where $\hat{\sigma}_p^{\pm} = (\hat{X}_p \mp i\hat{Y}_p)/2$ are the qubit raising/lowering operators and the $Z$-string $\prod_{q<p}\hat{Z}_q$ provides the fermionic sign (parity of all lower-indexed orbitals).

**Verification of anti-commutation:** Consider $\hat{a}_p^\dagger \hat{a}_q + \hat{a}_q \hat{a}_p^\dagger$ for $p \neq q$. The $Z$-strings on the overlapping qubits between $p$ and $q$ introduce a sign flip, yielding:

$$\{\hat{a}_p, \hat{a}_q^\dagger\}_{\text{JW}} = \delta_{pq}\hat{I}$$

as required.

**Number operator:** $\hat{n}_p = \hat{a}_p^\dagger \hat{a}_p \mapsto \frac{1}{2}(\hat{I} - \hat{Z}_p)$, which has eigenvalue 0 on $|0\rangle$ (unoccupied) and 1 on $|1\rangle$ (occupied).

**One-body terms under JW:**

$$\hat{a}_p^\dagger \hat{a}_q \mapsto \frac{1}{4}\bigl[(\hat{X}_p\hat{X}_q + \hat{Y}_p\hat{Y}_q) + i(\hat{X}_p\hat{Y}_q - \hat{Y}_p\hat{X}_q)\bigr]\prod_{k=\min(p,q)+1}^{\max(p,q)-1}\hat{Z}_k$$

For $p = q$: $\hat{a}_p^\dagger\hat{a}_p = \frac{1}{2}(\hat{I} - \hat{Z}_p)$.

The $Z$-string between indices $p$ and $q$ makes the operator weight (number of Pauli terms) scale as $O(N)$ for distant orbitals. This is the principal cost of JW: non-local operators.

**Two-body terms under JW:** Each $\hat{a}_p^\dagger\hat{a}_q^\dagger\hat{a}_s\hat{a}_r$ term generates a sum of Pauli strings of weight up to $O(N)$. The total number of Pauli terms in the qubit Hamiltonian scales as $O(M^4)$.

#### 2.2 Parity Mapping

The **parity mapping** (Bravyi-Kitaev parity basis) encodes the parity of the occupation number rather than the occupation itself. Qubit $p$ stores the parity $\bigoplus_{q \leq p} n_q$ (sum modulo 2 of occupations up to orbital $p$).

The creation operator becomes:

$$\hat{a}_p^\dagger \mapsto \frac{1}{2}\hat{X}_p \bigotimes_{q \in U(p)} \hat{X}_q \cdot \bigotimes_{k \in P(p)} \hat{Z}_k - \frac{i}{2}\hat{Y}_p \bigotimes_{q \in U(p)} \hat{X}_q$$

where $U(p)$ and $P(p)$ are the "update" and "parity" sets determined by the binary tree structure.

**Qubit reduction advantage:** In the parity basis, the last qubit stores the total parity of the electron number, and the second-to-last stores the total spin parity. For systems with fixed electron number and fixed spin (which is always the case in chemistry), these two qubits are redundant and can be removed, reducing the qubit count by 2:

$$M \text{ spin-orbitals} \to M - 2 \text{ qubits (parity mapping with tapering)}$$

This is significant for small molecules: H₂ (STO-3G) goes from 4 qubits (JW) to 2 qubits (parity + tapering).

#### 2.3 Worked Example: H₂ in STO-3G

**Setup:** H₂ has 2 electrons in a minimal basis of 2 spatial orbitals ($\sigma_g$, $\sigma_u$) giving 4 spin-orbitals: $\sigma_{g\alpha}$, $\sigma_{g\beta}$, $\sigma_{u\alpha}$, $\sigma_{u\beta}$. This maps to 4 qubits under JW.

**Electronic integrals (at equilibrium $R = 0.735$ Å):**

Using PySCF to compute the integrals and PennyLane's `qml.qchem` to perform the JW mapping, the qubit Hamiltonian for H₂ (STO-3G) is:

$$\hat{H} = g_0 \hat{I} + g_1 \hat{Z}_0 + g_2 \hat{Z}_1 + g_3 \hat{Z}_2 + g_4 \hat{Z}_3 + g_5 \hat{Z}_0\hat{Z}_1 + g_6 \hat{Z}_0\hat{Z}_2 + g_7 \hat{Z}_0\hat{Z}_3 + g_8 \hat{Z}_1\hat{Z}_2 + g_9 \hat{Z}_1\hat{Z}_3 + g_{10}\hat{Z}_2\hat{Z}_3 + g_{11}\hat{X}_0\hat{X}_1\hat{Y}_2\hat{Y}_3 + g_{12}\hat{X}_0\hat{Y}_1\hat{Y}_2\hat{X}_3 + g_{13}\hat{Y}_0\hat{X}_1\hat{X}_2\hat{Y}_3 + g_{14}\hat{Y}_0\hat{Y}_1\hat{X}_2\hat{X}_3$$

with numerical coefficients (at $R = 0.735$ Å):

| Term | Coefficient (Hartree) |
|:-----|:-----:|
| $\hat{I}$ | $-0.0988$ |
| $\hat{Z}_0$ | $+0.1711$ |
| $\hat{Z}_1$ | $+0.1711$ |
| $\hat{Z}_2$ | $-0.2232$ |
| $\hat{Z}_3$ | $-0.2232$ |
| $\hat{Z}_0\hat{Z}_1$ | $+0.1686$ |
| $\hat{Z}_0\hat{Z}_2$ | $+0.1205$ |
| $\hat{Z}_0\hat{Z}_3$ | $+0.1659$ |
| $\hat{Z}_1\hat{Z}_2$ | $+0.1659$ |
| $\hat{Z}_1\hat{Z}_3$ | $+0.1205$ |
| $\hat{Z}_2\hat{Z}_3$ | $+0.1743$ |
| $\hat{X}_0\hat{X}_1\hat{Y}_2\hat{Y}_3$ | $-0.0454$ |
| $\hat{X}_0\hat{Y}_1\hat{Y}_2\hat{X}_3$ | $+0.0454$ |
| $\hat{Y}_0\hat{X}_1\hat{X}_2\hat{Y}_3$ | $+0.0454$ |
| $\hat{Y}_0\hat{Y}_1\hat{X}_2\hat{X}_3$ | $-0.0454$ |

The 15-term qubit Hamiltonian acts on $2^4 = 16$-dimensional Hilbert space. The HF reference state is $|1100\rangle$ (first two spin-orbitals occupied). The FCI ground-state energy at equilibrium is $E_{\text{FCI}} = -1.1373$ Hartree.

With **parity mapping + 2-qubit tapering**, this reduces to a 2-qubit Hamiltonian:

$$\hat{H}_{\text{tapered}} = a_0\hat{I} + a_1\hat{Z}_0 + a_2\hat{Z}_1 + a_3\hat{Z}_0\hat{Z}_1 + a_4\hat{X}_0\hat{X}_1$$

with only 5 terms on 2 qubits — dramatically simpler to simulate.

#### 2.4 Qubit Count and Circuit Depth Comparison

| Molecule | Basis | Spatial Orbs | Spin-Orbs (M) | JW Qubits | Parity (tapered) | UCCSD Params | CNOT Depth (UCCSD) |
|:---------|:------|:----------:|:-------:|:-------:|:---------:|:--------:|:---------:|
| H₂ | STO-3G | 2 | 4 | 4 | 2 | 3 | ~16 |
| LiH | STO-3G | 6 | 12 | 12 | 10 | ~52 | ~1200 |
| H₂O | STO-3G | 7 | 14 | 14 | 12 | ~80 | ~3000 |
| H₂O | 6-31G* | 19 | 38 | 38 | 36 | ~2000 | ~80000 |

The CNOT depth scales roughly as $O(M^3 N_e)$ for UCCSD, making large molecules intractable on NISQ hardware without active space reduction.

---

### 3. Variational Quantum Eigensolver (VQE)

#### 3.1 The Variational Principle

The variational principle states that for any normalised trial state $|\psi(\boldsymbol{\theta})\rangle$:

$$E(\boldsymbol{\theta}) = \langle\psi(\boldsymbol{\theta})|\hat{H}|\psi(\boldsymbol{\theta})\rangle \geq E_0$$

where $E_0$ is the exact ground-state energy. This follows directly from expanding $|\psi(\boldsymbol{\theta})\rangle$ in the eigenbasis $\{|E_k\rangle\}$ of $\hat{H}$:

$$E(\boldsymbol{\theta}) = \sum_k |c_k|^2 E_k \geq E_0 \sum_k |c_k|^2 = E_0$$

since $E_k \geq E_0$ for all $k$ and $\sum_k|c_k|^2 = 1$.

This means any parameterised quantum circuit that prepares $|\psi(\boldsymbol{\theta})\rangle$ gives an **upper bound** on the ground-state energy. By minimising $E(\boldsymbol{\theta})$ over $\boldsymbol{\theta}$, we approach $E_0$ from above. The tightness of the bound depends on the expressibility of the ansatz — whether the parameterised family can represent (or closely approximate) the true ground state.

#### 3.2 The UCCSD Ansatz

**Classical Coupled Cluster (CC):** The coupled-cluster wavefunction is:

$$|\Psi_{\text{CC}}\rangle = e^{\hat{T}}|\Phi_{\text{HF}}\rangle$$

where $\hat{T} = \hat{T}_1 + \hat{T}_2 + \cdots$ is the cluster operator. For CCSD (Singles and Doubles):

$$\hat{T}_1 = \sum_{i \in \text{occ}} \sum_{a \in \text{virt}} t_i^a\, \hat{a}_a^\dagger \hat{a}_i$$

$$\hat{T}_2 = \sum_{i<j \in \text{occ}} \sum_{a<b \in \text{virt}} t_{ij}^{ab}\, \hat{a}_a^\dagger \hat{a}_b^\dagger \hat{a}_j \hat{a}_i$$

where $i, j$ label occupied orbitals, $a, b$ label virtual (unoccupied) orbitals, and $t_i^a$, $t_{ij}^{ab}$ are the cluster amplitudes (variational parameters).

**Problem for quantum computing:** The classical CC operator $e^{\hat{T}}$ is not unitary ($\hat{T}$ is not anti-Hermitian), so it cannot be directly implemented as a quantum circuit.

**Unitary Coupled Cluster (UCC):** Replace $\hat{T}$ with the anti-Hermitian operator $\hat{T} - \hat{T}^\dagger$:

$$|\Psi_{\text{UCC}}\rangle = e^{\hat{T} - \hat{T}^\dagger}|\Phi_{\text{HF}}\rangle$$

This is manifestly unitary: $\hat{U} = e^{\hat{T} - \hat{T}^\dagger}$ satisfies $\hat{U}^\dagger\hat{U} = \hat{I}$ since $(\hat{T} - \hat{T}^\dagger)^\dagger = -(\hat{T} - \hat{T}^\dagger)$.

For UCCSD:

$$\hat{T} - \hat{T}^\dagger = \sum_{i,a} t_i^a(\hat{a}_a^\dagger\hat{a}_i - \hat{a}_i^\dagger\hat{a}_a) + \sum_{i<j,a<b} t_{ij}^{ab}(\hat{a}_a^\dagger\hat{a}_b^\dagger\hat{a}_j\hat{a}_i - \hat{a}_i^\dagger\hat{a}_j^\dagger\hat{a}_b\hat{a}_a)$$

The number of parameters is $n_{\text{singles}} = N_{\text{occ}} \times N_{\text{virt}}$ plus $n_{\text{doubles}} = \binom{N_{\text{occ}}}{2}\binom{N_{\text{virt}}}{2}$.

#### 3.3 Trotterisation of UCCSD into Gates

The UCC unitary $e^{\hat{T}-\hat{T}^\dagger}$ cannot be decomposed exactly into single/two-qubit gates because the individual excitation operators do not commute. We use the first-order Suzuki-Trotter decomposition:

$$e^{\hat{T}-\hat{T}^\dagger} \approx \prod_{\mu} e^{\theta_\mu \hat{\tau}_\mu}$$

where $\hat{\tau}_\mu$ are individual excitation generators (single or double), $\theta_\mu$ are the corresponding amplitudes, and the product is over all excitations.

**Single excitation gate:** For the excitation $i \to a$:

$$e^{\theta(a_a^\dagger a_i - a_i^\dagger a_a)}$$

Under JW, this maps to a rotation involving $X$, $Y$ Pauli strings with $Z$-strings between orbital indices $i$ and $a$. PennyLane provides `qml.SingleExcitation(theta, wires=[i, a])` which implements this as a Givens rotation:

$$\text{SingleExcitation}(\theta) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos(\theta/2) & -\sin(\theta/2) & 0 \\ 0 & \sin(\theta/2) & \cos(\theta/2) & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

acting on the $\{|01\rangle, |10\rangle\}$ subspace of wires $[i, a]$.

**Double excitation gate:** For $i, j \to a, b$:

$$e^{\theta(a_a^\dagger a_b^\dagger a_j a_i - a_i^\dagger a_j^\dagger a_b a_a)}$$

PennyLane provides `qml.DoubleExcitation(theta, wires=[i, j, a, b])` which implements a 4-qubit rotation in the $\{|1100\rangle, |0011\rangle\}$ subspace. This gate decomposes into $\sim 16$ CNOT gates.

**Trotterisation error:** The first-order Trotter formula has error $O(\sum_{\mu < \nu} \|[\hat{\tau}_\mu, \hat{\tau}_\nu]\| \cdot |\theta_\mu||\theta_\nu|)$. In practice, for VQE, the Trotter error is absorbed into the variational parameters — the optimizer compensates for the discretisation, so the ground-state energy can still be reached if the Trotterised ansatz is expressive enough.

#### 3.4 Hardware-Efficient Ansatz

An alternative to UCCSD is the **hardware-efficient ansatz (HEA)**, which consists of layers of native single-qubit rotations and entangling gates that are directly executable on a given hardware topology:

$$|\psi(\boldsymbol{\theta})\rangle = \prod_{l=1}^{L}\left[\hat{U}_{\text{ent}} \cdot \prod_{i=1}^{N_q} R_Y(\theta_i^{(l)}) R_Z(\phi_i^{(l)})\right] |\Phi_{\text{HF}}\rangle$$

where $\hat{U}_{\text{ent}}$ is a layer of CNOT (or CZ) gates following the hardware connectivity.

**Advantages over UCCSD:**
- Much shallower circuits (fewer CNOTs per layer).
- No chemistry-specific structure needed — can be used for any Hamiltonian.
- Circuit depth is a tunable hyperparameter (number of layers $L$).

**Disadvantages:**
- No physical motivation: the circuit does not encode any knowledge of the electronic structure problem.
- Suffers from barren plateaus more readily as $N_q$ and $L$ increase.
- Does not preserve electron number or spin symmetry unless explicitly constrained.
- Convergence to the ground state is not guaranteed (the landscape may have many local minima).

For small molecules (H₂, LiH), HEA with 1–3 layers can reach chemical accuracy. For larger molecules, UCCSD (or adaptive variants like ADAPT-VQE) is generally superior.

#### 3.5 The Full VQE Loop

1. **State preparation:** Prepare $|\psi(\boldsymbol{\theta})\rangle = \hat{U}(\boldsymbol{\theta})|\Phi_{\text{HF}}\rangle$ on the quantum computer.

2. **Measurement:** Measure $\langle\psi(\boldsymbol{\theta})|\hat{H}|\psi(\boldsymbol{\theta})\rangle$ by decomposing $\hat{H}$ into Pauli terms:
$$\hat{H} = \sum_k c_k \hat{P}_k$$
Measure each $\langle\hat{P}_k\rangle$ from repeated preparations + measurements. Group commuting Pauli terms to reduce the number of distinct measurement circuits.

3. **Classical optimization:** Update $\boldsymbol{\theta}$ using a classical optimizer (gradient-based or gradient-free) to minimise $E(\boldsymbol{\theta})$.

4. **Convergence:** Repeat until $|E(\boldsymbol{\theta}^{(k)}) - E(\boldsymbol{\theta}^{(k-1)})| < \varepsilon$ or maximum iterations reached.

The hybrid quantum-classical nature of VQE keeps quantum circuit depth shallow (only state preparation + measurement) while offloading the optimisation to classical hardware. The bottleneck is the number of circuit evaluations needed for sufficient shot statistics and optimiser convergence.

---

### 4. Potential Energy Surface (PES)

#### 4.1 Definition

The **potential energy surface** (PES) is the electronic energy $E_0(\mathbf{R})$ as a function of nuclear coordinates $\mathbf{R}$:

$$\text{PES}: \mathbb{R}^{3N_n - 6} \to \mathbb{R}$$

(subtracting 3 translational and 3 rotational degrees of freedom, or $3N_n - 5$ for linear molecules).

For a diatomic molecule like H₂, the PES is a one-dimensional curve $E(R)$ where $R$ is the internuclear distance. For polyatomics, the PES is a high-dimensional hypersurface.

#### 4.2 What the PES Reveals

- **Equilibrium geometry $R_e$:** The minimum of $E(R)$, where $dE/dR = 0$ and $d^2E/dR^2 > 0$. For H₂, $R_e \approx 0.74$ Å.

- **Dissociation energy $D_e$:** The energy difference between the minimum and the asymptotic separated-atom limit: $D_e = E(R \to \infty) - E(R_e)$. For H₂, $D_e \approx 4.75$ eV (exact: $4.7467$ eV).

- **Bond stretching frequency $\omega_e$:** Related to the curvature at the minimum: $\omega_e = \sqrt{k/\mu}$ where $k = d^2E/dR^2|_{R_e}$ and $\mu$ is the reduced mass.

- **Transition states:** Saddle points on the PES (for polyatomics) correspond to transition states connecting reactants and products.

#### 4.3 Accuracy Hierarchy

The hierarchy of computational chemistry methods, in increasing accuracy and cost:

| Method | Scaling | Correlation | PES Quality |
|:-------|:-------:|:----------:|:----------:|
| HF | $O(M^4)$ | None | Fails at dissociation |
| MP2 | $O(M^5)$ | Perturbative | Fair near equilibrium, fails at dissociation |
| CISD | $O(M^6)$ | Variational, not size-consistent | Good near equilibrium |
| CCSD(T) | $O(M^7)$ | "Gold standard" of single-ref methods | Excellent near equilibrium |
| FCI | Exponential | Exact (within basis) | Exact reference |
| VQE (UCCSD) | Quantum poly | Variational, size-consistent | Approaches FCI |

**Chemical accuracy** is defined as $\Delta E \leq 1.6$ mHartree $\approx 1$ kcal/mol $\approx 43$ meV. Achieving chemical accuracy across the entire PES (including dissociation) requires a method that correctly handles strong correlation — this is where VQE with a multi-reference-capable ansatz can outperform CCSD(T).

#### 4.4 H₂ Dissociation Curve: The Benchmark

The H₂ molecule is the standard test case for VQE. The exact (FCI/STO-3G) dissociation curve has:

- $R_e = 0.735$ Å, $E(R_e) = -1.1373$ Hartree
- $E(R \to \infty) = -0.9997$ Hartree (two separated hydrogen atoms)
- $D_e = 0.1376$ Hartree $= 3.745$ eV (in STO-3G; the basis set limit value is higher)

**RHF failure:** At $R = 3.0$ Å, the RHF energy is $\sim -0.90$ Hartree while the exact energy is $\sim -1.00$ Hartree — an error of $\sim 0.1$ Hartree ($\sim 63$ kcal/mol), completely unacceptable.

**VQE target:** Reproduce the FCI/STO-3G curve to within **chemical accuracy** (1.6 mHartree) at all bond lengths $R \in [0.3, 3.0]$ Å. With the UCCSD ansatz on 4 qubits (or 2 qubits with parity + tapering), this is achievable in noiseless simulation and is the primary validation target.

---

### 5. Molecular Orbitals and Wavefunction Visualization

#### 5.1 LCAO Expansion

Each molecular orbital is expanded in the AO basis:

$$\phi_p(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu p}\, \chi_\mu(\mathbf{r})$$

where $\chi_\mu$ are Gaussian-type basis functions centred on nuclei:

$$\chi_\mu(\mathbf{r}) = N_\mu\, (x - A_x)^{l_x}(y - A_y)^{l_y}(z - A_z)^{l_z}\, e^{-\alpha_\mu |\mathbf{r} - \mathbf{A}|^2}$$

with angular momentum quantum numbers $(l_x, l_y, l_z)$, exponent $\alpha_\mu$, centre $\mathbf{A}$, and normalisation $N_\mu$. The MO coefficients $C_{\mu p}$ come from the converged HF calculation.

#### 5.2 Evaluating $\psi(\mathbf{r})$ on a 3D Grid

To visualise an MO, evaluate $\phi_p(\mathbf{r})$ on a uniform 3D Cartesian grid:

1. Define the grid: $\mathbf{r}_{ijk} = (x_{\min} + i\Delta x, \, y_{\min} + j\Delta y, \, z_{\min} + k\Delta z)$ for $i = 0, \ldots, N_x-1$, etc. Typical grid: $\Delta x = \Delta y = \Delta z = 0.1$ Bohr, extent $\pm 5$ Bohr from molecular centre, giving $\sim 100^3 = 10^6$ grid points.

2. For each grid point, compute each AO value:
$$\chi_\mu(\mathbf{r}_{ijk}) = N_\mu \prod_d (r_d - A_d)^{l_d} \exp(-\alpha_\mu|\mathbf{r}_{ijk} - \mathbf{A}|^2)$$

3. Contract with MO coefficients:
$$\phi_p(\mathbf{r}_{ijk}) = \sum_\mu C_{\mu p}\, \chi_\mu(\mathbf{r}_{ijk})$$

This produces a **volumetric scalar field** $V_p[i,j,k] = \phi_p(\mathbf{r}_{ijk})$ of shape $(N_x, N_y, N_z)$.

**Optimisation:** For GTOs, the exponential decay makes most basis functions negligible beyond a few Bohr from their centre. Use a cutoff radius per AO (e.g., where $\chi_\mu < 10^{-8}$) and only evaluate within that radius. The evaluation is vectorised by forming a displacement tensor $\Delta_{g\mu d} = r_{gd} - A_{\mu d}$ of shape $(N_{\text{grid}}, K, 3)$, computing squared distances $|\Delta_{g\mu}|^2 = \sum_d \Delta_{g\mu d}^2$, and evaluating all Gaussians in a single broadcasted exponential $G_{g\mu} = \exp(-\alpha_\mu |\Delta_{g\mu}|^2)$. The angular polynomial prefactors are applied element-wise, and the final MO values are obtained by contracting with MO coefficients: $\phi_p(\mathbf{r}_g) = \sum_\mu C_{\mu p}\, [\text{angular}]_{g\mu}\, G_{g\mu}$.

*(Full implementation: `report2_scripts/orbital_grid.py`)*

#### 5.3 Electron Density and the One-Particle Reduced Density Matrix

The **electron density** is the sum over occupied MOs:

$$\rho(\mathbf{r}) = \sum_{i \in \text{occ}} n_i\, |\phi_i(\mathbf{r})|^2$$

where $n_i$ is the occupation number (2 for RHF doubly-occupied orbitals). More generally, for a correlated wavefunction with one-particle reduced density matrix (1-RDM) $\gamma_{pq}$:

$$\rho(\mathbf{r}) = \sum_{p,q} \gamma_{pq}\, \phi_p^*(\mathbf{r})\phi_q(\mathbf{r})$$

where $\gamma_{pq} = \langle\Psi|\hat{a}_p^\dagger\hat{a}_q|\Psi\rangle$. For VQE, the 1-RDM can be measured on the quantum computer by measuring $\langle\hat{a}_p^\dagger\hat{a}_q\rangle$ in the Pauli basis (each $\hat{a}_p^\dagger\hat{a}_q$ maps to a sum of Pauli strings under JW).

#### 5.4 Isosurface Extraction

An **isosurface** is the set of points where a scalar field equals a constant value:

$$S_c = \{\mathbf{r} : \phi_p(\mathbf{r}) = c\}$$

For MO visualisation, the conventional isovalues are $c = \pm 0.05$ a.u. (atomic units). The positive lobe ($c = +0.05$) shows where the wavefunction is positive; the negative lobe ($c = -0.05$) shows where it's negative. These lobes are typically coloured blue/red respectively.

**Why $\pm 0.05$ a.u.?** This isovalue encloses approximately 90% of the electron density associated with the orbital, providing a visually useful representation of the orbital shape without being too diffuse (smaller $|c|$) or too compact (larger $|c|$).

The **marching cubes algorithm** (Lorensen & Cline, 1987) extracts a triangle mesh approximating the isosurface from the volumetric data. For each voxel (cube formed by 8 neighbouring grid points):

1. Classify each corner as inside ($\phi > c$) or outside ($\phi < c$) the isosurface.
2. The 8 binary corner classifications give $2^8 = 256$ possible configurations, reduced to 15 unique topological cases by symmetry.
3. Look up the corresponding triangle configuration in a precomputed table (the **marching cubes lookup table**).
4. Interpolate vertex positions along cube edges where the isosurface crosses (linear interpolation between corner values).
5. Compute face normals from the gradient of the scalar field at each vertex (or from the triangle vertices directly).

---

## Part 2 — Software Architecture

### Project Structure

```
vqe_explorer/
├── main.py                  # Orchestration, CLI, entry point
├── molecule.py              # Molecule definition, geometry, basis set
├── classical_hf.py          # PySCF interface: HF, integrals, MO coefficients
├── hamiltonian_builder.py   # Second-quantized → qubit Hamiltonian
├── ansatz.py                # UCCSD and hardware-efficient ansätze
├── vqe_solver.py            # VQE optimization loop
├── pes_scanner.py           # PES sweep over bond lengths
├── orbital_grid.py          # MO evaluation on 3D grid
├── isosurface.py            # Marching cubes isosurface extraction
├── renderer.py              # OpenGL 3D renderer
├── hardware_runner.py       # IBM / IonQ hardware submission
├── benchmarker.py           # Method comparison (HF, CISD, CCSD(T), FCI, VQE)
├── requirements.txt
├── config.yaml              # Default configuration
└── tests/
    ├── test_molecule.py
    ├── test_hamiltonian.py
    ├── test_vqe.py
    ├── test_pes.py
    ├── test_orbital.py
    ├── test_isosurface.py
    └── test_renderer.py
```

### Data Flow Diagram

```
                        ┌──────────────┐
                        │ config.yaml  │
                        └──────┬───────┘
                               │
                        ┌──────▼───────┐
                        │   main.py    │   CLI / orchestration
                        └──────┬───────┘
                               │
              ┌────────────────┼────────────────────┐
              │                │                     │
     ┌────────▼───────┐ ┌─────▼──────────┐ ┌───────▼──────────┐
     │  molecule.py   │ │ pes_scanner.py │ │  benchmarker.py  │
     │ (geometry,     │ │ (bond length   │ │ (HF/CISD/CCSD(T) │
     │  basis set)    │ │  loop)         │ │  /FCI comparison) │
     └────────┬───────┘ └─────┬──────────┘ └───────┬──────────┘
              │               │                     │
     ┌────────▼───────┐      │                     │
     │ classical_hf.py│◄─────┘─────────────────────┘
     │ (PySCF: HF,    │
     │  integrals, MO)│
     └────────┬───────┘
              │
     ┌────────▼──────────────┐
     │ hamiltonian_builder.py│
     │ (JW/Parity mapping,  │
     │  QuTiP cross-check)  │
     └────────┬──────────────┘
              │
    ┌─────────┼──────────────────┐
    │         │                  │
┌───▼───┐ ┌──▼──────────┐ ┌─────▼─────────────┐
│ansatz │ │ vqe_solver  │ │ hardware_runner.py│
│ .py   │ │   .py       │ │ (IBM/IonQ device) │
└───┬───┘ └──┬──────────┘ └─────┬─────────────┘
    │        │                   │
    └────────┼───────────────────┘
             │
    ┌────────▼───────┐     ┌───────────────┐
    │ orbital_grid.py│────►│isosurface.py  │
    │ (3D MO grid)   │     │(marching cubes)│
    └────────────────┘     └───────┬───────┘
                                   │
                           ┌───────▼───────┐
                           │  renderer.py  │
                           │ (OpenGL:      │
                           │  isosurface + │
                           │  PES curve +  │
                           │  mol. geom.)  │
                           └───────────────┘
```

### Module Specifications

#### 1. `molecule.py`

**Purpose:** Define molecular systems with geometry, basis set, and charge/multiplicity.

This module provides two core data classes and a set of factory/utility functions:

- **`Atom`** — Holds an element symbol (e.g., `'H'`, `'O'`) and a Cartesian position vector $(x, y, z)$ in Ångströms.
- **`MoleculeSpec`** — Aggregates a list of `Atom` objects together with basis set name (default `sto-3g`), net charge (default 0), spin multiplicity $2S+1$ (default 1 = singlet), and a human-readable name for logging.

| Function | Description | Key Parameters | Returns |
|:---------|:------------|:---------------|:--------|
| `create_h2(bond_length)` | Places two H atoms symmetrically along the $z$-axis at $\pm R/2$. | $R$ in Å (default 0.735) | `MoleculeSpec` for H₂ |
| `create_lih(bond_length)` | Places Li at the origin and H at distance $R$ along $z$. | $R$ in Å (default 1.595) | `MoleculeSpec` for LiH |
| `create_h2o(oh_length, angle_deg)` | Constructs a bent H₂O geometry with O at the origin and two H atoms placed at bond length $R_{\text{OH}}$ and opening angle $\alpha$. | $R_{\text{OH}}$ in Å (default 0.957), $\alpha$ in degrees (default 104.5) | `MoleculeSpec` for H₂O |
| `set_bond_length(mol, a, b, new_length)` | Returns a deep copy of `mol` with the distance between atoms `a` and `b` set to `new_length`, keeping the bond midpoint fixed. The direction vector $\hat{\mathbf{d}} = (\mathbf{r}_b - \mathbf{r}_a)/\|\mathbf{r}_b - \mathbf{r}_a\|$ is preserved while both atoms are displaced symmetrically. | Atom indices, new $R$ in Å | New `MoleculeSpec` |

```python
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np


@dataclass
class Atom:
    """A single atom with element symbol and Cartesian position (Å)."""
    symbol: str
    position: np.ndarray  # shape: (3,)


@dataclass
class MoleculeSpec:
    """Full specification of a molecular system."""
    atoms: list[Atom]
    basis: str = "sto-3g"
    charge: int = 0
    multiplicity: int = 1
    name: str = ""


def create_h2(bond_length: float = 0.735) -> MoleculeSpec:
    """Place two H atoms symmetrically along z at ±R/2."""
    ...


def create_lih(bond_length: float = 1.595) -> MoleculeSpec:
    """Place Li at origin and H at distance R along z."""
    ...


def create_h2o(
    oh_length: float = 0.957,
    angle_deg: float = 104.5,
) -> MoleculeSpec:
    """Construct bent H₂O: O at origin, two H at R_OH and angle α."""
    ...


def set_bond_length(
    mol: MoleculeSpec,
    atom_idx_a: int,
    atom_idx_b: int,
    new_length: float,
) -> MoleculeSpec:
    """Return a deep copy with bond a–b set to new_length (Å).

    Preserves bond midpoint and direction unit vector.
    """
    ...
```

*(Full implementation: `report2_scripts/molecule.py`)*

#### 2. `classical_hf.py`

**Purpose:** Interface to PySCF for Hartree-Fock calculations, integral extraction, and MO coefficients.

This module wraps PySCF to provide all the classical quantum-chemistry quantities needed downstream. It exposes a result container and four functions:

- **`HFResult`** — Data class bundling everything produced by a converged RHF calculation:
  - Total HF energy $E_{\text{HF}}$ and nuclear repulsion energy $E_{\text{nuc}}$ (both in Hartree).
  - MO coefficient matrix $\mathbf{C}$ of shape $(K, M)$ where $K$ = number of AO basis functions and $M$ = number of MOs.
  - Orbital energies $\{\epsilon_p\}$ of shape $(M,)$.
  - One-electron integrals $h_{pq}$ in the MO basis, shape $(M, M)$.
  - Two-electron integrals $h_{pqrs}$ in the MO basis, shape $(M, M, M, M)$.
  - AO overlap matrix $\mathbf{S}$ of shape $(K, K)$, electron count $N_e$, orbital count $M$, and AO labels.

| Function | What it computes | Returns |
|:---------|:-----------------|:--------|
| `run_hartree_fock(mol_spec)` | Builds a PySCF `Mole` object from the `MoleculeSpec`, runs RHF SCF, then extracts the converged Fock matrix $\mathbf{F}$, MO coefficients $\mathbf{C}$, orbital energies $\{\epsilon_p\}$, and transforms the one- and two-electron integrals to the MO basis via $h_{pq} = \mathbf{C}^T \mathbf{h}^{\text{core}} \mathbf{C}$ and the four-index AO-to-MO transformation for $h_{pqrs}$. | `HFResult` |
| `get_fci_energy(mol_spec)` | Runs RHF followed by PySCF's Full CI solver to obtain the exact ground-state energy within the given basis. | $E_{\text{FCI}}$ (float, Hartree) |
| `get_cisd_energy(mol_spec)` | Runs RHF + Configuration Interaction Singles and Doubles. | $E_{\text{CISD}}$ (float, Hartree) |
| `get_ccsd_t_energy(mol_spec)` | Runs RHF + CCSD, then adds the perturbative triples correction (T). | $E_{\text{CCSD(T)}}$ (float, Hartree) |

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import ao2mo, gto, scf, fci, ci, cc


@dataclass
class HFResult:
    """Results from a converged RHF calculation."""
    energy_hf: float                        # E_HF (Hartree)
    energy_nuc: float                       # nuclear repulsion (Hartree)
    mo_coefficients: np.ndarray             # shape: (K, M) — AO-to-MO matrix
    mo_energies: np.ndarray                 # shape: (M,)
    one_electron_integrals: np.ndarray      # h_pq in MO basis, shape: (M, M)
    two_electron_integrals: np.ndarray      # h_pqrs in MO basis, shape: (M, M, M, M)
    overlap_matrix: np.ndarray              # S_μν in AO basis, shape: (K, K)
    n_electrons: int
    n_orbitals: int
    ao_labels: list[str]


def run_hartree_fock(mol_spec: MoleculeSpec) -> HFResult:
    """Build PySCF Mole, run RHF, extract MO-basis integrals.

    Transforms one-electron integrals via C^T h^core C and
    two-electron integrals via ao2mo.full.
    """
    ...


def get_fci_energy(mol_spec: MoleculeSpec) -> float:
    """Run RHF + Full CI; return E_FCI (Hartree)."""
    ...


def get_cisd_energy(mol_spec: MoleculeSpec) -> float:
    """Run RHF + CISD; return E_CISD (Hartree)."""
    ...


def get_ccsd_t_energy(mol_spec: MoleculeSpec) -> float:
    """Run RHF + CCSD(T); return E_CCSD(T) (Hartree)."""
    ...
```

*(Full implementation: `report2_scripts/classical_hf.py`)*

#### 3. `hamiltonian_builder.py`

**Purpose:** Build the qubit Hamiltonian from molecular integrals using Jordan-Wigner or Parity mapping. Cross-validate with QuTiP exact diagonalisation.

This module converts classical molecular integrals into a qubit-operator representation suitable for VQE, and provides an independent validation path.

- **`QubitHamiltonian`** — Data class holding the PennyLane Hamiltonian object, qubit count $N_q$, number of Pauli terms, the mapping name (`'jordan_wigner'` or `'parity'`), and an optional exact ground-state energy for validation.

| Function | Description |
|:---------|:------------|
| `build_hamiltonian(h_pq, h_pqrs, n_electrons, E_nuc, mapping, ...)` | Assembles the second-quantised fermionic Hamiltonian $\hat{H} = \sum_{pq}h_{pq}\hat{a}_p^\dagger\hat{a}_q + \tfrac{1}{2}\sum_{pqrs}h_{pqrs}\hat{a}_p^\dagger\hat{a}_q^\dagger\hat{a}_s\hat{a}_r + E_{\text{nuc}}$ from the supplied MO-basis integrals, then applies the chosen fermion-to-qubit mapping (Jordan-Wigner or parity) via PennyLane. Optionally restricts to an active space of the specified number of active electrons and orbitals. Returns a `QubitHamiltonian`. |
| `validate_with_qutip(qubit_hamiltonian)` | Converts the PennyLane Hamiltonian to a QuTiP `Qobj`, performs exact diagonalisation, and returns the ground-state energy $E_0$. This provides an independent check that the fermion-to-qubit mapping is correct. |
| `pennylane_to_qutip(hamiltonian, n_qubits)` | Maps each Pauli string $c_k \hat{P}_k$ in the PennyLane Hamiltonian to a tensor product of QuTiP Pauli matrices ($\hat{\sigma}_x$, $\hat{\sigma}_y$, $\hat{\sigma}_z$, $\hat{I}$), accumulating into a single QuTiP `Qobj` of dimension $(2^{N_q}, 2^{N_q})$. |

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pennylane as qml
import qutip


@dataclass
class QubitHamiltonian:
    """Container for a qubit Hamiltonian and metadata."""
    hamiltonian: qml.Hamiltonian
    n_qubits: int
    n_paulis: int
    mapping: str
    fci_energy: float | None = None


def build_hamiltonian(
    one_electron_integrals: np.ndarray,
    two_electron_integrals: np.ndarray,
    n_electrons: int,
    nuclear_repulsion: float,
    mapping: str = "jordan_wigner",
    active_electrons: int | None = None,
    active_orbitals: int | None = None,
) -> QubitHamiltonian:
    """Assemble fermionic Hamiltonian from MO integrals and map to qubits.

    Supports Jordan-Wigner and parity mappings via PennyLane.
    Optionally restricts to an active space.
    """
    ...


def validate_with_qutip(qubit_hamiltonian: QubitHamiltonian) -> float:
    """Convert to QuTiP Qobj, diagonalise, return ground-state energy E₀."""
    ...


def pennylane_to_qutip(
    hamiltonian: qml.Hamiltonian,
    n_qubits: int,
) -> qutip.Qobj:
    """Map PennyLane Pauli Hamiltonian to a QuTiP Qobj of dim (2^N, 2^N)."""
    ...
```

*(Full implementation: `report2_scripts/hamiltonian_builder.py`)*

#### 4. `ansatz.py`

**Purpose:** Construct UCCSD and hardware-efficient ansätze using PennyLane.

This module provides two ansatz constructors, each returning a callable circuit function and a parameter count.

> **`build_uccsd_ansatz(n_qubits, n_electrons, hf_state=None)`** constructs the Unitary Coupled Cluster Singles and Doubles ansatz. It first enumerates all spin-allowed single excitations $\{(i \to a)\}$ and double excitations $\{(i,j \to a,b)\}$ using PennyLane's `qchem.excitations`. The UCCSD unitary is:
> $$|\Psi(\boldsymbol{\theta})\rangle = e^{\hat{T}(\boldsymbol{\theta}) - \hat{T}^\dagger(\boldsymbol{\theta})}|\Phi_{\text{HF}}\rangle$$
> where $\hat{T} = \hat{T}_1 + \hat{T}_2$ with $\hat{T}_1 = \sum_{i,a}\theta_i^a \hat{a}_a^\dagger \hat{a}_i$ and $\hat{T}_2 = \sum_{i<j,a<b}\theta_{ij}^{ab}\hat{a}_a^\dagger\hat{a}_b^\dagger\hat{a}_j\hat{a}_i$. PennyLane's `qml.UCCSD` (or individual `SingleExcitation` / `DoubleExcitation` gates) compiles this into a gate sequence via first-order Trotterisation. The HF reference state fills the first $N_e$ spin-orbitals (or uses a user-supplied binary vector). Returns `(ansatz_fn, n_params, excitation_list)`.

> **`build_hardware_efficient_ansatz(n_qubits, n_layers=2, entangler='CNOT')`** builds a hardware-efficient ansatz consisting of $L$ layers, each comprising $R_Y(\theta)$ and $R_Z(\phi)$ rotations on every qubit followed by a ladder of CNOT (or CZ) entangling gates on nearest-neighbour pairs. The total parameter count is $2 N_q L$. Returns `(ansatz_fn, n_params)`.

```python
from __future__ import annotations

from typing import Callable

import numpy as np
import pennylane as qml


def build_uccsd_ansatz(
    n_qubits: int,
    n_electrons: int,
    hf_state: np.ndarray | None = None,
) -> tuple[Callable, int, list]:
    """Construct UCCSD ansatz from PennyLane qchem excitations.

    Returns (ansatz_fn, n_params, excitation_list).
    ansatz_fn accepts a parameter array of shape (n_params,).
    """
    ...


def build_hardware_efficient_ansatz(
    n_qubits: int,
    n_layers: int = 2,
    entangler: str = "CNOT",
) -> tuple[Callable, int]:
    """Build L-layer hardware-efficient ansatz (RY, RZ + entangler).

    Total parameters: 2 * n_qubits * n_layers.
    Returns (ansatz_fn, n_params).
    """
    ...
```

*(Full implementation: `report2_scripts/ansatz.py`)*

#### 5. `vqe_solver.py`

**Purpose:** VQE optimisation loop with multiple optimiser backends.

- **`VQEResult`** — Data class storing the best energy $E^*$ (Hartree), optimal parameters $\boldsymbol{\theta}^*$, full trajectory $(\boldsymbol{\theta}^{(k)}, E^{(k)})$ at each iteration, total iteration count, convergence flag, and gradient norms.

> **`run_vqe(hamiltonian, ansatz_fn, n_params, n_qubits, ...)`** executes the VQE optimisation loop. A PennyLane QNode wraps the ansatz and evaluates $E(\boldsymbol{\theta}) = \langle\psi(\boldsymbol{\theta})|\hat{H}|\psi(\boldsymbol{\theta})\rangle$. The solver supports four optimiser backends:
>
> - **Adam / GradientDescent**: PennyLane built-in optimisers using automatic differentiation (backprop or parameter-shift).
> - **COBYLA**: SciPy's derivative-free constrained optimiser.
> - **SPSA** (Simultaneous Perturbation Stochastic Approximation): estimates the gradient from two function evaluations per step via random perturbation $\hat{g}_k = \frac{f(\boldsymbol{\theta}+c_k\boldsymbol{\Delta}) - f(\boldsymbol{\theta}-c_k\boldsymbol{\Delta})}{2c_k\boldsymbol{\Delta}}$ with decaying step sizes $a_k = a_0/(k+A)^\alpha$ and $c_k = c_0/k^\gamma$.
>
> Convergence is declared when $|E^{(k)} - E^{(k-1)}| < \varepsilon$. Initial parameters default to zero with a small random perturbation to break symmetry. An optional callback is invoked each iteration with `(iteration, params, energy)`.

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pennylane as qml


@dataclass
class VQEResult:
    """Container for VQE optimisation results."""
    optimal_energy: float
    optimal_params: np.ndarray
    trajectory: list[tuple[np.ndarray, float]]
    n_iterations: int
    converged: bool
    gradient_norms: list[float]


def run_vqe(
    hamiltonian: qml.Hamiltonian,
    ansatz_fn: Callable,
    n_params: int,
    n_qubits: int,
    init_params: np.ndarray | None = None,
    optimizer: str = "Adam",
    learning_rate: float = 0.05,
    maxiter: int = 200,
    tol: float = 1e-6,
    device_name: str = "default.qubit",
    diff_method: str = "backprop",
    callback: Callable | None = None,
) -> VQEResult:
    """Run the VQE optimisation loop.

    Supports Adam, GradientDescent, COBYLA, and SPSA backends.
    Convergence declared when |ΔE| < tol.
    """
    ...
```

*(Full implementation: `report2_scripts/vqe_solver.py`)*

#### 6. `pes_scanner.py`

**Purpose:** Sweep over bond lengths, run VQE at each geometry, collect PES data.

- **`PESData`** — Data class holding arrays of bond lengths (Å), VQE energies, HF energies, FCI energies, and optionally CISD and CCSD(T) energies, plus the list of optimal VQE parameter vectors at each geometry.

> **`scan_pes(molecule_factory, bond_lengths, ...)`** iterates over an array of bond lengths $\{R_1, \ldots, R_N\}$. At each $R_k$ it:
>
> 1. Calls `molecule_factory(R_k)` to build the geometry.
> 2. Runs Hartree-Fock via `run_hartree_fock` to extract integrals.
> 3. Builds the qubit Hamiltonian via `build_hamiltonian`.
> 4. Constructs the chosen ansatz (UCCSD or hardware-efficient).
> 5. Runs VQE, optionally **warm-starting** from the previous geometry's optimal parameters $\boldsymbol{\theta}^*_{k-1}$ to accelerate convergence along the scan.
> 6. Optionally computes HF, CISD, CCSD(T), and FCI energies via `classical_hf` for comparison.
>
> The warm-start strategy exploits the fact that molecular orbitals change smoothly with geometry, so $\boldsymbol{\theta}^*_{k-1}$ is a good initial guess for $\boldsymbol{\theta}^*_k$.

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class PESData:
    """Container for PES scan results."""
    bond_lengths: np.ndarray
    energies_vqe: np.ndarray
    energies_hf: np.ndarray
    energies_fci: np.ndarray
    energies_cisd: np.ndarray | None = None
    energies_ccsd_t: np.ndarray | None = None
    vqe_params: list[np.ndarray] | None = None


def scan_pes(
    molecule_factory: Callable[[float], MoleculeSpec],
    bond_lengths: np.ndarray,
    basis: str = "sto-3g",
    mapping: str = "jordan_wigner",
    ansatz_type: str = "uccsd",
    optimizer: str = "Adam",
    maxiter: int = 200,
    use_prev_params: bool = True,
    compute_classical: bool = True,
    callback: Callable | None = None,
) -> PESData:
    """Sweep bond lengths, run VQE at each geometry, collect PES data.

    Warm-starts from previous geometry's optimal parameters when
    use_prev_params is True and parameter count is unchanged.
    """
    ...
```

*(Full implementation: `report2_scripts/pes_scanner.py`)*

#### 7. `orbital_grid.py`

**Purpose:** Evaluate molecular orbital wavefunctions on a 3D Cartesian grid.

- **`OrbitalVolumeData`** — Data class holding the 3D scalar field array of shape $(N_x, N_y, N_z)$, grid origin, grid spacing (all in Bohr), grid shape, and the MO index.

> **`evaluate_orbital_on_grid(mo_coefficients, orbital_index, atom_coords, atom_symbols, basis_name, grid_extent=5.0, grid_spacing=0.1)`** computes the MO wavefunction $\phi_p(\mathbf{r})$ on a uniform Cartesian grid. The grid extends $\pm$`grid_extent` Bohr from the molecular centre in each direction with spacing `grid_spacing`. At each grid point $\mathbf{r}_g$, PySCF's `eval_gto` evaluates all AO basis functions $\{\chi_\mu(\mathbf{r}_g)\}$ efficiently (exploiting Gaussian sparsity), then contracts with the $p$-th column of the MO coefficient matrix:
> $$\phi_p(\mathbf{r}_g) = \sum_{\mu} C_{\mu p}\, \chi_\mu(\mathbf{r}_g)$$
> The result is reshaped into a 3D volume of shape $(N_x, N_y, N_z)$.

> **`evaluate_electron_density(mo_coefficients, n_occupied, ...)`** computes the total RHF electron density on the same grid type by summing over occupied spatial orbitals:
> $$\rho(\mathbf{r}) = \sum_{i=0}^{N_{\text{occ}}-1} 2\,|\phi_i(\mathbf{r})|^2$$
> The factor of 2 accounts for double occupation in RHF. Returns an `OrbitalVolumeData` with `orbital_index = -1` to indicate a density rather than a single orbital.

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyscf import gto


@dataclass
class OrbitalVolumeData:
    """Volumetric data for a molecular orbital on a 3D grid."""
    values: np.ndarray           # shape: (N_x, N_y, N_z)
    grid_origin: tuple[float, float, float]    # in Bohr
    grid_spacing: tuple[float, float, float]   # in Bohr
    grid_shape: tuple[int, int, int]
    orbital_index: int


def evaluate_orbital_on_grid(
    mo_coefficients: np.ndarray,
    orbital_index: int,
    atom_coords: np.ndarray,
    atom_symbols: list[str],
    basis_name: str,
    grid_extent: float = 5.0,
    grid_spacing: float = 0.1,
) -> OrbitalVolumeData:
    """Evaluate MO φ_p(r) on a uniform 3D Cartesian grid.

    Uses PySCF eval_gto for efficient AO evaluation, then
    contracts with MO coefficient column p.
    """
    ...


def evaluate_electron_density(
    mo_coefficients: np.ndarray,
    n_occupied: int,
    atom_coords: np.ndarray,
    atom_symbols: list[str],
    basis_name: str,
    grid_extent: float = 5.0,
    grid_spacing: float = 0.1,
) -> OrbitalVolumeData:
    """Compute total RHF electron density ρ(r) = Σ 2|φ_i(r)|².

    Returns OrbitalVolumeData with orbital_index = -1.
    """
    ...
```

*(Full implementation: `report2_scripts/orbital_grid.py`)*

#### 8. `isosurface.py`

**Purpose:** Marching cubes isosurface extraction from volumetric data.

- **`IsosurfaceMesh`** — Data class holding a triangle mesh: vertex positions $(N_v, 3)$, vertex normals $(N_v, 3)$, face index triples $(N_f, 3)$, and a per-vertex sign array ($+1$ for positive lobe, $-1$ for negative lobe) used for colouring.

> **`marching_cubes(volume, isovalue, grid_origin, grid_spacing)`** extracts the isosurface $S_c = \{\mathbf{r} : \phi(\mathbf{r}) = c\}$ from a 3D scalar field using scikit-image's implementation of the Lorensen & Cline (1987) algorithm. Each voxel's 8 corners are classified as inside/outside the isosurface, yielding an 8-bit index into a lookup table of 256 triangulation cases (15 unique topologies). Edge-crossing positions are linearly interpolated, and vertex normals are computed from the scalar field gradient. Vertices are translated to world coordinates by adding the grid origin.

> **`extract_dual_isosurface(volume, isovalue, grid_origin, grid_spacing)`** calls `marching_cubes` twice — once at $+c$ and once at $-c$ — then concatenates the two meshes. Vertices from the positive isosurface are tagged with sign $= +1$; negative with sign $= -1$. Face indices of the negative mesh are offset by the vertex count of the positive mesh before concatenation.

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from skimage.measure import marching_cubes as sk_marching_cubes


@dataclass
class IsosurfaceMesh:
    """Triangle mesh from isosurface extraction."""
    vertices: np.ndarray    # shape: (N_v, 3)
    normals: np.ndarray     # shape: (N_v, 3)
    faces: np.ndarray       # shape: (N_f, 3), uint32
    signs: np.ndarray       # shape: (N_v,), +1 or -1


def marching_cubes(
    volume: np.ndarray,
    isovalue: float,
    grid_origin: tuple[float, float, float],
    grid_spacing: tuple[float, float, float],
) -> IsosurfaceMesh:
    """Extract isosurface {r : φ(r) = c} via scikit-image marching cubes.

    Vertices are translated to world coordinates by adding grid_origin.
    """
    ...


def extract_dual_isosurface(
    volume: np.ndarray,
    isovalue: float,
    grid_origin: tuple[float, float, float],
    grid_spacing: tuple[float, float, float],
) -> IsosurfaceMesh:
    """Extract +c and -c isosurfaces; concatenate with sign tags.

    Positive-lobe vertices: sign = +1. Negative: sign = -1.
    Face indices of the negative mesh are offset by |V_pos|.
    """
    ...
```

*(Full implementation: `report2_scripts/isosurface.py`)*

#### 9. `renderer.py`

**Purpose:** OpenGL renderer for orbital isosurfaces, PES curves, and molecular geometry.

This module provides the full 3D rendering pipeline using GLFW + OpenGL 3.3 Core Profile.

| Function / Class | Description |
|:-----------------|:------------|
| `init_opengl_context(width, height, title)` | Initialises a GLFW window with an OpenGL 3.3 core-profile context, enables depth testing, 4× MSAA, and alpha blending. Returns the GLFW window handle. |
| `create_orbital_vao(mesh)` | Uploads an `IsosurfaceMesh` to the GPU. Vertex layout: `[x, y, z, nx, ny, nz, sign]` (7 float32 values, 28 bytes per vertex). Position at attribute location 0, normal at location 1, sign at location 2. Uses `GL_DYNAMIC_DRAW` since the mesh changes with bond length. Returns `(VAO, VBO, EBO, n_indices)`. |
| `update_orbital_vao(vbo, ebo, mesh)` | Orphans the old buffer and uploads new mesh data to the same VBO/EBO — a double-buffering strategy that avoids GPU sync stalls. Returns the new index count. |
| `create_pes_curve_vao(bond_lengths, energies_dict)` | Creates a VAO for rendering multiple PES curves as coloured line strips. Each method (HF, CISD, VQE, FCI, etc.) maps to a distinct RGBA colour. Vertices store `[R, E, 0, r, g, b, a]`. Returns `(VAO, VBO)`. |
| `draw_molecule(atom_positions, atom_symbols, bonds, shader_program)` | Renders atoms as spheres and bonds as cylinders using instanced drawing of unit meshes. Atom colours follow the CPK convention (H = white, O = red, Li = purple, etc.) and radii are scaled van der Waals radii ($\times 0.3$). |
| `ArcballCamera` | Camera class supporting arcball rotation (left-drag), zoom (scroll), and heading/pitch angles. Computes view ($\mathbf{V} = \text{lookAt}$) and perspective ($\mathbf{P}$) matrices for the shader pipeline. |

```python
from __future__ import annotations

import ctypes

import glfw
import numpy as np
from OpenGL.GL import *


def init_opengl_context(
    width: int = 1400,
    height: int = 900,
    title: str = "VQE Molecular Explorer",
) -> object:
    """Init GLFW window with OpenGL 3.3 core profile, depth, MSAA, blending."""
    ...


def create_orbital_vao(
    mesh: IsosurfaceMesh,
) -> tuple[int, int, int, int]:
    """Upload IsosurfaceMesh to GPU.

    Vertex layout: [x, y, z, nx, ny, nz, sign] — 7 floats, 28 bytes.
    Returns (VAO, VBO, EBO, n_indices).
    """
    ...


def update_orbital_vao(
    vbo: int,
    ebo: int,
    mesh: IsosurfaceMesh,
) -> int:
    """Orphan-and-reupload new mesh data to existing VBO/EBO.

    Returns new index count.
    """
    ...


def create_pes_curve_vao(
    bond_lengths: np.ndarray,
    energies_dict: dict[str, np.ndarray],
) -> tuple[int, int]:
    """Create VAO for PES curves as coloured line strips.

    Vertex layout: [R, E, 0, r, g, b, a] — 7 floats per vertex.
    Returns (VAO, VBO).
    """
    ...


def draw_molecule(
    atom_positions: np.ndarray,
    atom_symbols: list[str],
    bonds: list[tuple[int, int]],
    shader_program: int,
) -> None:
    """Render atoms as spheres (CPK colours, 0.3× vdW radii) and bonds as cylinders."""
    ...


class ArcballCamera:
    """Arcball camera: left-drag rotates, scroll zooms.

    Provides view (lookAt) and perspective (45°, near=0.1, far=200)
    matrices for the shader pipeline.
    """

    def __init__(
        self,
        position: np.ndarray | None = None,
        target: np.ndarray | None = None,
    ) -> None: ...

    def on_mouse_drag(self, dx: float, dy: float) -> None: ...

    def on_scroll(self, delta: float) -> None: ...

    def get_view_matrix(self) -> np.ndarray: ...

    def get_projection_matrix(self, aspect: float, fov: float = 45.0) -> np.ndarray: ...
```

*(Full implementation: `report2_scripts/renderer.py`)*

#### 10. `hardware_runner.py`

**Purpose:** Submit VQE circuits to IBM Quantum or IonQ via PennyLane device plugins.

- **`HardwareVQEResult`** — Data class storing the measured energy, raw bitstring counts, shot count, hardware backend name, and the number of commuting Pauli groups measured.

| Function | Description |
|:---------|:------------|
| `create_ibm_device(n_qubits, backend_name, shots, ibm_token)` | Creates a PennyLane device backed by the PennyLane-Qiskit plugin. Authenticates with an IBM Quantum API token (from argument or `IBMQ_TOKEN` environment variable). Default backend: `ibm_brisbane`, default shots: 8192. |
| `create_ionq_device(n_qubits, backend, shots, api_key)` | Creates a PennyLane device backed by the PennyLane-IonQ plugin. Default backend: `ionq_harmony` (simulator); switch to `ionq.qpu` for real hardware. Authenticates via argument or `IONQ_API_KEY` environment variable. |
| `estimate_shot_budget(hamiltonian, target_precision)` | Estimates the total number of shots $S$ required to achieve standard deviation $\leq \epsilon$ on the energy estimate, using the variance bound $\text{Var}[\hat{E}] \leq \sum_k |c_k|^2 / S$ (assuming independent Pauli measurements). This gives $S \geq \sum_k |c_k|^2 / \epsilon^2$. For H₂ with $\sum|c_k|^2 \approx 0.3$ and $\epsilon = 1.6 \times 10^{-3}$ Ha, the bound yields $S \approx 117{,}000$. The result is clamped to $[1000, 10^7]$. |

```python
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pennylane as qml


@dataclass
class HardwareVQEResult:
    """Results from a hardware VQE execution."""
    energy: float
    raw_counts: dict[str, int]
    shots: int
    backend_name: str
    n_pauli_groups: int


def create_ibm_device(
    n_qubits: int,
    backend_name: str = "ibm_brisbane",
    shots: int = 8192,
    ibm_token: str | None = None,
) -> qml.Device:
    """Create PennyLane device via qiskit.ibmq plugin.

    Authenticates with ibm_token or IBMQ_TOKEN env var.
    """
    ...


def create_ionq_device(
    n_qubits: int,
    backend: str = "ionq_harmony",
    shots: int = 1024,
    api_key: str | None = None,
) -> qml.Device:
    """Create PennyLane device via ionq.simulator plugin.

    Authenticates with api_key or IONQ_API_KEY env var.
    """
    ...


def estimate_shot_budget(
    hamiltonian: qml.Hamiltonian,
    target_precision: float = 1.6e-3,
) -> int:
    """Estimate shots S ≥ Σ|c_k|² / ε² for energy precision ε.

    Result clamped to [1000, 10⁷].
    """
    ...
```

*(Full implementation: `report2_scripts/hardware_runner.py`)*

#### 11. `benchmarker.py`

**Purpose:** Compare VQE results against classical methods across the PES.

- **`BenchmarkResult`** — Data class holding scanned bond lengths, list of method names, a dictionary of energy arrays keyed by method, correlation energy recovery ratios $r_{\text{corr}} = (E_{\text{method}} - E_{\text{HF}}) / (E_{\text{FCI}} - E_{\text{HF}})$, and absolute errors $|E_{\text{method}} - E_{\text{FCI}}|$ for each method at every bond length.

> **`run_benchmark(pes_data)`** takes a `PESData` object (from `scan_pes`) and computes, for every available method, the correlation energy recovery ratio and the absolute error relative to FCI. Division-by-zero is guarded at geometries where $E_{\text{FCI}} \approx E_{\text{HF}}$.

> **`print_benchmark_table(result, bond_lengths_subset=None)`** formats the benchmark into a human-readable table showing $R$, each method's energy, and the VQE-vs-FCI error. If no subset is specified, five evenly spaced bond lengths are selected automatically.

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BenchmarkResult:
    """Results from a multi-method PES benchmark."""
    bond_lengths: np.ndarray
    methods: list[str]
    energies: dict[str, np.ndarray]
    correlation_recovery: dict[str, np.ndarray]
    errors_vs_fci: dict[str, np.ndarray]


def run_benchmark(pes_data: PESData) -> BenchmarkResult:
    """Compute correlation recovery ratio and |E - E_FCI| for each method.

    Division-by-zero guarded where |E_FCI - E_HF| < 1e-10.
    """
    ...


def print_benchmark_table(
    result: BenchmarkResult,
    bond_lengths_subset: np.ndarray | None = None,
) -> str:
    """Format benchmark as a human-readable table.

    Selects 5 evenly spaced bond lengths if no subset given.
    """
    ...
```

*(Full implementation: `report2_scripts/benchmarker.py`)*

#### 12. `main.py`

**Purpose:** CLI entry point and workflow orchestration.

> **`parse_args()`** defines the CLI interface via `argparse`. Supported flags include `--molecule` (`h2`, `lih`, `h2o`), `--basis` (e.g., `sto-3g`, `6-31g*`), `--bond-range` (`start,stop,n_points`), `--ansatz` (`uccsd` or `hardware_efficient`), `--optimizer` (`Adam`, `GradientDescent`, `COBYLA`, `SPSA`), `--maxiter`, `--hardware` / `--backend` for real-device execution, `--render` to enable the OpenGL viewer, `--benchmark` to compute all classical reference methods, `--orbital` for the MO index to visualise, and `--config` for a YAML configuration file.

> **`run_pipeline(args)`** orchestrates the full workflow:
>
> 1. Parses the bond-range string into a `linspace` array.
> 2. Selects the appropriate molecule factory function.
> 3. Calls `scan_pes` to sweep bond lengths with VQE.
> 4. If `--benchmark` is set, calls `run_benchmark` and `print_benchmark_table`.
> 5. If `--render` is set, runs HF at equilibrium, evaluates the selected MO on a 3D grid, extracts the dual isosurface, and launches the OpenGL render loop with arcball camera interaction.
> 6. Prints a summary: equilibrium VQE energy, bond length, FCI energy, and absolute error.

```python
from __future__ import annotations

import argparse

import numpy as np


def parse_args() -> argparse.Namespace:
    """Define CLI: --molecule, --basis, --bond-range, --ansatz, --optimizer,
    --maxiter, --hardware, --backend, --render, --benchmark, --orbital, --config.
    """
    ...


def run_pipeline(args: argparse.Namespace) -> None:
    """Orchestrate: parse bond range → select factory → scan PES →
    benchmark (optional) → orbital render (optional) → print summary.
    """
    ...
```

*(Full implementation: `report2_scripts/main.py`)*

---

## Part 3 — Implementation Guidelines (Instructional Pseudocode)

### Module 1: `molecule.py`

```python
# --- molecule.py pseudocode ---

def create_h2(bond_length=0.735):
    # Place two H on z-axis at ±R/2
    atoms = [
        Atom("H", np.array([0.0, 0.0, -bond_length / 2])),
        Atom("H", np.array([0.0, 0.0, +bond_length / 2])),
    ]
    return MoleculeSpec(atoms, basis="sto-3g", charge=0,
                        multiplicity=1, name=f"H2_R{bond_length:.3f}")


def create_lih(bond_length=1.595):
    # Li at origin, H at (0,0,R)
    atoms = [
        Atom("Li", np.array([0.0, 0.0, 0.0])),
        Atom("H",  np.array([0.0, 0.0, bond_length])),
    ]
    return MoleculeSpec(atoms, basis="sto-3g", name=f"LiH_R{bond_length:.3f}")


def create_h2o(oh_length=0.957, angle_deg=104.5):
    alpha = np.radians(angle_deg)
    # H positions: (±R*sin(α/2), 0, -R*cos(α/2))
    h1 = np.array([+oh_length * np.sin(alpha / 2), 0.0,
                    -oh_length * np.cos(alpha / 2)])
    h2 = np.array([-oh_length * np.sin(alpha / 2), 0.0,
                    -oh_length * np.cos(alpha / 2)])
    atoms = [Atom("O", np.zeros(3)), Atom("H", h1), Atom("H", h2)]
    return MoleculeSpec(atoms, basis="sto-3g", name=f"H2O_R{oh_length:.3f}")


def set_bond_length(mol, atom_idx_a, atom_idx_b, new_length):
    new_mol = copy.deepcopy(mol)
    pa = new_mol.atoms[atom_idx_a].position  # shape: (3,)
    pb = new_mol.atoms[atom_idx_b].position  # shape: (3,)
    midpoint = (pa + pb) / 2.0               # shape: (3,)
    d_hat = (pb - pa)                        # direction vector
    d_hat = d_hat / np.linalg.norm(d_hat)    # unit vector
    # Reposition symmetrically about midpoint
    new_mol.atoms[atom_idx_a].position = midpoint - d_hat * new_length / 2
    new_mol.atoms[atom_idx_b].position = midpoint + d_hat * new_length / 2
    return new_mol
```

*(Full implementation: `report2_scripts/molecule.py`)*

### Module 2: `classical_hf.py`

```python
# --- classical_hf.py pseudocode ---

def _build_pyscf_mol(mol_spec):
    """Shared helper: MoleculeSpec → PySCF Mole."""
    # Format atom string: "H 0.0 0.0 0.0; H 0.0 0.0 0.735"
    atom_str = "; ".join(
        f"{a.symbol} {a.position[0]:.6f} {a.position[1]:.6f} {a.position[2]:.6f}"
        for a in mol_spec.atoms
    )
    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = mol_spec.basis
    mol.charge = mol_spec.charge
    mol.spin = mol_spec.multiplicity - 1  # PySCF uses 2S, not 2S+1
    mol.unit = "Angstrom"
    mol.build()
    return mol


def run_hartree_fock(mol_spec):
    mol = _build_pyscf_mol(mol_spec)

    # Run restricted Hartree-Fock
    mf = scf.RHF(mol)
    mf.kernel()  # converge SCF

    mo_coeff = mf.mo_coeff      # shape: (n_ao, n_mo) — AO→MO matrix
    mo_energy = mf.mo_energy    # shape: (n_mo,)
    n_mo = mo_coeff.shape[1]

    # One-electron integrals: AO → MO via C^T h^core C
    h1_ao = mf.get_hcore()                  # shape: (n_ao, n_ao)
    h1_mo = mo_coeff.T @ h1_ao @ mo_coeff   # shape: (n_mo, n_mo)

    # Two-electron integrals: efficient AO→MO with ao2mo.full
    h2_mo = ao2mo.full(mol, mo_coeff)       # packed triangular form
    h2_mo = ao2mo.restore(1, h2_mo, n_mo)   # shape: (n_mo, n_mo, n_mo, n_mo)

    overlap = mol.intor("int1e_ovlp")       # shape: (n_ao, n_ao)
    e_nuc = mol.energy_nuc()                # scalar (Hartree)

    return HFResult(
        energy_hf=mf.e_tot, energy_nuc=e_nuc,
        mo_coefficients=mo_coeff, mo_energies=mo_energy,
        one_electron_integrals=h1_mo, two_electron_integrals=h2_mo,
        overlap_matrix=overlap,
        n_electrons=mol.nelectron, n_orbitals=n_mo,
        ao_labels=[str(l) for l in mol.ao_labels()],
    )


def get_fci_energy(mol_spec):
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol).run()       # RHF as starting point
    cisolver = fci.FCI(mf)
    e_fci, _ = cisolver.kernel()   # exact ground-state energy
    return float(e_fci)


def get_cisd_energy(mol_spec):
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol).run()
    myci = ci.CISD(mf).run()
    return float(myci.e_tot)       # E_CISD (Hartree)


def get_ccsd_t_energy(mol_spec):
    mol = _build_pyscf_mol(mol_spec)
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf).run()
    # Perturbative triples correction
    e_t = mycc.ccsd_t()            # ΔE_(T)
    return float(mycc.e_tot + e_t)  # E_CCSD(T) = E_CCSD + ΔE_(T)
```

*(Full implementation: `report2_scripts/classical_hf.py`)*

### Module 3: `hamiltonian_builder.py`

```python
# --- hamiltonian_builder.py pseudocode ---

def build_hamiltonian(h1, h2, n_electrons, e_nuc,
                      mapping="jordan_wigner",
                      active_electrons=None, active_orbitals=None):
    n_mo = h1.shape[0]
    # n_spin_orbs = 2 * n_mo  (each spatial orbital → 2 spin-orbitals)

    # --- Stage 1: Build fermionic Hamiltonian ---
    # One-electron terms: h_ij a+_{2i+σ} a_{2j+σ}  (spin-diagonal)
    fermi_h = 0
    for i in range(n_mo):
        for j in range(n_mo):
            if abs(h1[i, j]) < 1e-12:
                continue
            for sigma in (0, 1):  # alpha, beta
                p, q = 2 * i + sigma, 2 * j + sigma
                fermi_h += h1[i, j] * FermiWord({(0, p): "+", (1, q): "-"})

    # Two-electron terms: ½ h_ijkl a+_p a+_q a_s a_r (all spin combos, skip p==q)
    for i in range(n_mo):
        for j in range(n_mo):
            for k in range(n_mo):
                for l in range(n_mo):
                    if abs(h2[i, j, k, l]) < 1e-12:
                        continue
                    for s1 in (0, 1):
                        for s2 in (0, 1):
                            p = 2 * i + s1
                            q = 2 * k + s2
                            r = 2 * j + s1
                            s = 2 * l + s2
                            if p == q:
                                continue
                            fermi_h += 0.5 * h2[i, j, k, l] * FermiWord(
                                {(0, p): "+", (1, q): "+", (2, s): "-", (3, r): "-"}
                            )

    # Nuclear repulsion: constant * Identity
    # fermi_h += e_nuc * Identity

    # --- Stage 2: Fermion-to-qubit mapping ---
    if mapping == "jordan_wigner":
        qubit_h = qml.jordan_wigner(fermi_h)   # JW: a+_p → ½(X-iY)∏Z
    else:
        qubit_h = qml.parity_transform(fermi_h)  # parity encoding
    qubit_h = qml.simplify(qubit_h)  # merge & cancel redundant Pauli terms
    # ... wrap into QubitHamiltonian dataclass ...

    # Alternative (simpler): one-call PennyLane interface
    # H, n_qubits = qml.qchem.molecular_hamiltonian(
    #     symbols, coords, basis=basis, mapping=mapping)


def pennylane_to_qutip(hamiltonian, n_qubits):
    """Cross-validate by converting to QuTiP and diagonalising."""
    pauli_map = {"I": qutip.qeye(2), "X": qutip.sigmax(),
                 "Y": qutip.sigmay(), "Z": qutip.sigmaz()}
    H_qt = 0
    for coeff, obs in zip(hamiltonian.coeffs, hamiltonian.ops):
        # Decompose obs into per-qubit Pauli labels
        factors = [pauli_map["I"]] * n_qubits
        for wire, pauli_label in _decompose_pauli(obs):
            factors[wire] = pauli_map[pauli_label]
        # Build tensor product: ⊗_q σ_q
        term = qutip.tensor(factors)        # dim: (2^N, 2^N)
        H_qt = H_qt + coeff * term
    # Exact diagonalisation
    eigenvalues = H_qt.eigenenergies()      # sorted ascending
    e0 = float(eigenvalues[0])              # ground-state energy
    return e0
```

*(Full implementation: `report2_scripts/hamiltonian_builder.py`)*

### Module 4: `ansatz.py`

```python
# --- ansatz.py pseudocode ---

def build_uccsd_ansatz(n_qubits, n_electrons, hf_state=None):
    # Prepare HF reference state: first n_electrons orbitals occupied
    if hf_state is None:
        hf_state = np.zeros(n_qubits, dtype=int)
        hf_state[:n_electrons] = 1   # e.g. |1100⟩ for H₂

    # Enumerate spin-allowed excitations via PennyLane qchem
    singles, doubles = qml.qchem.excitations(n_electrons, n_qubits)
    # H₂ example: singles = [[0,2], [1,3]], doubles = [[0,1,2,3]]
    excitations = singles + doubles
    n_params = len(excitations)             # = |singles| + |doubles|

    def ansatz_fn(params):                  # params shape: (n_params,)
        qml.BasisState(hf_state, wires=range(n_qubits))

        idx = 0
        for exc in singles:
            # exp(iθ (a+_a a_i − h.c.)) via first-order Trotter
            qml.SingleExcitation(params[idx], wires=exc)
            idx += 1
        for exc in doubles:
            # exp(iθ (a+_a a+_b a_j a_i − h.c.))
            qml.DoubleExcitation(params[idx], wires=exc)
            idx += 1

    # NOTE: Excitation ordering determines Trotter decomposition.
    # Optimizer compensates for Trotter error.
    # For H₂: by symmetry, singles are zero at the optimum;
    # only the double excitation matters (1 effective parameter).
    return ansatz_fn, n_params, excitations


def build_hardware_efficient_ansatz(n_qubits, n_layers=2, entangler="CNOT"):
    n_params = 2 * n_qubits * n_layers      # RY + RZ per qubit per layer

    def ansatz_fn(params):                   # params shape: (n_params,)
        p = params.reshape(n_layers, n_qubits, 2)  # shape: (L, N_q, 2)

        for l in range(n_layers):
            # Single-qubit rotation layer
            for q in range(n_qubits):
                qml.RY(p[l, q, 0], wires=q)
                qml.RZ(p[l, q, 1], wires=q)
            # Entangling ladder (nearest-neighbour)
            for q in range(n_qubits - 1):
                if entangler == "CNOT":
                    qml.CNOT(wires=[q, q + 1])
                elif entangler == "CZ":
                    qml.CZ(wires=[q, q + 1])

    return ansatz_fn, n_params
```

*(Full implementation: `report2_scripts/ansatz.py`)*

### Module 5: `vqe_solver.py`

```python
# --- vqe_solver.py pseudocode ---

def run_vqe(hamiltonian, ansatz_fn, n_params, n_qubits,
            init_params=None, optimizer="Adam", learning_rate=0.05,
            maxiter=200, tol=1e-6, device_name="default.qubit",
            diff_method="backprop", callback=None):

    dev = qml.device(device_name, wires=n_qubits)

    @qml.qnode(dev, diff_method=diff_method)
    def cost_fn(params):                    # params shape: (n_params,)
        ansatz_fn(params)
        return qml.expval(hamiltonian)      # <ψ(θ)|H|ψ(θ)>

    # Initialise with small random perturbation to break symmetry
    if init_params is None:
        init_params = np.random.normal(0, 0.01, size=n_params)  # shape: (n_params,)
    params = init_params.copy()
    trajectory = []                         # list of (params_copy, energy)
    gradient_norms = []

    if optimizer in ("Adam", "GradientDescent"):
        # PennyLane built-in optimiser (autograd / backprop)
        opt = qml.AdamOptimizer(stepsize=learning_rate)  # or GradientDescentOptimizer
        for k in range(maxiter):
            params, prev_e = opt.step_and_cost(cost_fn, params)
            energy = float(cost_fn(params))
            grad = qml.grad(cost_fn)(params)          # shape: (n_params,)
            g_norm = float(np.linalg.norm(grad))
            trajectory.append((params.copy(), energy))
            gradient_norms.append(g_norm)
            if callback:
                callback(k, params, energy)
            # Convergence: |ΔE| < tol
            if len(trajectory) > 1 and abs(trajectory[-1][1] - trajectory[-2][1]) < tol:
                break

    elif optimizer == "COBYLA":
        # SciPy derivative-free optimiser
        def objective(p):
            e = float(cost_fn(p))
            trajectory.append((p.copy(), e))
            return e
        scipy.optimize.minimize(objective, params, method="COBYLA",
                                options={"maxiter": maxiter, "rhobeg": 0.5})

    elif optimizer == "SPSA":
        # Simultaneous Perturbation Stochastic Approximation
        a0, c0, A, alpha, gamma = 0.1, 0.1, 10, 0.602, 0.101
        for k in range(1, maxiter + 1):
            a_k = a0 / (k + A) ** alpha       # step size (decaying)
            c_k = c0 / k ** gamma             # perturbation size (decaying)
            delta = np.random.choice([-1, 1], size=n_params)  # ±1 vector
            f_plus  = float(cost_fn(params + c_k * delta))
            f_minus = float(cost_fn(params - c_k * delta))
            # Gradient estimate: 2 circuit evals regardless of n_params
            g_hat = (f_plus - f_minus) / (2 * c_k * delta)  # shape: (n_params,)
            params = params - a_k * g_hat
            trajectory.append((params.copy(), f_plus))
            gradient_norms.append(float(np.linalg.norm(g_hat)))

    # Select result with lowest energy from full trajectory
    best_idx = int(np.argmin([e for _, e in trajectory]))
    best_params, best_energy = trajectory[best_idx]
    # ... return VQEResult ...
```

*(Full implementation: `report2_scripts/vqe_solver.py`)*

### Module 6: `pes_scanner.py`

```python
# --- pes_scanner.py pseudocode ---

def scan_pes(molecule_factory, bond_lengths, basis="sto-3g",
             mapping="jordan_wigner", ansatz_type="uccsd",
             optimizer="Adam", maxiter=200,
             use_prev_params=True, compute_classical=True, callback=None):

    n_points = len(bond_lengths)
    energies_vqe = np.zeros(n_points)       # shape: (N,)
    energies_hf  = np.zeros(n_points)       # shape: (N,)
    energies_fci = np.zeros(n_points)       # shape: (N,)
    prev_params = None

    for idx, R in enumerate(bond_lengths):
        # 1. Build geometry at bond length R_k
        mol_spec = molecule_factory(R)
        mol_spec.basis = basis

        # 2. Run Hartree-Fock → MO integrals
        hf_result = run_hartree_fock(mol_spec)
        energies_hf[idx] = hf_result.energy_hf

        # 3. Build qubit Hamiltonian (JW or parity)
        qham = build_hamiltonian(
            hf_result.one_electron_integrals,    # shape: (M, M)
            hf_result.two_electron_integrals,    # shape: (M, M, M, M)
            hf_result.n_electrons,
            hf_result.energy_nuc,
            mapping=mapping,
        )

        # 4. Build ansatz (UCCSD or hardware-efficient)
        if ansatz_type == "uccsd":
            ansatz_fn, n_params, _ = build_uccsd_ansatz(
                qham.n_qubits, hf_result.n_electrons)
        else:
            ansatz_fn, n_params = build_hardware_efficient_ansatz(
                qham.n_qubits, n_layers=2)

        # 5. Warm-start: reuse previous geometry's optimal params
        #    (only if param count unchanged; smooth PES → smooth θ*)
        init_p = None
        if use_prev_params and prev_params is not None:
            if len(prev_params) == n_params:
                init_p = prev_params  # warm-start halves iterations

        # 6. Run VQE
        vqe_result = run_vqe(
            qham.hamiltonian, ansatz_fn, n_params, qham.n_qubits,
            init_params=init_p, optimizer=optimizer, maxiter=maxiter)
        energies_vqe[idx] = vqe_result.optimal_energy
        prev_params = vqe_result.optimal_params

        # 7. Classical benchmarks (optional)
        if compute_classical:
            energies_fci[idx] = get_fci_energy(mol_spec)
            # ... similarly for CISD, CCSD(T) ...

    # Caveat: at large R, RHF may converge to a broken-symmetry solution.
    # Active space must remain fixed across the scan for warm-start.
    return PESData(bond_lengths=bond_lengths, energies_vqe=energies_vqe,
                   energies_hf=energies_hf, energies_fci=energies_fci, ...)
```

*(Full implementation: `report2_scripts/pes_scanner.py`)*

### Module 7: `orbital_grid.py`

```python
# --- orbital_grid.py pseudocode ---

def evaluate_orbital_on_grid(mo_coefficients, orbital_index,
                              atom_coords, atom_symbols, basis_name,
                              grid_extent=5.0, grid_spacing=0.1):
    # 1. Build PySCF Mole for AO evaluation
    mol = gto.Mole()
    atom_str = "; ".join(
        f"{sym} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}"
        for sym, c in zip(atom_symbols, atom_coords)
    )
    mol.atom = atom_str
    mol.basis = basis_name
    mol.unit = "Bohr"           # PySCF uses Bohr internally
    mol.build()

    # 2. Construct uniform 3D grid centred on molecular centroid
    center = np.mean(atom_coords, axis=0)    # shape: (3,)
    x = np.arange(center[0] - grid_extent,
                   center[0] + grid_extent, grid_spacing)
    y = np.arange(...)                        # analogous for y, z
    z = np.arange(...)
    nx, ny, nz = len(x), len(y), len(z)      # typically ~100 each

    # Flatten to (N_grid, 3) for batch evaluation
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    coords = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)
    # coords shape: (nx*ny*nz, 3)

    # 3. Evaluate all K AO basis functions at all grid points
    ao_values = mol.eval_gto("GTOval_cart", coords)
    # ao_values shape: (nx*ny*nz, K) — exploits Gaussian sparsity internally

    # 4. Contract with MO coefficient column p: φ_p(r) = Σ_μ C_{μp} χ_μ(r)
    mo_col = mo_coefficients[:, orbital_index]   # shape: (K,)
    mo_values = ao_values @ mo_col               # shape: (nx*ny*nz,)

    # 5. Reshape to 3D volume
    volume = mo_values.reshape(nx, ny, nz)       # shape: (N_x, N_y, N_z)

    return OrbitalVolumeData(
        values=volume,
        grid_origin=(float(x[0]), float(y[0]), float(z[0])),
        grid_spacing=(grid_spacing, grid_spacing, grid_spacing),
        grid_shape=(nx, ny, nz),
        orbital_index=orbital_index,
    )


def evaluate_electron_density(mo_coefficients, n_occupied, ...):
    # Electron density: ρ(r) = Σ_{i=0}^{n_occ-1} 2 * |φ_i(r)|²
    # Factor of 2: RHF double occupation
    density = np.zeros((nx, ny, nz))             # shape: (N_x, N_y, N_z)
    for i in range(n_occupied):
        orb_data = evaluate_orbital_on_grid(mo_coefficients, i, ...)
        density += 2.0 * orb_data.values ** 2
    return OrbitalVolumeData(values=density, ..., orbital_index=-1)

    # UNIT CAVEAT: If atom_coords are in Ångström, convert first:
    # atom_coords_bohr = atom_coords * 1.8897259886
```

*(Full implementation: `report2_scripts/orbital_grid.py`)*

### Module 8: `isosurface.py`

```python
# --- isosurface.py pseudocode ---

def marching_cubes(volume, isovalue, grid_origin, grid_spacing):
    """
    Lorensen & Cline (1987) via scikit-image.
    1. Each voxel (2×2×2 corners) classified inside/outside → 8-bit index
    2. Index selects triangulation from 256-entry lookup table (15 unique)
    3. Edge crossings linearly interpolated between corner values
    4. Vertex normals from scalar field gradient
    """
    # scikit-image returns (verts, faces, normals, values)
    verts, faces, normals, _ = sk_marching_cubes(
        volume,                          # shape: (N_x, N_y, N_z)
        level=isovalue,
        spacing=grid_spacing,            # (dx, dy, dz) in Bohr
    )
    # verts shape: (N_v, 3) in grid-spacing units
    # faces shape: (N_f, 3) triangle vertex indices
    # normals shape: (N_v, 3)

    # Translate to world coordinates
    origin = np.array(grid_origin)           # shape: (3,)
    verts = verts + origin[None, :]          # broadcast add

    signs = np.ones(len(verts), dtype=np.float32)  # +1 for positive lobe
    return IsosurfaceMesh(
        vertices=verts.astype(np.float32),
        normals=normals.astype(np.float32),
        faces=faces.astype(np.uint32),
        signs=signs,
    )


def extract_dual_isosurface(volume, isovalue, grid_origin, grid_spacing):
    """Extract both +c and -c surfaces for positive/negative lobes."""
    # Positive lobe (sign = +1)
    mesh_pos = marching_cubes(volume, +isovalue, grid_origin, grid_spacing)
    # Negative lobe (sign = -1)
    mesh_neg = marching_cubes(volume, -isovalue, grid_origin, grid_spacing)
    mesh_neg.signs[:] = -1.0                 # tag as negative lobe

    n_verts_pos = len(mesh_pos.vertices)

    # Concatenate meshes; offset negative face indices by |V_pos|
    combined_verts   = np.concatenate([mesh_pos.vertices, mesh_neg.vertices])
    combined_normals = np.concatenate([mesh_pos.normals, mesh_neg.normals])
    combined_signs   = np.concatenate([mesh_pos.signs, mesh_neg.signs])
    neg_faces_offset = mesh_neg.faces + n_verts_pos   # index offset
    combined_faces   = np.concatenate([mesh_pos.faces, neg_faces_offset])

    return IsosurfaceMesh(vertices=combined_verts, normals=combined_normals,
                          faces=combined_faces, signs=combined_signs)

    # EDGE CASE: marching_cubes fails if isovalue outside data range.
    # For negative surface, need volume.min() < -c AND volume.max() > -c.
    # Very small isovalues near 0 can produce degenerate triangles.
```

*(Full implementation: `report2_scripts/isosurface.py`)*

### Module 9: `renderer.py`

```python
# --- renderer.py pseudocode ---

def init_opengl_context(width=1400, height=900, title="VQE Molecular Explorer"):
    glfw.init()
    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
    glfw.window_hint(glfw.SAMPLES, 4)        # 4× MSAA
    window = glfw.create_window(width, height, title, None, None)
    glfw.make_context_current(window)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_MULTISAMPLE)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    return window


def create_orbital_vao(mesh):
    """Upload IsosurfaceMesh to GPU.
    Vertex layout: [x, y, z, nx, ny, nz, sign] — 7 float32, 28 bytes/vertex
    Position: location 0, Normal: location 1, Sign: location 2.
    Uses GL_DYNAMIC_DRAW (mesh changes with bond length).
    """
    n_verts = len(mesh.vertices)
    vertices = np.zeros((n_verts, 7), dtype=np.float32)
    vertices[:, 0:3] = mesh.vertices         # position (x, y, z)
    vertices[:, 3:6] = mesh.normals          # normal (nx, ny, nz)
    vertices[:, 6]   = mesh.signs            # lobe sign (±1)

    vao = glGenVertexArrays(1)
    vbo = glGenBuffers(1)
    ebo = glGenBuffers(1)
    glBindVertexArray(vao)

    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_DYNAMIC_DRAW)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.faces.nbytes, mesh.faces, GL_DYNAMIC_DRAW)

    stride = 7 * 4                           # 28 bytes
    # Attribute 0: position (3 floats, offset 0)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
    glEnableVertexAttribArray(0)
    # Attribute 1: normal (3 floats, offset 12)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(12))
    glEnableVertexAttribArray(1)
    # Attribute 2: sign (1 float, offset 24)
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(24))
    glEnableVertexAttribArray(2)

    glBindVertexArray(0)
    return vao, vbo, ebo, len(mesh.faces) * 3


def update_orbital_vao(vbo, ebo, mesh):
    """Orphan old buffer → upload new mesh data (avoids GPU sync stalls)."""
    # ... rebuild vertex array same layout as create_orbital_vao ...
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_DYNAMIC_DRAW)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.faces.nbytes, mesh.faces, GL_DYNAMIC_DRAW)
    return len(mesh.faces) * 3


def create_pes_curve_vao(bond_lengths, energies_dict):
    """Each method (HF, CISD, VQE, FCI) → coloured GL_LINE_STRIP.
    Vertex layout: [R, E, 0, r, g, b, a] — 7 floats per vertex.
    """
    METHOD_COLORS = {"HF": (0.5, 0.5, 0.5, 1), "CISD": (0, 0.8, 0, 1),
                     "VQE": (0.2, 0.4, 0.9, 1), "FCI": (1, 0.84, 0, 1),
                     "CCSD(T)": (0.6, 0.2, 0.8, 1)}
    all_verts = []
    for method, energies in energies_dict.items():
        rgba = METHOD_COLORS.get(method, (0.7, 0.7, 0.7, 1))
        for R, E in zip(bond_lengths, energies):
            all_verts.append([R, E, 0.0, *rgba])  # shape per vert: (7,)
    # ... upload to VBO ...


def draw_molecule(atom_positions, atom_symbols, bonds, shader_program):
    # CPK colours and scaled vdW radii
    ATOM_COLORS = {"H": (1,1,1), "C": (0.5,0.5,0.5), "O": (1,0,0), "Li": (0.6,0.2,0.8)}
    VDW_RADII   = {"H": 1.20, "C": 1.70, "O": 1.52, "Li": 1.82}
    for pos, sym in zip(atom_positions, atom_symbols):
        radius = VDW_RADII.get(sym, 1.5) * 0.3   # scaled ×0.3
        model = _make_sphere_model_matrix(pos, radius)
        # ... set uniforms, draw instanced icosphere ...
    for i, j in bonds:
        model = _make_cylinder_model_matrix(atom_positions[i], atom_positions[j], 0.1)
        # ... draw unit cylinder ...


class ArcballCamera:
    """Arcball: yaw/pitch from mouse drag, zoom from scroll."""
    def __init__(self, position=None, target=None):
        self.yaw, self.pitch = -90.0, 20.0
        self.zoom = 10.0
        # ...

    def on_mouse_drag(self, dx, dy):
        self.yaw   += dx * 0.3               # Δyaw = 0.3 * dx
        self.pitch -= dy * 0.3               # Δpitch = -0.3 * dy
        self.pitch  = np.clip(self.pitch, -89, 89)
        self._update()

    def on_scroll(self, delta):
        self.zoom -= delta * 0.5
        self.zoom = max(1.0, self.zoom)
        self._update()

    def _update(self):
        # Recompute position in spherical coords relative to target
        yr, pr = np.radians(self.yaw), np.radians(self.pitch)
        self.position = self.target + self.zoom * np.array([
            np.cos(pr) * np.cos(yr), np.sin(pr), np.cos(pr) * np.sin(yr)])

    def get_view_matrix(self):
        return _look_at(self.position, self.target, up=[0, 1, 0])

    def get_projection_matrix(self, aspect, fov=45.0):
        return _perspective(fov, aspect, near=0.1, far=200.0)
```

*(Full implementation: `report2_scripts/renderer.py`)*

### Module 10: `hardware_runner.py`

```python
# --- hardware_runner.py pseudocode ---

def create_ibm_device(n_qubits, backend_name="ibm_brisbane",
                      shots=8192, ibm_token=None):
    token = ibm_token or os.environ.get("IBMQ_TOKEN")
    dev = qml.device(
        "qiskit.ibmq",
        wires=n_qubits,
        backend=backend_name,
        shots=shots,
        ibmqx_token=token,
    )
    return dev


def create_ionq_device(n_qubits, backend="ionq_harmony",
                       shots=1024, api_key=None):
    key = api_key or os.environ.get("IONQ_API_KEY")
    dev = qml.device(
        "ionq.simulator",           # or "ionq.qpu" for real hardware
        wires=n_qubits,
        shots=shots,
        api_key=key,
    )
    return dev


def estimate_shot_budget(hamiltonian, target_precision=1.6e-3):
    # Variance bound: Var[E] ≤ Σ|c_k|² / S  (independent measurements)
    coeffs = np.array(hamiltonian.coeffs)    # shape: (n_paulis,)
    variance_bound = np.sum(coeffs ** 2)
    # Need: √Var ≤ ε ⇒ S ≥ Σ|c_k|² / ε²
    shots_needed = int(np.ceil(variance_bound / target_precision ** 2))
    shots_needed = max(1000, min(shots_needed, 10_000_000))  # clamp
    return shots_needed
    # H₂ (STO-3G): Σ|c_k|² ≈ 0.3, ε = 1.6e-3 ⇒ S ≈ 117,000
    # Pauli grouping into commuting cliques reduces actual requirement.
```

*(Full implementation: `report2_scripts/hardware_runner.py`)*

### Module 11: `benchmarker.py`

```python
# --- benchmarker.py pseudocode ---

def run_benchmark(pes_data):
    e_hf  = pes_data.energies_hf       # shape: (N,)
    e_fci = pes_data.energies_fci      # shape: (N,)
    denom = e_fci - e_hf               # shape: (N,) — negative values

    methods = ["HF", "VQE"]
    energies = {"HF": e_hf, "VQE": pes_data.energies_vqe, "FCI": e_fci}
    if pes_data.energies_cisd is not None:
        methods.append("CISD")
        energies["CISD"] = pes_data.energies_cisd
    if pes_data.energies_ccsd_t is not None:
        methods.append("CCSD(T)")
        energies["CCSD(T)"] = pes_data.energies_ccsd_t

    correlation_recovery = {}
    errors_vs_fci = {}
    for method in methods:
        e_m = energies[method]                     # shape: (N,)
        # Correlation energy recovery: r = (E_method - E_HF) / (E_FCI - E_HF)
        # Guard division-by-zero where |E_FCI - E_HF| < 1e-10
        recovery = np.where(
            np.abs(denom) > 1e-10,
            (e_m - e_hf) / denom,
            0.0,
        )                                          # shape: (N,)
        correlation_recovery[method] = recovery
        # Absolute error vs FCI: |E_method - E_FCI|
        errors_vs_fci[method] = np.abs(e_m - e_fci)  # shape: (N,)

    return BenchmarkResult(
        bond_lengths=pes_data.bond_lengths,
        methods=methods, energies=energies,
        correlation_recovery=correlation_recovery,
        errors_vs_fci=errors_vs_fci,
    )


def print_benchmark_table(result, bond_lengths_subset=None):
    # Select 5 evenly spaced bond lengths if no subset given
    if bond_lengths_subset is None:
        indices = np.linspace(0, len(result.bond_lengths) - 1, 5, dtype=int)
    else:
        indices = [np.argmin(np.abs(result.bond_lengths - bl))
                   for bl in bond_lengths_subset]

    # Print header: R (Å) | HF | CISD | VQE | FCI | |VQE-FCI|
    header = f"{'R (Å)':>8}"
    for m in result.methods + ["FCI"]:
        header += f" | {m:>12}"
    header += f" | {'|VQE-FCI|':>12}"
    print(header)
    # ... print rows for each bond length in indices ...
```

*(Full implementation: `report2_scripts/benchmarker.py`)*

### Module 12: `main.py`

```python
# --- main.py pseudocode ---

def parse_args():
    parser = argparse.ArgumentParser(description="VQE PES Explorer")
    parser.add_argument("--molecule", default="h2", choices=["h2", "lih", "h2o"])
    parser.add_argument("--basis", default="sto-3g")
    parser.add_argument("--bond-range", default="0.3,3.0,20")
    parser.add_argument("--ansatz", default="uccsd", choices=["uccsd", "hardware_efficient"])
    parser.add_argument("--optimizer", default="Adam",
                        choices=["Adam", "GradientDescent", "COBYLA", "SPSA"])
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--hardware", action="store_true")
    parser.add_argument("--backend", default="ibm_brisbane")
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--orbital", type=int, default=0)    # MO index to visualise
    parser.add_argument("--config", default="config.yaml")
    return parser.parse_args()


def run_pipeline(args):
    # 1. Parse bond range → linspace array
    start, stop, n_pts = args.bond_range.split(",")
    bond_lengths = np.linspace(float(start), float(stop), int(n_pts))

    # 2. Select molecule factory
    factories = {"h2": create_h2, "lih": create_lih, "h2o": create_h2o}
    mol_factory = factories[args.molecule]

    # 3. PES scan
    pes_data = scan_pes(
        mol_factory, bond_lengths,
        basis=args.basis, ansatz_type=args.ansatz,
        optimizer=args.optimizer, maxiter=args.maxiter,
        compute_classical=args.benchmark,
    )

    # 4. Benchmark (optional)
    if args.benchmark:
        bench = run_benchmark(pes_data)
        print_benchmark_table(bench)

    # 5. 3D orbital render (optional)
    if args.render:
        mol_eq = mol_factory()               # equilibrium geometry
        hf_eq = run_hartree_fock(mol_eq)
        # Convert atom coords: Å → Bohr
        atom_coords_bohr = np.array([a.position * 1.8897 for a in mol_eq.atoms])
        orb_data = evaluate_orbital_on_grid(
            hf_eq.mo_coefficients, args.orbital,
            atom_coords_bohr, [a.symbol for a in mol_eq.atoms],
            mol_eq.basis,
        )
        mesh = extract_dual_isosurface(
            orb_data.values, isovalue=0.05,
            orb_data.grid_origin, orb_data.grid_spacing,
        )
        window = init_opengl_context()
        vao, vbo, ebo, n_idx = create_orbital_vao(mesh)
        camera = ArcballCamera()
        # ... enter render loop with arcball interaction ...

    # 6. Print summary
    eq_idx = np.argmin(pes_data.energies_vqe)
    print(f"VQE min energy: {pes_data.energies_vqe[eq_idx]:.6f} Ha "
          f"at R = {bond_lengths[eq_idx]:.3f} Å")
```

*(Full implementation: `report2_scripts/main.py`)*

---

## Part 4 — OpenGL Visualization Design

### 1. Molecular Orbital Isosurface

#### VAO/VBO Layout

Each vertex stores 7 float32 values (28 bytes):

| Offset | Attribute | Size | Description |
|:------:|:---------:|:----:|:------------|
| 0 | Position | 3 floats | $(x, y, z)$ in Bohr |
| 12 | Normal | 3 floats | Surface normal for lighting |
| 24 | Sign | 1 float | $+1$ (positive lobe) or $-1$ (negative lobe) |

The EBO stores triangle indices as uint32.

#### GLSL Shaders for the Isosurface

**Vertex Shader (`orbital.vert`):**

The vertex shader transforms each mesh vertex from object space to clip space via the standard $\mathbf{P} \cdot \mathbf{V} \cdot \mathbf{M}$ pipeline. It accepts three vertex attributes — position $(x, y, z)$, normal $(n_x, n_y, n_z)$, and lobe sign ($\pm 1$) — bound at locations 0, 1, and 2 respectively. The world-space position and the transformed normal (via the inverse-transpose of the model matrix, $\mathbf{n}' = (\mathbf{M}^{-T})\,\mathbf{n}$) are passed as varyings to the fragment stage.

*(Full shader: `report2_scripts/shaders/orbital.vert`)*

**Fragment Shader (`orbital.frag`):**

The fragment shader implements Phong shading with two-sided lighting (using $|\hat{\mathbf{n}} \cdot \hat{\mathbf{l}}|$ instead of $\max(0, \hat{\mathbf{n}} \cdot \hat{\mathbf{l}})$, so that back-faces of orbital lobes are visible). The lobe sign interpolated from the vertex stage selects the base colour: positive lobes are rendered in blue $(0.2, 0.4, 0.9)$ and negative lobes in red $(0.9, 0.2, 0.2)$, blended via `mix(negColor, posColor, (sign + 1)/2)`. A specular highlight with shininess exponent 32 and strength 0.4 is added. The final fragment is output with $\alpha = 0.7$ (semi-transparent) to allow visual inspection of internal structure.

#### Mesh Update on Bond Length Change

When the user moves a bond-length slider:

1. Recompute the MO at the new geometry: `evaluate_orbital_on_grid(...)` with updated atom positions.
2. Re-extract the isosurface: `extract_dual_isosurface(...)`.
3. Upload new mesh data via `update_orbital_vao(vbo, ebo, new_mesh)`.

To avoid stalls, use **double-buffered VBOs**: allocate two VBOs and alternate between them. When a bond-length change triggers recomputation (which can run on a background thread), the new mesh is uploaded to the inactive VBO. Once the upload completes, the active and inactive VBOs are swapped, giving a seamless transition with no pipeline stall.

*(See `update_orbital_vao` in `report2_scripts/renderer.py`)*

### 2. Potential Energy Surface Curve

The PES is rendered as a 2D curve in a 3D viewport, positioned to the side of the molecular orbital. For each method (HF, CISD, VQE, FCI, etc.), bond lengths are mapped to the $x$-axis and energies to the $y$-axis, with $z = 0$. Each method is drawn as a coloured `GL_LINE_STRIP` with a distinct RGBA colour (grey for HF, green for CISD, blue for VQE, gold for FCI, purple for CCSD(T)). Vertices store `[x, y, z, r, g, b, a]` (7 floats). The VBO is uploaded once per PES scan.

*(See `create_pes_curve_vao` in `report2_scripts/renderer.py`)*

A **moving marker** (small sphere) at the current bond length indicates which point on the PES corresponds to the displayed orbital.

### 3. Molecular Geometry

Atoms are rendered as spheres (using a tessellated icosphere VAO) and bonds as cylinders (using a tessellated tube VAO):

- Atom colours follow CPK convention (H=white, O=red, Li=purple, C=grey).
- Radii are scaled van der Waals radii ($\times 0.3$).
- Bond cylinders connect bonded atom pairs with radius $\sim 0.1$ Bohr.
- As the PES scan progresses, atom positions update and the molecular geometry animates smoothly.

### 4. Camera, Lighting, and UI Layout

**Layout:** The window is divided conceptually:
- Left 2/3: 3D viewport with orbital isosurface and molecular geometry.
- Right 1/3: PES curve plot and info panel (via Dear ImGui).

**Lighting:** A single directional light from above-right, plus ambient. The light position is fixed in camera space so it rotates with the view.

**Camera:** Arcball camera (see `ArcballCamera` class) with:
- Left-drag: rotate
- Scroll: zoom
- Middle-drag: pan

**ImGui Panel:** The Dear ImGui overlay renders a side panel displaying the molecule name, basis set, current bond length, VQE and FCI energies, VQE–FCI error, and orbital index. A `slider_float` widget allows the user to drag the bond length between 0.3 and 3.0 Å, triggering `update_visualization(current_R)` whenever the value changes.

*(See `report2_scripts/renderer.py` for the full render loop and ImGui integration.)*

---

## Part 5 — Hardware Execution Plan

### 1. IBM Quantum Submission via PennyLane

#### 1.1 Circuit Transpilation and Qubit Layout

The PennyLane-Qiskit device handles transpilation automatically. To control the process, pass a `transpile_options` dictionary with `optimization_level` (0–3, where 3 performs the most aggressive gate cancellation and routing), and an `initial_layout` specifying which physical qubits to use (e.g., 4 qubits forming a line or star in the coupling map with minimal CNOT error). Use the `select_optimal_qubits` strategy from the QAOA report.

#### 1.2 Pauli Grouping and Measurement Circuits

The H₂ Hamiltonian has 15 Pauli terms. Not all can be measured simultaneously. Commuting Pauli terms can be grouped:

**Commuting clique partitioning:** Group Paulis that mutually commute (can be simultaneously diagonalised). For H₂:

- Group 1 (all Z-type): $\{I, Z_0, Z_1, Z_2, Z_3, Z_0Z_1, Z_0Z_2, Z_0Z_3, Z_1Z_2, Z_1Z_3, Z_2Z_3\}$ — all commute, measured in the computational basis.
- Group 2: $\{X_0X_1Y_2Y_3, Y_0Y_1X_2X_3\}$ — commute.
- Group 3: $\{X_0Y_1Y_2X_3, Y_0X_1X_2Y_3\}$ — commute.

This requires only 3 distinct measurement circuits instead of 15 separate ones.

PennyLane performs this grouping automatically when using `qml.expval(H)` with a hardware device. The `grouping_type` parameter (e.g., `"qwc"` for qubit-wise commuting) on `qml.Hamiltonian` controls the strategy.

#### 1.3 Shot Budget for Chemical Accuracy

For the H₂ Hamiltonian with $\sum_k |c_k|^2 \approx 0.3$ Hartree$^2$:

$$S \geq \frac{\sum_k |c_k|^2}{\epsilon^2} = \frac{0.3}{(1.6 \times 10^{-3})^2} \approx 117{,}000 \text{ shots}$$

With Pauli grouping (3 groups), distribute shots proportionally to the variance contribution of each group. In practice, 100,000–200,000 total shots (across all groups) suffices for chemical accuracy on H₂.

### 2. Error Mitigation for Chemistry

#### 2.1 ZNE for VQE Energy

ZNE is particularly well-suited for VQE because:
- The output is a scalar (energy), and ZNE extrapolates scalar values.
- The noise in VQE circuits is dominated by incoherent errors (after Pauli twirling), which scale approximately linearly or polynomially with noise factor.

Apply ZNE exactly as described in the QAOA report (global unitary folding at scale factors 1, 3, 5 with Richardson extrapolation). For the H₂ UCCSD circuit (depth ~16 CNOTs):
- $\lambda = 1$: 16 CNOTs
- $\lambda = 3$: 48 CNOTs 
- $\lambda = 5$: 80 CNOTs

At 80 CNOTs on ibm_brisbane ($\sim 0.5\%$ CNOT error), the total error probability is $\sim 33\%$, which is marginal but still within the regime where ZNE works. For deeper circuits (LiH, H₂O), ZNE may not suffice.

#### 2.2 Symmetry Verification

Chemistry problems have symmetries that the exact ground state obeys:
- **Particle number:** $\langle\hat{N}\rangle = N_e$
- **Spin:** $\langle\hat{S}_z\rangle = 0$ (for singlet states), $\langle\hat{S}^2\rangle = 0$

Post-selection: discard measurement outcomes that violate these symmetries. This reduces the effective shot count but removes unphysical contributions. For H₂, filter the raw bitstring samples to keep only those whose Hamming weight equals $N_e$ (the electron count), then recompute the energy expectation from the filtered sample set only.

#### 2.3 Expected Performance

| Configuration | $|E_{\text{VQE}} - E_{\text{FCI}}|$ | Chemical accuracy? |
|:--|:-:|:-:|
| H₂, noiseless sim | $< 10^{-6}$ Ha | Yes |
| H₂, noisy sim (1% depol) | $\sim 5$ mHa | No |
| H₂, ibm_brisbane (raw) | $\sim 10$–$30$ mHa | No |
| H₂, ibm_brisbane + ZNE | $\sim 2$–$5$ mHa | Marginal |
| H₂, ibm_brisbane + ZNE + symmetry | $\sim 1$–$3$ mHa | Possible |

### 3. Comparison Checklist

| Metric | Sim (noiseless) | Sim (noisy) | Hardware | HW + ZNE |
|:-------|:---:|:---:|:---:|:---:|
| $E_{\text{VQE}}$ at $R_e$ | ✓ | ✓ | ✓ | ✓ |
| $|E_{\text{VQE}} - E_{\text{FCI}}|$ | ✓ | ✓ | ✓ | ✓ |
| Correlation recovery ratio | ✓ | ✓ | ✓ | ✓ |
| PES shape (qualitative) | ✓ | ✓ | ✓ | ✓ |
| Dissociation energy $D_e$ | ✓ | ✓ | ✓ | ✓ |
| Equilibrium $R_e$ | ✓ | ✓ | ✓ | ✓ |
| VQE iterations to converge | ✓ | ✓ | — | — |
| Total circuit evaluations | ✓ | ✓ | ✓ | ✓ |

**Plots to produce:**
1. PES: $E(R)$ for HF, CISD, CCSD(T), VQE, FCI on same axes.
2. Error: $|E_{\text{method}} - E_{\text{FCI}}|$ vs $R$ for each method.
3. Correlation recovery ratio vs $R$.
4. VQE convergence: energy vs iteration at $R_e$ and at $R = 3.0$ Å.

---

## Part 6 — Benchmarking & Validation

### 1. Validation Against PySCF FCI

To validate VQE energies, a helper function builds a PySCF `Mole` from the `MoleculeSpec`, runs RHF + FCI to obtain the exact ground-state energy $E_{\text{FCI}}$, and computes the absolute error $\Delta E = |E_{\text{VQE}} - E_{\text{FCI}}|$. The error is reported in Hartree, milliHartree ($\times 1000$), and kcal/mol ($\times 627.509$). Chemical accuracy ($\Delta E < 1.6 \times 10^{-3}$ Ha) is flagged with a boolean.

*(Full implementation: `report2_scripts/benchmarker.py` and `report2_scripts/classical_hf.py`)*

### 2. Metrics to Measure and Plot

**Correlation energy recovered:**

$$r_{\text{corr}} = \frac{E_{\text{VQE}} - E_{\text{HF}}}{E_{\text{FCI}} - E_{\text{HF}}}$$

A value of $r_{\text{corr}} = 1.0$ means VQE recovers all correlation energy.

**State fidelity:** If the VQE state $|\psi_{\text{VQE}}\rangle$ and exact ground state $|\psi_0\rangle$ are both available (in simulation):

$$F = |\langle\psi_{\text{VQE}}|\psi_0\rangle|^2$$

In PennyLane, extract the statevector via `qml.state()` and compute overlap with the FCI eigenvector.

**Circuit depth vs accuracy:**

| UCCSD layers | CNOT depth (H₂) | $|E - E_{\text{FCI}}|$ |
|:--:|:--:|:--:|
| 1 Trotter step | ~16 | $< 10^{-6}$ |
| 2 Trotter steps | ~32 | $< 10^{-8}$ |

For H₂, a single Trotter step of UCCSD is exact (the ansatz is expressible). For LiH and H₂O, multiple Trotter steps improve accuracy.

### 3. Results Table Template

| $R$ (Å) | $E_{\text{HF}}$ | $E_{\text{CISD}}$ | $E_{\text{CCSD(T)}}$ | $E_{\text{VQE}}$ | $E_{\text{FCI}}$ | $|E_{\text{VQE}} - E_{\text{FCI}}|$ |
|:------:|:----------:|:-----------:|:-------------:|:---------:|:---------:|:------------------:|
| 0.50 | | | | | | |
| 0.75 | | | | | | |
| 1.00 | | | | | | |
| 1.50 | | | | | | |
| 2.50 | | | | | | |

---

## Part 7 — Research Contribution Angle

### 1. Three Directions Toward a Publishable Result

**Direction 1: Adaptive Ansatz PES — ADAPT-VQE vs UCCSD Across Dissociation**

Implement ADAPT-VQE (Grimsley et al., 2019), which builds the ansatz one operator at a time by selecting the operator with the largest gradient at each step. Compare the PES from ADAPT-VQE vs fixed UCCSD across the full dissociation curve for H₂, LiH, and BeH₂. Key questions: Does ADAPT-VQE use fewer parameters? Is the PES smoother? Does it handle the multi-reference region (near dissociation) better? Benchmark circuit depth and parameter count at each bond length.

**Direction 2: Noise-Aware VQE for Robust PES**

Study how noisy VQE distorts the PES and propose corrections. At each bond length, run VQE under varying noise strengths and measure: (a) equilibrium bond length shift, (b) dissociation energy error, (c) PES curvature distortion (affects vibrational frequencies). Propose a noise-aware cost function $\tilde{E}(\theta; \lambda) = E_{\text{noisy}}(\theta) - f(\lambda)$ where $f(\lambda)$ is a noise correction learned from ZNE. Show that this gives a more accurate PES than raw noisy VQE.

**Direction 3: SHARC-VQE Partitioned Hamiltonian Benchmarking**

Implement the SHARC (Singlet-Triplet Hamiltonian Arithmetic) partitioning from Ralli et al. (2023), which splits the molecular Hamiltonian into fragments that can be solved independently on smaller quantum devices. Benchmark the partitioned approach vs monolithic UCCSD for H₂O and small organic molecules (methane, ethylene). Measure the accuracy-vs-resource tradeoff and identify the optimal partitioning strategy for NISQ hardware.

### 2. Most Relevant 2024–2025 Papers

1. **Grimsley et al. (2024).** "ADAPT-VQE is insensitive to rough parameter landscapes and barren plateaus." *npj Quantum Information*. — Robustness of adaptive ansätze.

2. **Motta et al. (2024).** "Bridging quantum computing and quantum chemistry: Applications of VQE to molecular systems." *WIREs Computational Molecular Science*. — Comprehensive VQE review.

3. **Ralli et al. (2023–2024).** "SHARC-VQE: Simplified Hamiltonian Approach with Resource-efficient Circuits." *Physical Review Research*. — Hamiltonian partitioning for VQE.

4. **Gonthier et al. (2024).** "Identifying challenges towards practical quantum advantage through resource estimation." *Communications Physics*. — What it takes for quantum advantage in chemistry.

5. **Lee et al. (2024).** "Is there evidence for exponential quantum advantage in quantum chemistry?" *Nature Communications*. — Critical analysis of quantum chemistry advantage claims.

6. **Cerezo et al. (2024).** "Variational quantum algorithms: fundamental limits, robustness, and practical implementations." *Reviews of Modern Physics*. — Definitive VQA review.

7. **Tilly et al. (2022, updated 2024).** "The Variational Quantum Eigensolver: A review of methods and best practices." *Physics Reports*. — Canonical VQE reference.

### 3. Tractable Molecules on Today's Hardware

| Molecule | Qubits (STO-3G) | Qubits (6-31G*) | Active Space | Tractable? |
|:---------|:---:|:---:|:---:|:---:|
| H₂ | 4 (2 tapered) | 8 | Full | Yes (demonstrated) |
| LiH | 12 (10 tapered) | 22 | (2e, 5o) → 10 qubits | Yes with active space |
| BeH₂ | 14 | 26 | (4e, 6o) → 12 qubits | Marginal |
| H₂O | 14 (12 tapered) | 38 | (4e, 4o) → 8 qubits | Yes with active space |
| NH₃ | 16 | 40 | (6e, 6o) → 12 qubits | Marginal |
| N₂ | 20 | 48 | (6e, 6o) → 12 qubits | Marginal (strong corr.) |

**Why active space is essential:** The full orbital space is far too large for current hardware. Active space methods freeze core electrons and discard high-energy virtual orbitals, reducing the qubit count to a tractable number. The active space must be carefully chosen to include the orbitals involved in bond-breaking (the valence space).

**N₂ is the holy grail:** The triple bond in N₂ involves strong static correlation across 6 electrons in 6 orbitals. Accurately computing the N₂ PES is a long-standing challenge for both classical and quantum methods. Demonstrating VQE on N₂ with an active space of (6e, 6o) on 12 qubits would be a significant milestone.
