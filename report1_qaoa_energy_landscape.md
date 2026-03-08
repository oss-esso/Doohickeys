# QAOA Energy Landscape Explorer: Full Research & Implementation Report

---

## Part 1 — Theoretical Foundations

### 1. Combinatorial Optimization & the Ising Model

#### 1.1 Derivation of the Ising Hamiltonian from First Principles

The classical Ising model originates in statistical mechanics as a model of interacting magnetic spins on a lattice. Consider $N$ binary spin variables $s_i \in \{-1, +1\}$, each sitting on a vertex $i$ of a graph $G = (V, E)$. The most general pairwise Ising Hamiltonian with external fields is:

$$H_{\text{Ising}}(\mathbf{s}) = -\sum_{(i,j) \in E} J_{ij}\, s_i s_j - \sum_{i \in V} h_i\, s_i$$

where $J_{ij}$ are coupling constants (positive = ferromagnetic, negative = anti-ferromagnetic) and $h_i$ are local field strengths. The ground state of this Hamiltonian corresponds to the spin configuration $\mathbf{s}^*$ that minimises $H_{\text{Ising}}$.

To promote this to a quantum mechanical operator, replace each classical spin $s_i$ with the Pauli-$Z$ operator $\hat{Z}_i$ acting on the $i$-th qubit:

$$\hat{H}_{\text{Ising}} = -\sum_{(i,j) \in E} J_{ij}\, \hat{Z}_i \hat{Z}_j - \sum_{i \in V} h_i\, \hat{Z}_i$$

This is diagonal in the computational basis $\{|0\rangle, |1\rangle\}^{\otimes N}$ since each $\hat{Z}_i$ has eigenvalues $\pm 1$ on $|0\rangle$ and $|1\rangle$ respectively. The eigenvalues of $\hat{H}_{\text{Ising}}$ on a computational basis state $|z_1 z_2 \cdots z_N\rangle$ reproduce the classical Ising energy $H_{\text{Ising}}(s_1, \ldots, s_N)$ under the identification $s_i = 1 - 2z_i$ (i.e., $|0\rangle \mapsto s=+1$, $|1\rangle \mapsto s=-1$).

#### 1.2 QUBO ↔ Ising Equivalence

A Quadratic Unconstrained Binary Optimization (QUBO) problem is defined over binary variables $x_i \in \{0, 1\}$:

$$H_{\text{QUBO}}(\mathbf{x}) = \sum_{i} Q_{ii}\, x_i + \sum_{i < j} Q_{ij}\, x_i x_j$$

where $Q$ is a real symmetric (or upper-triangular) matrix. To convert QUBO → Ising, substitute $x_i = \frac{1 - s_i}{2}$:

$$x_i x_j = \frac{(1-s_i)(1-s_j)}{4} = \frac{1 - s_i - s_j + s_i s_j}{4}$$

$$x_i = \frac{1 - s_i}{2}$$

Substituting into $H_{\text{QUBO}}$:

$$H_{\text{QUBO}} = \sum_i Q_{ii} \frac{1 - s_i}{2} + \sum_{i<j} Q_{ij} \frac{1 - s_i - s_j + s_i s_j}{4}$$

Collecting terms:

$$H_{\text{Ising}}(\mathbf{s}) = \sum_{i<j} \frac{Q_{ij}}{4}\, s_i s_j - \sum_i \left(\frac{Q_{ii}}{2} + \sum_{j \neq i} \frac{Q_{ij}}{4}\right) s_i + \text{const}$$

Hence the Ising couplings and fields are:

$$J_{ij} = \frac{Q_{ij}}{4}, \qquad h_i = \frac{Q_{ii}}{2} + \sum_{j \neq i} \frac{Q_{ij}}{4}$$

(with appropriate sign conventions matching the $-J_{ij}$ in the Ising Hamiltonian).

#### 1.3 Worked Example: MaxCut on a 4-Node Graph

Consider the cycle graph $C_4$ on vertices $\{0, 1, 2, 3\}$ with edges $E = \{(0,1), (1,2), (2,3), (3,0)\}$, all with unit weight.

**MaxCut objective:** Maximise the number of edges crossing the partition. For edge $(i,j)$, the contribution is $\frac{1 - s_i s_j}{2}$ (equals 1 if $s_i \neq s_j$, 0 otherwise). Total:

$$C(\mathbf{s}) = \sum_{(i,j) \in E} \frac{1 - s_i s_j}{2} = \frac{|E|}{2} - \frac{1}{2}\sum_{(i,j) \in E} s_i s_j$$

To maximise $C$, we minimise $\sum_{(i,j)} s_i s_j$. The cost Hamiltonian (to be minimised) in qubit form is:

$$\hat{H}_C = \frac{1}{2}\sum_{(i,j) \in E} (\hat{I} - \hat{Z}_i \hat{Z}_j) = \frac{|E|}{2}\hat{I} - \frac{1}{2}\sum_{(i,j) \in E} \hat{Z}_i \hat{Z}_j$$

For $C_4$:

$$\hat{H}_C = 2\hat{I} - \frac{1}{2}(\hat{Z}_0\hat{Z}_1 + \hat{Z}_1\hat{Z}_2 + \hat{Z}_2\hat{Z}_3 + \hat{Z}_3\hat{Z}_0)$$

Since we maximise $\langle \hat{H}_C \rangle$, equivalently we minimise $-\hat{H}_C$. The maximum cut on $C_4$ is 4 (e.g., $|0101\rangle$ or $|1010\rangle$, putting alternating vertices in different sets), giving $C_{\max} = 4$.

**Explicit Ising form (to be maximised):** $J_{ij} = -\frac{1}{2}$ for all edges, $h_i = 0$. The ground state of $\hat{H}_{\text{cost}} = -\hat{H}_C$ has eigenvalue $-4$.

#### 1.4 The Energy Landscape: Geometric Definition and Significance

Given a parameterised quantum state $|\psi(\boldsymbol{\theta})\rangle$ and a cost Hamiltonian $\hat{H}_C$, the **energy landscape** is the scalar field:

$$\mathcal{L}(\boldsymbol{\theta}) = \langle \psi(\boldsymbol{\theta}) | \hat{H}_C | \psi(\boldsymbol{\theta}) \rangle$$

mapping the parameter space $\boldsymbol{\theta} \in \mathbb{R}^d$ to the real line. For QAOA with depth $p$, $\boldsymbol{\theta} = (\gamma_1, \ldots, \gamma_p, \beta_1, \ldots, \beta_p) \in \mathbb{R}^{2p}$.

When $p = 1$, $\mathcal{L}(\gamma, \beta)$ is a two-dimensional surface embedded in $\mathbb{R}^3$, directly visualisable. Geometrically, this surface encodes:

- **Global minima**: the optimal QAOA parameters yielding the best approximation ratio.
- **Local minima**: suboptimal parameter traps that gradient-based optimisers may converge to.
- **Saddle points**: critical points where the Hessian has mixed signature — these can slow convergence.
- **Flat regions (barren plateaus)**: zones where the gradient vanishes exponentially, making navigation impossible.

The topology and curvature of this landscape directly determine the trainability of the variational algorithm, the success of classical optimisers, and the quality of the approximate solution. Understanding landscape structure is therefore a prerequisite for QAOA performance analysis.

---

### 2. QAOA: Full Derivation

#### 2.1 From Quantum Approximate Optimization to the QAOA Circuit

The Quantum Approximate Optimization Algorithm (Farhi, Goldstone, Gutmann, 2014) constructs a sequence of parameterised unitaries designed to prepare a state with high overlap on the ground state of a cost Hamiltonian $\hat{H}_C$.

**Setup:** Given a combinatorial optimization problem encoded in a cost Hamiltonian $\hat{H}_C$ (diagonal in the computational basis), define:

- **Cost unitary**: $\hat{U}_C(\gamma) = e^{-i\gamma \hat{H}_C}$
- **Mixer unitary**: $\hat{U}_B(\beta) = e^{-i\beta \hat{H}_B}$ where $\hat{H}_B = \sum_{i=1}^{N} \hat{X}_i$ is the transverse-field mixer.

**Initial state:** The uniform superposition $|+\rangle^{\otimes N} = \hat{H}^{\otimes N}|0\rangle^{\otimes N}$, which is the ground state of $\hat{H}_B$.

**QAOA state at depth $p$:**

$$|\boldsymbol{\gamma}, \boldsymbol{\beta}\rangle = \hat{U}_B(\beta_p)\hat{U}_C(\gamma_p) \cdots \hat{U}_B(\beta_1)\hat{U}_C(\gamma_1)|+\rangle^{\otimes N}$$

Equivalently:

$$|\boldsymbol{\gamma}, \boldsymbol{\beta}\rangle = \prod_{k=1}^{p} \bigl[\hat{U}_B(\beta_k)\,\hat{U}_C(\gamma_k)\bigr] |+\rangle^{\otimes N}$$

where the product is time-ordered (rightmost applied first).

#### 2.2 Cost Hamiltonian $\hat{H}_C$ and Mixer Hamiltonian $\hat{H}_B$

For MaxCut on a graph $G = (V, E)$ with edge weights $w_{ij}$:

$$\hat{H}_C = \sum_{(i,j) \in E} \frac{w_{ij}}{2}(\hat{I} - \hat{Z}_i \hat{Z}_j)$$

The mixer Hamiltonian drives transitions between computational basis states:

$$\hat{H}_B = \sum_{i=1}^{N} \hat{X}_i$$

The choice of mixer is not unique. The $X$-mixer preserves no symmetries and explores the full Hilbert space. For constrained problems, one may use XY-mixers, Grover mixers, or other symmetry-preserving alternatives (Hadfield et al., 2019).

#### 2.3 Explicit Gate Decomposition of the Full Unitary

**Cost unitary decomposition:**

Since $\hat{H}_C = \sum_{(i,j)} \frac{w_{ij}}{2}(\hat{I} - \hat{Z}_i\hat{Z}_j)$, and all $\hat{Z}_i\hat{Z}_j$ terms commute mutually (they are all diagonal), we can decompose:

$$\hat{U}_C(\gamma) = e^{-i\gamma \hat{H}_C} = \prod_{(i,j) \in E} e^{-i\gamma \frac{w_{ij}}{2}(\hat{I} - \hat{Z}_i\hat{Z}_j)}$$

Each factor:

$$e^{-i\gamma \frac{w_{ij}}{2}(\hat{I} - \hat{Z}_i\hat{Z}_j)} = e^{-i\gamma w_{ij}/2} \cdot e^{i\gamma \frac{w_{ij}}{2} \hat{Z}_i\hat{Z}_j}$$

The global phase $e^{-i\gamma w_{ij}/2}$ is unphysical. The $ZZ$-rotation is implemented as:

$$e^{i\gamma \frac{w_{ij}}{2} \hat{Z}_i \hat{Z}_j} = \text{CNOT}_{i,j} \cdot R_Z\bigl(-\gamma w_{ij}\bigr)_j \cdot \text{CNOT}_{i,j}$$

where $R_Z(\theta) = e^{-i\theta \hat{Z}/2}$. Alternatively, using the native $R_{ZZ}$ gate:

$$R_{ZZ}(\theta) = e^{-i\frac{\theta}{2}\hat{Z}\otimes\hat{Z}}$$

so the cost layer is a sequence of $R_{ZZ}(-\gamma w_{ij})$ gates on all edge pairs.

**Mixer unitary decomposition:**

Since all $\hat{X}_i$ operators commute when acting on distinct qubits:

$$\hat{U}_B(\beta) = e^{-i\beta \sum_i \hat{X}_i} = \prod_{i=1}^{N} e^{-i\beta \hat{X}_i} = \prod_{i=1}^{N} R_X(2\beta)_i$$

where $R_X(\theta) = e^{-i\theta \hat{X}/2}$.

**Full circuit:** For depth $p$, the circuit is:

1. Apply $\hat{H}^{\otimes N}$ to $|0\rangle^{\otimes N}$
2. For $k = 1, \ldots, p$:
   - Apply $R_{ZZ}(-\gamma_k w_{ij})$ for each edge $(i,j) \in E$
   - Apply $R_X(2\beta_k)$ to each qubit $i \in V$
3. Measure in the computational basis

#### 2.4 The $p \to \infty$ Limit and Adiabatic Connection

QAOA is formally equivalent to a digitised (Trotterised) quantum annealing schedule. Consider the time-dependent Hamiltonian:

$$\hat{H}(t) = \left(1 - \frac{t}{T}\right)\hat{H}_B + \frac{t}{T}\hat{H}_C, \qquad t \in [0, T]$$

The adiabatic theorem guarantees that if $T$ is large enough (specifically $T \sim O(1/\Delta_{\min}^2)$ where $\Delta_{\min}$ is the minimum spectral gap), the system tracks the instantaneous ground state and ends in the ground state of $\hat{H}_C$.

Trotterise this evolution into $p$ steps with $\Delta t = T/p$:

$$\hat{U} \approx \prod_{k=1}^{p} e^{-i\Delta t (1 - t_k/T)\hat{H}_B} \, e^{-i\Delta t (t_k/T)\hat{H}_C}$$

This is exactly the QAOA unitary with:

$$\beta_k = \Delta t\left(1 - \frac{k}{p}\right), \qquad \gamma_k = \Delta t \cdot \frac{k}{p}$$

So for fixed $T$, as $p \to \infty$, QAOA reproduces the adiabatic evolution exactly (up to Trotter error $O(\Delta t^2)$). Since QAOA optimises over all $(\boldsymbol{\gamma}, \boldsymbol{\beta})$ and is not restricted to linear annealing schedules, it is strictly more powerful than Trotterised adiabatic evolution: it can exploit shortcuts to adiabaticity and counterdiabatic terms implicitly.

#### 2.5 Analytical Expectation Value for $p=1$ on $C_4$

For the cycle graph $C_4$ with $\hat{H}_C = \sum_{(i,j) \in E} \frac{1}{2}(\hat{I} - \hat{Z}_i \hat{Z}_j)$ (4 edges, 4 qubits), the $p=1$ QAOA expectation value can be computed analytically.

Start with $|+\rangle^{\otimes 4}$. Apply $\hat{U}_C(\gamma) = \prod_{(i,j)} e^{i\gamma/2 \hat{Z}_i\hat{Z}_j}$ (dropping global phases), then $\hat{U}_B(\beta) = \prod_i R_X(2\beta)_i$.

By the linearity of $\langle \hat{H}_C \rangle$, we need:

$$\langle \hat{H}_C \rangle = \sum_{(i,j) \in E} \frac{1}{2}\bigl(1 - \langle \hat{Z}_i \hat{Z}_j \rangle\bigr)$$

By the symmetry of the cycle graph, all edge expectation values are identical. Define $f(\gamma, \beta) = \langle \hat{Z}_0 \hat{Z}_1 \rangle$. Then $\langle \hat{H}_C \rangle = 4 \cdot \frac{1}{2}(1 - f) = 2(1 - f)$.

For the cycle on 4 qubits (ring topology), the $p=1$ expectation value is:

$$f(\gamma, \beta) = \langle +|^{\otimes 4}\, \hat{U}_C^\dagger(\gamma)\, \hat{U}_B^\dagger(\beta)\, \hat{Z}_0 \hat{Z}_1\, \hat{U}_B(\beta)\, \hat{U}_C(\gamma)\, |+\rangle^{\otimes 4}$$

Using the fact that $R_X(2\beta)^\dagger \hat{Z} R_X(2\beta) = \cos(2\beta)\hat{Z} - \sin(2\beta)\hat{Y}$ and propagating through the $ZZ$ rotations, one obtains (after considerable algebra — see Farhi et al. Appendix or Wang et al., PRX 2018):

$$f(\gamma,\beta) = \cos^2(2\beta)\bigl[\cos(2\gamma) - \sin^2(2\gamma)\cos(2\gamma)\bigr] + \frac{1}{2}\sin^2(2\beta)\sin^2(2\gamma)\bigl[1 + \cos^2(2\gamma)\bigr]\sin(4\beta)(\ldots)$$

The full closed-form for a $d$-regular graph at $p=1$ is given by the Farhi-Goldstone-Gutmann formula. For the $d$-regular case with edge $(i,j)$ where vertex $i$ has $d_i - 1$ other neighbours and vertex $j$ has $d_j - 1$ other neighbours:

$$\langle \hat{Z}_i \hat{Z}_j \rangle = \frac{1}{2}\sin(4\beta)\sin(2\gamma)\bigl[\cos^{d_i-1}(2\gamma) + \cos^{d_j-1}(2\gamma)\bigr] - \frac{1}{4}\sin^2(2\beta)\sin^2(2\gamma)\Bigl[\cos^{n_{ij}}(2\gamma)\Bigr]$$

where $n_{ij} = d_i + d_j - 2 - |N(i) \cap N(j)|$ accounts for shared neighbours. For the 2-regular cycle $C_4$, $d_i = d_j = 2$, and $|N(0) \cap N(1)| = 0$, so $n_{01} = 2$. Substituting:

$$\langle \hat{Z}_0\hat{Z}_1\rangle_{C_4} = \frac{1}{2}\sin(4\beta)\sin(2\gamma)\bigl[2\cos(2\gamma)\bigr] - \frac{1}{4}\sin^2(2\beta)\sin^2(2\gamma)\cos^2(2\gamma)$$

$$= \sin(4\beta)\sin(2\gamma)\cos(2\gamma) - \frac{1}{4}\sin^2(2\beta)\sin^2(2\gamma)\cos^2(2\gamma)$$

Then:

$$\langle \hat{H}_C \rangle = 2\bigl(1 - f(\gamma,\beta)\bigr) = 2 - 2\Bigl[\sin(4\beta)\sin(4\gamma)/2 - \frac{1}{4}\sin^2(2\beta)\sin^2(4\gamma)/4\Bigr]$$

The maximum of $\langle \hat{H}_C \rangle$ over $(\gamma, \beta)$ gives the QAOA $p=1$ approximation ratio $r_1 = \langle \hat{H}_C \rangle_{\max} / C_{\max}$, where $C_{\max} = 4$.

Numerically, for $C_4$: $r_1 \approx 0.75$ (3 out of 4 edges), achieved near $\gamma^* \approx \pi/8$, $\beta^* \approx \pi/8$.

---

### 3. The $(\beta, \gamma)$ Parameter Landscape

#### 3.1 Periodicity

The QAOA landscape is periodic in both parameters. The periods depend on the Hamiltonian structure:

**$\gamma$-periodicity:** The cost unitary $e^{-i\gamma \hat{H}_C}$ is periodic in $\gamma$ with period determined by the spectral structure of $\hat{H}_C$. For MaxCut with integer weights, the eigenvalues of $\hat{H}_C$ are integers (or half-integers), and the cost unitary is $2\pi$-periodic. More precisely, if every eigenvalue difference of $\hat{H}_C$ is an integer, the period in $\gamma$ is $2\pi$. For the standard MaxCut formulation $\hat{H}_C = \sum_{(i,j)} \frac{1}{2}(I - Z_iZ_j)$, the eigenvalue differences are integers, so:

$$\text{Period in } \gamma = 2\pi$$

However, exploiting the parity symmetry $\hat{H}_C(-\gamma) = \hat{H}_C(\gamma)$ and the fact that the landscape has a reflection symmetry, it suffices to search $\gamma \in [0, \pi]$.

**$\beta$-periodicity:** The mixer unitary $e^{-i\beta \sum_i X_i}$: since $X_i$ has eigenvalues $\pm 1$, the $R_X(2\beta)$ gate is $2\pi$-periodic in $2\beta$, i.e., $\pi$-periodic in $\beta$:

$$\text{Period in } \beta = \pi$$

Again, the symmetry $\langle H_C \rangle(\beta) = \langle H_C \rangle(-\beta)$ allows restriction to $\beta \in [0, \pi/2]$.

**Combined search space for $p=1$:** $(\gamma, \beta) \in [0, \pi] \times [0, \pi/2]$ (or further reducible depending on graph symmetries).

#### 3.2 Landscape Structure

**Concentration of parameters (Brandao et al., 2018; Wurtz & Love, 2021):**

For $d$-regular graphs in the large $N$ limit, the optimal QAOA parameters $(\boldsymbol{\gamma}^*, \boldsymbol{\beta}^*)$ at fixed $p$ concentrate — they become independent of the specific graph instance and depend only on $d$ and $p$. This means:

- One can precompute "universal" optimal parameters on small instances and transfer them to larger graphs of the same regularity class.
- The variance of the optimal energy across random $d$-regular graphs decays as $O(1/N)$.

**Symmetry properties:**

1. **Inversion symmetry:** $\langle H_C \rangle(\boldsymbol{\gamma}, \boldsymbol{\beta}) = \langle H_C \rangle(-\boldsymbol{\gamma}, -\boldsymbol{\beta})$ — the landscape is symmetric under simultaneous negation.
2. **Time-reversal:** For real Hamiltonians, $\langle H_C \rangle(\boldsymbol{\gamma}, \boldsymbol{\beta}) = \langle H_C \rangle(\boldsymbol{\gamma}^R, \boldsymbol{\beta}^R)$ where $R$ denotes reversal of parameter ordering.
3. **Graph automorphisms** induce symmetries in the landscape.

**Good initial parameter strategies:**

- **INTERP (Zhou et al., 2020):** Optimise at depth $p$, then linearly interpolate the parameters to initialise depth $p+1$. Specifically, for $p \to p+1$:

$$\gamma_k^{(p+1)} = \frac{k-1}{p}\gamma_{k-1}^{(p)} + \frac{p+1-k}{p}\gamma_k^{(p)}, \quad k = 1, \ldots, p+1$$

(similarly for $\beta$). This exploits the smoothness of optimal parameter schedules.

- **FOURIER (Zhou et al., 2020):** Parameterise the QAOA angles in a Fourier basis:

$$\gamma_k = \sum_{q=1}^{Q} u_q \sin\!\left(\frac{(2k-1)q\pi}{2p}\right), \qquad \beta_k = \sum_{q=1}^{Q} v_q \cos\!\left(\frac{(2k-1)q\pi}{2p}\right)$$

where $Q \leq p$ is the number of Fourier components. Optimisation over $(u_1, \ldots, u_Q, v_1, \ldots, v_Q)$ is in a reduced $2Q$-dimensional space instead of $2p$, avoiding local minima.

- **TQA (Trotterised Quantum Annealing):** Use the linear schedule from the adiabatic connection (Section 2.4) as initialisation.

#### 3.3 Barren Plateaus

**Formal definition:** A variational circuit exhibits a barren plateau if the variance of the cost function partial derivatives vanishes exponentially with system size:

$$\text{Var}\!\left[\frac{\partial \langle \hat{H}_C \rangle}{\partial \theta_k}\right] \leq F(N)$$

where $F(N) \in O(1/b^N)$ for some $b > 1$.

For QAOA specifically, the situation is nuanced:

1. **Hardware-efficient random ansätze** (McClean et al., 2018): For circuits that form approximate 2-designs, the gradient variance decays as $O(2^{-N})$, guaranteeing barren plateaus.

2. **QAOA with fixed $p$:** QAOA at constant depth $p$ does NOT suffer from barren plateaus in the same way because the circuit is highly structured (alternating cost and mixer unitaries) and shallow. The gradient variance scales polynomially (often as $\Theta(1/\text{poly}(N))$), which is trainable.

3. **QAOA with $p = O(\text{poly}(N))$:** As $p$ grows polynomially with $N$, the landscape can develop barren plateau features depending on the cost function locality. For local cost functions (each term acts on $O(1)$ qubits), Cerezo et al. (2021) showed that shallow circuits have well-behaved gradients, while deep circuits suffer barren plateaus regardless.

**Practical consequence:** For the typical QAOA regime ($p \leq 10$, $N \leq 30$), barren plateaus are not the primary obstacle. Instead, the proliferation of local minima and the non-convexity of the landscape are the main challenges.

#### 3.4 Phase Transition Behaviour with Local Field Strength

The 2025 result in *Quantum* (Boulebnane & Montanaro, 2025; also Kremenetski et al., 2025) studies the QAOA landscape on anti-ferromagnetic Ising models with longitudinal fields:

$$\hat{H}_C = -\sum_{(i,j)} J_{ij}\hat{Z}_i\hat{Z}_j + h\sum_i \hat{Z}_i$$

where $J_{ij} < 0$ (anti-ferromagnetic). Three distinct regimes are identified as $h$ varies:

1. **Trivial regime** ($|h| \gg |J|$): The ground state is the fully polarised state (all spins aligned with the field). QAOA finds this trivially even at $p=1$ — the landscape has a single dominant minimum.

2. **Anti-ferromagnetic regime** ($|h| \ll |J|$): The ground state exhibits Néel-type anti-ferromagnetic ordering. The landscape at low $p$ is highly structured with multiple degenerate minima corresponding to the broken symmetry. The QAOA performance depends strongly on the graph geometry (bipartite vs. frustrated).

3. **Critical regime** ($|h| \sim |J|$, near the quantum phase transition): The spectral gap $\Delta$ closes polynomially or logarithmically with $N$, the landscape becomes complex with many nearly degenerate minima, and the approximation ratio achievable at fixed $p$ drops sharply. This is where QAOA most struggles — and where the landscape structure is richest from a visualisation perspective.

The key finding is that QAOA performance undergoes a sharp transition at the critical field strength, and the landscape topology (number of local minima, gradient magnitudes, Hessian structure) changes qualitatively across the phase boundary.

---

### 4. Classical Optimization of QAOA Parameters

#### 4.1 Parameter-Shift Rule: Full Derivation

For a parameterised unitary of the form $U(\theta) = e^{-i\theta G/2}$ where $G$ is a generator with eigenvalues $\pm r$ (typically $r = 1$ for Pauli generators), the parameter-shift rule states:

$$\frac{\partial}{\partial \theta}\langle \psi(\theta) | \hat{O} | \psi(\theta)\rangle = \frac{r}{2}\Bigl[\langle \hat{O}\rangle_{\theta + \pi/(2r)} - \langle \hat{O}\rangle_{\theta - \pi/(2r)}\Bigr]$$

**Derivation:** Let $|\psi(\theta)\rangle = U(\theta)|\psi_0\rangle$ with $U(\theta) = e^{-i\theta G/2}$. Decompose $|\psi_0\rangle$ in the eigenbasis of $G$: $G = r|+\rangle\langle +| - r|-\rangle\langle -|$.

$$|\psi_0\rangle = c_+|+\rangle + c_-|-\rangle$$

$$|\psi(\theta)\rangle = c_+ e^{-ir\theta/2}|+\rangle + c_- e^{ir\theta/2}|-\rangle$$

Let $f(\theta) = \langle\psi(\theta)|\hat{O}|\psi(\theta)\rangle$. Computing:

$$f(\theta) = |c_+|^2 O_{++} + |c_-|^2 O_{--} + c_+^* c_- e^{ir\theta} O_{+-} + c_- ^* c_+ e^{-ir\theta} O_{-+}$$

where $O_{ab} = \langle a|\hat{O}|b\rangle$. This is a sinusoidal function with frequency $r$:

$$f(\theta) = A + B\cos(r\theta) + C\sin(r\theta)$$

Therefore:

$$f'(\theta) = -Br\sin(r\theta) + Cr\cos(r\theta) = r\bigl[f(\theta + \pi/(2r)) - f(\theta - \pi/(2r))\bigr] / 2$$

which follows from the identity:

$$\frac{f(\theta + s) - f(\theta - s)}{2} = B\bigl[\cos r(\theta+s) - \cos r(\theta-s)\bigr]/2 + C\bigl[\sin r(\theta + s) - \sin r(\theta-s)\bigr]/2$$

Setting $s = \pi/(2r)$:

$$= -B\sin(r\theta)\sin(\pi/2) + C\cos(r\theta)\sin(\pi/2) = -B\sin(r\theta) + C\cos(r\theta) = f'(\theta)/r$$

Hence $f'(\theta) = r[f(\theta + \pi/(2r)) - f(\theta - \pi/(2r))]/2$. For standard Pauli generators, $r = 1$ and the shift is $\pm \pi/2$:

$$\frac{\partial \langle \hat{O} \rangle}{\partial \theta} = \frac{\langle \hat{O}\rangle_{\theta + \pi/2} - \langle \hat{O}\rangle_{\theta - \pi/2}}{2}$$

**For QAOA:** The cost unitary parameters $\gamma_k$ multiply $\hat{H}_C$, which has a richer spectrum than $\pm 1$. The $ZZ$-rotation gates $e^{-i\gamma w_{ij} \hat{Z}_i \hat{Z}_j/2}$ individually satisfy the two-eigenvalue condition with $r = w_{ij}$, so the parameter-shift rule applies gate-by-gate. However, the full circuit expectation value as a function of $\gamma_k$ involves multiple $ZZ$-gates at the same $\gamma_k$, leading to a multi-frequency sinusoidal. In practice, one can either:
- Decompose the cost layer and compute gradients one gate at a time (more circuits).
- Use the general parameter-shift rule for multi-frequency functions (Wierichs et al., 2022).

#### 4.2 Optimizer Comparison

| Optimizer | Type | Gradient-free? | Noisy-tolerant? | Typical use case |
|-----------|------|---------------|-----------------|------------------|
| **COBYLA** | Trust-region, gradient-free | Yes | Moderate | Small $p$, few parameters, hardware runs |
| **L-BFGS-B** | Quasi-Newton | No (needs gradient) | Poor | Noiseless simulation, smooth landscapes |
| **SPSA** | Stochastic gradient approximation | Effectively yes (2 evals/step) | **Excellent** | Hardware, noisy estimates, many parameters |
| **Adam** | Momentum-based gradient descent | No (needs gradient) | Moderate | Differentiable simulators (autograd) |

**COBYLA** (Constrained Optimization BY Linear Approximations): Uses a simplex-like trust region without gradients. Requires $O(d)$ function evaluations per step where $d$ is the parameter dimension. Works well for $d \leq 20$ but scales poorly. On noisy hardware, the trust region can be destabilised by shot noise.

**L-BFGS-B** (Limited-memory BFGS with Box constraints): Approximates the Hessian from gradient history. Converges superlinearly near optima in noiseless settings. Completely unusable on hardware because noise in the gradient estimates corrupts the Hessian approximation, causing divergence.

**SPSA** (Simultaneous Perturbation Stochastic Approximation): Estimates the gradient using only 2 function evaluations (regardless of dimension $d$) by perturbing all parameters simultaneously with random $\pm c_k$ shifts:

$$\hat{g}_k = \frac{f(\boldsymbol{\theta} + c_k \boldsymbol{\Delta}) - f(\boldsymbol{\theta} - c_k \boldsymbol{\Delta})}{2c_k} \boldsymbol{\Delta}^{-1}$$

where $\boldsymbol{\Delta}$ is a random perturbation vector (Rademacher). The gain schedule $a_k, c_k$ must be tuned. SPSA is the go-to optimiser for hardware because it needs constant circuit evaluations per step and is inherently noise-tolerant.

**Adam:** Adaptive moment estimation with bias correction. Works brilliantly with exact gradients from autograd simulators (PennyLane's backpropagation or adjoint differentiation). Not directly usable on hardware unless combined with parameter-shift gradients, which adds overhead.

#### 4.3 Gradient-Free vs Gradient-Based on Noisy Hardware

The fundamental tension:

- **Gradient-based** methods (parameter-shift + Adam/L-BFGS-B) give accurate search directions but require $O(d)$ circuit evaluations per gradient (2 per parameter via parameter-shift), each of which has shot noise. The total gradient estimation error scales as $O(d/\sqrt{S})$ where $S$ is the number of shots.

- **Gradient-free** methods (COBYLA, Nelder-Mead, SPSA) avoid the $O(d)$ overhead but rely on noisy function evaluations to infer descent directions, which may require more iterations to converge.

**Practical guidance:**
- For simulators: use autograd (PennyLane's `diff_method="backprop"`) with Adam. This is exact and fast.
- For hardware with $d \leq 10$ (QAOA $p \leq 5$): COBYLA or SPSA.
- For hardware with $d > 10$: SPSA is the only viable option due to its constant evaluation cost per step.
- Hybrid: optimise on simulator to find a good initialisation, then fine-tune on hardware with SPSA.

---

## Part 2 — Software Architecture

### Project Structure

```
qaoa_explorer/
├── main.py                  # Orchestration, CLI, entry point
├── graph_generator.py       # Graph construction, adjacency → Ising
├── qaoa_circuit.py          # QAOA circuit construction (PennyLane)
├── optimizer.py             # Classical optimization loop
├── landscape_sampler.py     # Grid sampling of (β, γ) space
├── noise_model.py           # Noise simulation (PennyLane + QuTiP)
├── hardware_runner.py       # IBM Quantum / IQM hardware submission
├── error_mitigation.py      # ZNE implementation
├── renderer.py              # OpenGL 3D landscape visualization
├── requirements.txt
├── config.yaml              # Default configuration
└── tests/
    ├── test_graph.py
    ├── test_circuit.py
    ├── test_optimizer.py
    ├── test_sampler.py
    └── test_renderer.py
```

### Data Flow Diagram

```
                     ┌───────────────┐
                     │  config.yaml  │
                     └──────┬────────┘
                            │
                     ┌──────▼────────┐
                     │   main.py     │  CLI / orchestration
                     └──────┬────────┘
                            │
              ┌─────────────┼─────────────────┐
              │             │                  │
     ┌────────▼──────┐ ┌───▼────────────┐ ┌───▼──────────────┐
     │graph_generator│ │landscape_sampler│ │   optimizer.py   │
     │   .py         │ │   .py           │ │                  │
     └────────┬──────┘ └───┬────────────┘ └───┬──────────────┘
              │            │                   │
              │   ┌────────▼─────────┐         │
              └──►│ qaoa_circuit.py  │◄────────┘
                  └────────┬─────────┘
                           │
           ┌───────────────┼───────────────────┐
           │               │                    │
  ┌────────▼──────┐ ┌──────▼───────┐  ┌────────▼──────────┐
  │noise_model.py │ │hardware_     │  │error_mitigation.py│
  │(default.mixed)│ │runner.py     │  │  (ZNE)            │
  └────────┬──────┘ │(IBM/IQM)     │  └────────┬──────────┘
           │        └──────┬───────┘           │
           └───────────────┼───────────────────┘
                           │
                    ┌──────▼───────┐
                    │ renderer.py  │  OpenGL surface
                    │ (β,γ,⟨H_C⟩) │
                    └──────────────┘
```

### Module Specifications

#### 1. `graph_generator.py`

**Purpose:** Construct graphs and convert adjacency representations to Ising Hamiltonians.

```python
from __future__ import annotations
import networkx as nx
import pennylane as qml


def create_graph(
    graph_type: str,
    n_nodes: int,
    seed: int | None = None,
    **kwargs,
) -> nx.Graph:
    """Create a NetworkX graph of the specified type.

    Supported types: 'cycle', 'complete', 'erdos_renyi', 'random_regular',
    'grid_2d', 'petersen', 'custom'. Stochastic families accept a random
    seed for reproducibility. Extra kwargs control family-specific
    parameters (e.g. p=0.5 for Erdos-Renyi, d=3 for random regular).

    Returns:
        A weighted undirected NetworkX graph with unit-weight edges
        by default.
    """
    ...


def adjacency_to_ising(
    graph: nx.Graph,
) -> tuple[list[float], list[tuple[int, int]], float]:
    """Convert a graph's adjacency structure to Ising Hamiltonian coefficients.

    Uses the MaxCut formulation:
        H_C = sum_{(i,j)} w_ij/2 (I - Z_i Z_j)

    For each edge (i,j) with weight w_ij, produces coefficient -w_ij/2
    paired with the Pauli-Z index pair (i,j), plus a constant energy
    offset sum_{(i,j)} w_ij/2.

    Returns:
        Tuple of (coefficients, observable_index_pairs, offset).
    """
    ...


def build_cost_hamiltonian(graph: nx.Graph) -> qml.Hamiltonian:
    """Wrap the Ising conversion into a PennyLane Hamiltonian.

    Each edge contributes a -w_ij/2 * Z_i @ Z_j term. The identity
    offset is included as an explicit qml.Identity term so that
    absolute energy values are preserved.

    Returns:
        PennyLane Hamiltonian ready for circuit evaluation.
    """
    ...
```

*(Full implementation: `report1_scripts/graph_generator.py`)*

#### 2. `qaoa_circuit.py`

**Purpose:** Construct and evaluate QAOA circuits using PennyLane.

```python
from __future__ import annotations
import numpy as np
import pennylane as qml


def cost_layer(
    gamma: float,
    graph_edges: list[tuple[int, int, float]],
) -> None:
    """Apply the cost unitary U_C(gamma) = exp(-i gamma H_C).

    Decomposes each weighted edge (i, j, w_ij) into an IsingZZ rotation:
        exp(-i gamma w_ij Z_i Z_j / 2)
    implemented via qml.IsingZZ(gamma * w_ij, wires=[i, j]).
    All ZZ terms commute, so ordering is irrelevant.
    """
    ...


def mixer_layer(beta: float, n_qubits: int) -> None:
    """Apply the mixer unitary U_B(beta) = exp(-i beta H_B).

    Applies independent single-qubit RX(2*beta) rotations on every
    qubit, since H_B = sum_i X_i factorises across qubits.
    """
    ...


def qaoa_ansatz(
    params: np.ndarray,
    cost_hamiltonian: qml.Hamiltonian,
    n_qubits: int,
    p: int,
) -> None:
    """Apply the full QAOA state-preparation sequence inside a QNode context.

    Prepares |+>^N via a Hadamard layer, then applies p alternating
    rounds of cost and mixer unitaries:
        prod_{k=1}^{p} [U_B(beta_k) U_C(gamma_k)]

    Args:
        params: Array of shape (2p,). First p entries are gamma_k,
                last p entries are beta_k.
    """
    ...


def build_qaoa_circuit(
    cost_hamiltonian: qml.Hamiltonian,
    n_qubits: int,
    p: int,
    device: qml.Device,
) -> qml.QNode:
    """Wrap the QAOA ansatz into a PennyLane QNode bound to the given device.

    The QNode accepts a parameter array of shape (2p,) -- the first p
    entries are the gamma_k values, the last p are the beta_k values --
    and returns <H_C>. Differentiation method is set to 'best'
    (backprop for simulators, parameter-shift for hardware).

    Returns:
        A PennyLane QNode returning the cost expectation value.
    """
    ...
```

*(Full implementation: `report1_scripts/qaoa_circuit.py`)*

#### 3. `optimizer.py`

**Purpose:** Classical optimization of QAOA parameters with multiple backends.

```python
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable
import numpy as np
import pennylane as qml


@dataclass
class OptimizationResult:
    """Container for QAOA optimization results.

    Attributes:
        best_params: Optimal parameter vector of shape (2p,).
        optimal_energy: Best <H_C> value found.
        trajectory: List of (params, energy) pairs at each iteration.
        n_evaluations: Total number of circuit evaluations.
        converged: Whether the optimizer reached the tolerance.
    """
    best_params: np.ndarray
    optimal_energy: float
    trajectory: list[tuple[np.ndarray, float]] = field(default_factory=list)
    n_evaluations: int = 0
    converged: bool = False


def optimize_qaoa(
    qnode: qml.QNode,
    p: int,
    method: str = "COBYLA",
    init_params: np.ndarray | None = None,
    maxiter: int = 500,
    tol: float = 1e-6,
    callback: Callable[[np.ndarray, float, int], None] | None = None,
) -> OptimizationResult:
    """Run classical optimisation of the 2p QAOA parameters.

    Wraps the QNode in a negated objective (QAOA maximises <H_C> but
    standard optimisers minimise).

    Supported methods:
        - 'COBYLA': Gradient-free trust-region via scipy.optimize.minimize.
          Best for few parameters (d <= 20).
        - 'L-BFGS-B': Quasi-Newton requiring exact gradients via the
          parameter-shift rule. Fast near optima in noiseless settings.
        - 'SPSA': Stochastic perturbation gradient approximation. Uses
          only 2 function evaluations per step regardless of dimension,
          with Rademacher perturbation vectors and standard gain schedule
          a_k = a0/(k+A)^alpha, c_k = c0/k^gamma.
        - 'Adam': PennyLane's adaptive moment-based gradient descent,
          suited for autograd-differentiated simulators.

    If init_params is None, gamma_k ~ Uniform(0, pi) and
    beta_k ~ Uniform(0, pi/2). The optional callback is invoked at
    each iteration with (params, energy, iteration).
    """
    ...


def parameter_shift_gradient(
    qnode: qml.QNode,
    params: np.ndarray,
    shift: float = np.pi / 2,
) -> np.ndarray:
    """Compute the gradient of the QNode via the parameter-shift rule.

    For each parameter theta_i, evaluates at theta_i +/- pi/2:
        d<H_C>/d(theta_i) = [f(theta_i + pi/2) - f(theta_i - pi/2)] / 2

    Requires 2d circuit evaluations for d parameters. Hardware-compatible
    (does not rely on autograd). Exact for gates whose generators have
    eigenvalues +/-1 (all Pauli-based rotations in the QAOA circuit).

    Returns:
        Gradient array of shape (2p,).
    """
    ...
```

*(Full implementation: `report1_scripts/optimizer.py`)*

#### 4. `landscape_sampler.py`

**Purpose:** Efficient grid sampling of the $(\beta, \gamma)$ landscape.

```python
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pennylane as qml


@dataclass
class LandscapeData:
    """Container for energy landscape sampling results.

    Attributes:
        gamma_vals: 1D array of gamma sample points, shape (N_gamma,).
        beta_vals: 1D array of beta sample points, shape (N_beta,).
        energies: 2D energy array of shape (N_gamma, N_beta).
        approx_ratios: 2D approximation ratios normalised by C_max,
                        shape (N_gamma, N_beta).
    """
    gamma_vals: np.ndarray
    beta_vals: np.ndarray
    energies: np.ndarray
    approx_ratios: np.ndarray


def sample_landscape(
    qnode: qml.QNode,
    p: int,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
    n_gamma: int,
    n_beta: int,
    fixed_params: np.ndarray | None = None,
    c_max: float = 1.0,
) -> LandscapeData:
    """Evaluate <H_C> on a uniform (N_gamma x N_beta) grid.

    For p=1 the two parameters are swept directly. For p>1, all
    parameters except the first gamma and first beta are held fixed
    at fixed_params, producing a 2D slice through the higher-dimensional
    landscape.

    Returns:
        LandscapeData with energies and approximation ratios.
    """
    ...


def sample_landscape_parallel(
    qnode: qml.QNode,
    p: int,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
    n_gamma: int,
    n_beta: int,
    c_max: float = 1.0,
) -> LandscapeData:
    """Batched grid evaluation using qml.execute over pre-constructed tapes.

    Builds parameter grid via np.meshgrid (indexing='ij'), flattens all
    (N_gamma * N_beta) parameter sets, executes in a single call to
    qml.execute, then reshapes into the 2D energy array. Avoids
    Python-loop overhead; significantly faster on simulators.

    Returns:
        LandscapeData with energies and approximation ratios.
    """
    ...
```

*(Full implementation: `report1_scripts/landscape_sampler.py`)*

#### 5. `noise_model.py`

**Purpose:** Noise simulation using PennyLane's `default.mixed` and QuTiP for Lindblad master equation analysis.

```python
from __future__ import annotations
import numpy as np
import pennylane as qml
import qutip


def create_noisy_device(
    n_qubits: int,
    depolarizing_rate: float = 0.001,
    thermal_relaxation_t1: float = 50e-6,
    thermal_relaxation_t2: float = 70e-6,
    gate_time_1q: float = 50e-9,
    gate_time_2q: float = 300e-9,
) -> qml.Device:
    """Return a PennyLane default.mixed device for mixed-state simulation.

    The device itself is noise-agnostic; noise channels (depolarising,
    thermal relaxation) are inserted explicitly inside the QNode via
    channel operations.

    Returns:
        PennyLane default.mixed device with n_qubits wires.
    """
    ...


def noisy_qaoa_qnode(
    cost_hamiltonian: qml.Hamiltonian,
    n_qubits: int,
    p: int,
    noise_params: dict,
) -> qml.QNode:
    """Construct a QNode on default.mixed that interleaves noise channels.

    After each single-qubit gate: DepolarizingChannel + ThermalRelaxationError.
    After each two-qubit IsingZZ: both qubits receive depolarising noise
    at ~10x the single-qubit rate (reflecting typical hardware) and
    thermal relaxation scaled by the two-qubit gate duration.

    Args:
        noise_params: Dict with keys 'depol_rate', 't1', 't2',
                      'gate_time_1q', 'gate_time_2q'.

    Returns:
        QNode producing <H_C> under gate-dependent incoherent noise.
    """
    ...


def lindblad_simulation(
    initial_state: qutip.Qobj,
    hamiltonian: qutip.Qobj,
    collapse_ops: list[qutip.Qobj],
    tlist: np.ndarray,
    observables: dict[str, qutip.Qobj],
) -> dict[str, np.ndarray]:
    """Run a full Lindblad master equation simulation via qutip.mesolve.

    Evolves the density matrix under:
        drho/dt = -i[H, rho] + sum_k (L_k rho L_k^dag - 1/2 {L_k^dag L_k, rho})

    Collapse operators encode amplitude damping (L = sqrt(gamma_1) sigma_-)
    and pure dephasing (L = sqrt(gamma_phi) Z/2).

    Returns:
        Dict mapping observable labels to time-series arrays of
        shape (len(tlist),).
    """
    ...


def compare_noise_channels(
    n_qubits: int,
    circuit_params: np.ndarray,
    noise_params: dict,
) -> dict[str, float]:
    """Compare PennyLane channel model vs QuTiP Lindblad for the same circuit.

    Runs the same QAOA circuit under both models and compares:
    - Energy expectation values
    - Inter-model state fidelity
    - Output-state purity Tr(rho^2)

    Returns:
        Dict with keys 'energy_pennylane', 'energy_qutip', 'fidelity',
        'purity_pl', 'purity_qt'.
    """
    ...
```

*(Full implementation: `report1_scripts/noise_model.py`)*

#### 6. `hardware_runner.py`

**Purpose:** Submit QAOA circuits to IBM Quantum or IQM via PennyLane device plugins.

```python
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pennylane as qml


@dataclass
class HardwareResult:
    """Container for hardware execution results.

    Attributes:
        expectation_value: Measured <H_C>.
        bitstring_counts: Raw bitstring counts dict.
        shots: Number of shots used.
        backend_name: Name of the quantum backend.
        wall_time: Wall-clock execution time in seconds.
    """
    expectation_value: float
    bitstring_counts: dict[str, int]
    shots: int
    backend_name: str
    wall_time: float


def create_ibm_device(
    n_qubits: int,
    backend_name: str = "ibm_brisbane",
    shots: int = 4096,
    ibm_token: str | None = None,
) -> qml.Device:
    """Return a PennyLane device targeting IBM Quantum hardware.

    Uses the PennyLane-Qiskit plugin (qml.device('qiskit.ibmq', ...)).
    The API token is read from the IBMQ_TOKEN environment variable if
    not supplied explicitly. PennyLane handles transpilation internally.

    Returns:
        PennyLane device bound to an IBM Quantum backend.
    """
    ...


def create_iqm_device(
    n_qubits: int,
    backend_url: str,
    shots: int = 4096,
) -> qml.Device:
    """Return a PennyLane device targeting IQM hardware.

    Uses the PennyLane-IQM plugin.

    Returns:
        PennyLane device bound to an IQM backend.
    """
    ...


def execute_on_hardware(
    qnode: qml.QNode,
    params: np.ndarray,
    cost_hamiltonian: qml.Hamiltonian,
) -> HardwareResult:
    """Run the QAOA circuit on hardware and measure <H_C>.

    Two measurement strategies:
        1. Direct: qml.expval(H_C) -- PennyLane decomposes the
           Hamiltonian into commuting Pauli groups.
        2. Sample-based: qml.sample() returns raw bitstrings of
           shape (S, N). Energy per sample computed classically via
           C(z) = sum_{(i,j)} w_ij/2 * (1 - s_i * s_j) where s_i = 1 - 2*z_i.

    Returns:
        HardwareResult with expectation value, counts, and metadata.
    """
    ...


def select_optimal_qubits(
    backend_name: str,
    n_qubits: int,
    ibm_token: str | None = None,
) -> list[int]:
    """Select the best-performing connected subgraph of physical qubits.

    Loads backend calibration data (T1, T2, CNOT error rates) and builds
    a weighted graph over physical qubits (edges weighted by 1/eps_CX).
    Searches for the size-N connected subgraph maximising
        sum_edges 1/eps_CX + alpha * sum_nodes T1*T2
    via BFS enumeration.

    Returns:
        List of physical qubit indices for the initial_layout.
    """
    ...
```

*(Full implementation: `report1_scripts/hardware_runner.py`)*

#### 7. `error_mitigation.py`

**Purpose:** Zero-noise extrapolation (ZNE) and supporting error mitigation techniques.

```python
from __future__ import annotations
from typing import Callable, Literal
import numpy as np
import pennylane as qml


def zne_extrapolate(
    qnode: qml.QNode,
    params: np.ndarray,
    scale_factors: list[float] = [1.0, 3.0, 5.0],
    extrapolation: Literal["linear", "richardson"] = "richardson",
) -> tuple[float, dict]:
    """Apply zero-noise extrapolation to a hardware or noisy-simulator QNode.

    For each scale factor lambda, amplifies circuit noise via global
    unitary folding and measures <H_C>_lambda. Extrapolates the
    (lambda, E(lambda)) data to lambda=0.

    Returns:
        Tuple of (extrapolated_value, details_dict) where details_dict
        contains 'noisy_values', 'scale_factors', and 'fit_coefficients'.
    """
    ...


def global_fold(
    circuit_fn: Callable,
    scale_factor: float,
) -> Callable:
    """Produce a noise-amplified circuit by global unitary folding.

    For scale factor lambda, replaces circuit U with U(U^dag U)^n where
    n = floor((lambda-1)/2). For non-integer scale factors, partial
    folding applies to the last frac * N_gates gates where
    frac = (lambda-1)/2 - n.

    Returns:
        A callable that applies the folded circuit.
    """
    ...


def richardson_extrapolation(
    scale_factors: np.ndarray,
    noisy_values: np.ndarray,
) -> float:
    """Fit a degree-(M-1) polynomial through M data points via Vandermonde.

    Solves V @ a = E for the coefficient vector a, where V is the
    Vandermonde matrix of the scale factors. Returns a_0 = E(0).

    Returns:
        Extrapolated zero-noise energy value.
    """
    ...


def linear_extrapolation(
    scale_factors: np.ndarray,
    noisy_values: np.ndarray,
) -> float:
    """Fit E(lambda) = a + b*lambda via least-squares.

    Returns:
        a = E(0), the extrapolated zero-noise energy value.
    """
    ...
```

*(Full implementation: `report1_scripts/error_mitigation.py`)*

#### 8. `renderer.py`

**Purpose:** OpenGL 3D surface rendering of the $(\beta, \gamma, \langle H_C \rangle)$ landscape.

```python
from __future__ import annotations
import numpy as np
import glfw
from OpenGL.GL import *


def init_opengl_context(
    width: int = 1280,
    height: int = 720,
    title: str = "QAOA Energy Landscape",
) -> glfw._GLFWwindow:
    """Initialise a GLFW window with OpenGL 3.3 core profile.

    Enables depth testing, 4x MSAA, and alpha blending.

    Returns:
        The GLFW window handle.
    """
    ...


def compile_shader_program(
    vertex_source: str,
    fragment_source: str,
) -> int:
    """Compile and link a vertex-fragment shader pair.

    Returns:
        OpenGL shader program ID.
    """
    ...


def create_surface_vbo(
    gamma_vals: np.ndarray,
    beta_vals: np.ndarray,
    energies: np.ndarray,
) -> tuple[int, int, int, int]:
    """Build a triangle-mesh VAO/VBO/EBO from the (N_gamma x N_beta) grid.

    Each vertex stores 7 float32 values (28 bytes): position
    (gamma, beta, <H_C>), surface normal (n_x, n_y, n_z) from
    finite-difference gradients, and normalised energy in [0,1].
    Two triangles per grid cell. Uses GL_DYNAMIC_DRAW for streaming.

    Returns:
        Tuple of (vao_id, vbo_id, ebo_id, n_indices).
    """
    ...


def update_surface_vbo(
    vbo_id: int,
    gamma_vals: np.ndarray,
    beta_vals: np.ndarray,
    energies: np.ndarray,
) -> None:
    """Overwrite existing VBO vertex data in-place via glBufferSubData."""
    ...


def draw_optimizer_trajectory(
    trajectory: list[tuple[float, float, float]],
    shader_program: int,
) -> None:
    """Render the optimiser's path as a GL_LINE_STRIP above the surface.

    Each vertex carries a normalised progress attribute t in [0,1]
    for blue-to-red colour interpolation. Slight positive z-offset
    (~0.01) prevents z-fighting. Uses GL_STREAM_DRAW.
    """
    ...


def render_side_panel(
    p: int,
    current_energy: float,
    approximation_ratio: float,
    gradient_norm: float,
    iteration: int,
) -> None:
    """Draw a text overlay with live optimisation statistics."""
    ...


class ArcballCamera:
    """Interactive 3D camera with arcball rotation, zoom, and pan.

    Mouse drag rotates around a focus point (yaw/pitch angles),
    scroll wheel zooms by adjusting camera-target distance,
    middle-click pans. Pitch clamped to +/-89 degrees.
    """

    def __init__(
        self,
        target: np.ndarray = np.zeros(3),
        distance: float = 5.0,
        yaw: float = -90.0,
        pitch: float = 30.0,
    ) -> None:
        ...

    def get_view_matrix(self) -> np.ndarray:
        """Return a 4x4 look-at view matrix (float32)."""
        ...

    def get_projection_matrix(
        self,
        aspect_ratio: float,
        fov: float = 45.0,
        near: float = 0.1,
        far: float = 100.0,
    ) -> np.ndarray:
        """Return a 4x4 perspective projection matrix (float32)."""
        ...
```

*(Full implementation: `report1_scripts/renderer.py`)*

#### 9. `main.py`

**Purpose:** CLI entry point and workflow orchestration.

```python
from __future__ import annotations
import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Define and parse CLI arguments.

    Key flags:
        --graph-type    Graph family (cycle, complete, erdos_renyi, etc.)
        --n-nodes       Number of graph nodes
        --p             QAOA depth (number of layers)
        --optimizer     Optimiser backend (COBYLA / L-BFGS-B / SPSA / Adam)
        --grid-size     Landscape sampling resolution
        --noise         Enable noisy simulation (default.mixed)
        --hardware      Run on a real quantum processor
        --backend       Backend name (e.g. ibm_brisbane)
        --zne           Enable zero-noise extrapolation
        --render        Launch OpenGL visualisation
        --config        Path to a YAML configuration file
        --seed          Random seed for reproducibility

    Returns:
        Parsed arguments namespace.
    """
    ...


def run_pipeline(args: argparse.Namespace) -> None:
    """Execute the full QAOA exploration pipeline.

    Steps:
        1. Graph generation -- create_graph + build_cost_hamiltonian;
           compute C_max via brute force (N <= 20).
        2. Circuit construction -- build QNode on appropriate backend
           (noiseless simulator, noisy default.mixed, or hardware).
        3. Landscape sampling -- for p=1 with rendering, sample the
           2D (gamma, beta) energy grid.
        4. Classical optimisation -- run optimize_qaoa, recording the
           trajectory for visualisation.
        5. Error mitigation -- if ZNE requested, apply zne_extrapolate
           at optimal parameters.
        6. Rendering -- if enabled, launch OpenGL window with surface,
           trajectory overlay, and status panel.
    """
    ...
```

*(Full implementation: `report1_scripts/main.py`)*

---

## Part 3 — Implementation Guidelines (Instructional Pseudocode)

### Module 1: `graph_generator.py`

```python
# --- graph_generator.py pseudocode ---

def create_graph(graph_type, n_nodes, seed, **kwargs):
    # Dispatch on graph_type string
    if graph_type == "cycle":
        G = nx.cycle_graph(n_nodes)
    elif graph_type == "complete":
        G = nx.complete_graph(n_nodes)
    elif graph_type == "erdos_renyi":
        G = nx.erdos_renyi_graph(n_nodes, kwargs["p"], seed=seed)
    elif graph_type == "random_regular":
        G = nx.random_regular_graph(kwargs["d"], n_nodes, seed=seed)
    elif graph_type == "grid_2d":
        rows, cols = kwargs.get("rows", n_nodes), kwargs.get("cols", n_nodes)
        G = nx.grid_2d_graph(rows, cols)
        G = nx.convert_node_labels_to_integers(G)  # relabel (i,j) -> int
    elif graph_type == "petersen":
        G = nx.petersen_graph()
    elif graph_type == "custom":
        G = nx.Graph()
        G.add_edges_from(kwargs["edges"])
    # ... other types ...

    # Assign unit weights to edges that don't have a weight
    for u, v in G.edges():
        if "weight" not in G[u][v]:
            G[u][v]["weight"] = 1.0
    return G


def adjacency_to_ising(graph):
    coeffs = []                         # list of float
    obs = []                            # list of (int, int) qubit pairs
    offset = 0.0
    for i, j, w in graph.edges(data="weight"):
        # MaxCut: H_C = sum w_ij/2 (I - Z_i Z_j)
        coeffs.append(-w / 2)
        obs.append((i, j))              # Pauli Z_i Z_j pair
        offset += w / 2                 # identity contribution
    return coeffs, obs, offset


def build_cost_hamiltonian(graph):
    coeffs, obs, offset = adjacency_to_ising(graph)
    # Build PennyLane Hamiltonian
    all_coeffs = coeffs + [offset]
    all_obs = [qml.PauliZ(i) @ qml.PauliZ(j) for i, j in obs] + [qml.Identity(0)]
    return qml.Hamiltonian(all_coeffs, all_obs)
```

*(Full implementation: `report1_scripts/graph_generator.py`)*

### Module 2: `qaoa_circuit.py`

```python
# --- qaoa_circuit.py pseudocode ---

def cost_layer(gamma, graph_edges):
    # graph_edges: list of (i, j, w_ij)
    for i, j, w in graph_edges:
        # IsingZZ implements exp(-i phi Z⊗Z / 2), so pass phi = gamma * w
        qml.IsingZZ(gamma * w, wires=[i, j])


def mixer_layer(beta, n_qubits):
    for q in range(n_qubits):
        # RX(theta) = exp(-i theta X / 2), so theta = 2*beta for exp(-i beta X)
        qml.RX(2 * beta, wires=q)


def extract_edges_from_hamiltonian(cost_hamiltonian):
    """Recover weighted edge list from qml.Hamiltonian."""
    edges = []
    for coeff, obs in zip(cost_hamiltonian.coeffs, cost_hamiltonian.ops):
        if isinstance(obs, qml.Identity):
            continue  # skip the offset term
        # obs is PauliZ(i) @ PauliZ(j); coeff is -w/2
        wires = obs.wires.tolist()       # [i, j]
        w = -2 * coeff                   # recover original weight
        edges.append((wires[0], wires[1], w))
    return edges


def qaoa_ansatz(params, cost_hamiltonian, n_qubits, p):
    # params shape: (2p,) -- [gamma_1, ..., gamma_p, beta_1, ..., beta_p]
    gammas = params[:p]                  # shape (p,)
    betas = params[p:]                   # shape (p,)
    edges = extract_edges_from_hamiltonian(cost_hamiltonian)

    # 1. Prepare |+>^N
    for q in range(n_qubits):
        qml.Hadamard(wires=q)

    # 2. Apply p layers of cost + mixer
    for k in range(p):
        cost_layer(gammas[k], edges)     # U_C(gamma_k)
        mixer_layer(betas[k], n_qubits)  # U_B(beta_k)


def build_qaoa_circuit(cost_hamiltonian, n_qubits, p, device):
    @qml.qnode(device, diff_method="best")
    def circuit(params):
        # params: np.ndarray of shape (2p,)
        qaoa_ansatz(params, cost_hamiltonian, n_qubits, p)
        return qml.expval(cost_hamiltonian)
    return circuit
```

*(Full implementation: `report1_scripts/qaoa_circuit.py`)*

### Module 3: `optimizer.py`

```python
# --- optimizer.py pseudocode ---

def optimize_qaoa(qnode, p, method, init_params, maxiter, tol, callback):
    d = 2 * p                         # total number of parameters
    eval_count = 0                    # mutable circuit evaluation counter
    trajectory = []                   # list of (params_copy, energy)

    # --- Initialisation ---
    if init_params is None:
        gammas = np.random.uniform(0, np.pi, size=p)       # shape (p,)
        betas = np.random.uniform(0, np.pi / 2, size=p)    # shape (p,)
        params = np.concatenate([gammas, betas])            # shape (d,)
    else:
        params = init_params.copy()

    # Negated objective: optimisers minimise, QAOA maximises <H_C>
    def objective(theta):
        nonlocal eval_count
        eval_count += 1
        val = float(qnode(theta))
        return -val

    # --- COBYLA ---
    if method == "COBYLA":
        from scipy.optimize import minimize
        result = minimize(
            objective, params, method="COBYLA",
            options={"maxiter": maxiter, "rhobeg": 0.5},
            callback=lambda x: trajectory.append((x.copy(), -objective(x))),
        )
        best_params = result.x

    # --- L-BFGS-B ---
    elif method == "L-BFGS-B":
        from scipy.optimize import minimize
        def jac(theta):
            return -parameter_shift_gradient(qnode, theta)
        result = minimize(
            objective, params, method="L-BFGS-B", jac=jac,
            options={"maxiter": maxiter, "ftol": tol},
        )
        best_params = result.x

    # --- SPSA ---
    elif method == "SPSA":
        # Gain schedule constants
        a0, c0 = 0.1, 0.1
        alpha, gamma_spsa = 0.602, 0.101
        A = 0.1 * maxiter

        for k in range(1, maxiter + 1):
            a_k = a0 / (k + A) ** alpha       # step size schedule
            c_k = c0 / k ** gamma_spsa         # perturbation size schedule

            # Rademacher random perturbation vector: each entry +1 or -1
            delta = np.random.choice([-1, 1], size=d)   # shape (d,)

            # Two-point gradient estimate (only 2 evaluations per step)
            f_plus = objective(params + c_k * delta)
            f_minus = objective(params - c_k * delta)
            g_hat = (f_plus - f_minus) / (2 * c_k) * (1.0 / delta)  # shape (d,)

            params = params - a_k * g_hat      # gradient descent step
            trajectory.append((params.copy(), -float(qnode(params))))

            if callback:
                callback(params, trajectory[-1][1], k)

        best_params = params

    # --- Adam ---
    elif method == "Adam":
        opt = qml.AdamOptimizer(stepsize=0.05)
        for k in range(maxiter):
            params, cost = opt.step_and_cost(lambda t: -qnode(t), params)
            energy = -cost
            trajectory.append((params.copy(), energy))
            if callback:
                callback(params, energy, k)
            # Convergence check
            if k > 0 and abs(trajectory[-1][1] - trajectory[-2][1]) < tol:
                break
        best_params = params

    return OptimizationResult(
        best_params=best_params,
        optimal_energy=-objective(best_params),
        trajectory=trajectory,
        n_evaluations=eval_count,
        converged=...,
    )


def parameter_shift_gradient(qnode, params, shift=np.pi / 2):
    d = len(params)                    # number of parameters
    grad = np.zeros(d)                 # shape (d,)
    for i in range(d):
        shifted_plus = params.copy()
        shifted_minus = params.copy()
        shifted_plus[i] += shift       # theta_i + pi/2
        shifted_minus[i] -= shift      # theta_i - pi/2
        # Exact parameter-shift formula (valid for eigenvalues +/-1)
        grad[i] = (float(qnode(shifted_plus)) - float(qnode(shifted_minus))) / 2.0
    # Total: 2*d circuit evaluations
    return grad
```

*(Full implementation: `report1_scripts/optimizer.py`)*

### Module 4: `landscape_sampler.py`

```python
# --- landscape_sampler.py pseudocode ---

def sample_landscape(qnode, p, gamma_range, beta_range, n_gamma, n_beta,
                     fixed_params, c_max):
    gamma_vals = np.linspace(gamma_range[0], gamma_range[1], n_gamma)  # (N_gamma,)
    beta_vals = np.linspace(beta_range[0], beta_range[1], n_beta)      # (N_beta,)
    energies = np.zeros((n_gamma, n_beta))                              # (N_gamma, N_beta)

    for i, gamma in enumerate(gamma_vals):
        for j, beta in enumerate(beta_vals):
            if p == 1:
                params = np.array([gamma, beta])                        # (2,)
            else:
                # 2D slice through higher-dimensional landscape
                params = fixed_params.copy()                            # (2p,)
                params[0] = gamma      # overwrite first gamma
                params[p] = beta       # overwrite first beta
            energies[i, j] = float(qnode(params))

    approx_ratios = energies / c_max                                    # (N_gamma, N_beta)
    return LandscapeData(gamma_vals, beta_vals, energies, approx_ratios)


def sample_landscape_parallel(qnode, p, gamma_range, beta_range,
                              n_gamma, n_beta, c_max):
    gamma_vals = np.linspace(gamma_range[0], gamma_range[1], n_gamma)  # (N_gamma,)
    beta_vals = np.linspace(beta_range[0], beta_range[1], n_beta)      # (N_beta,)

    # Build full parameter grid via meshgrid with 'ij' indexing
    gamma_grid, beta_grid = np.meshgrid(gamma_vals, beta_vals, indexing='ij')
    # gamma_grid shape: (N_gamma, N_beta)
    # beta_grid  shape: (N_gamma, N_beta)

    # Flatten into 1D arrays of length N_gamma * N_beta
    gamma_flat = gamma_grid.ravel()                                     # (N_gamma*N_beta,)
    beta_flat = beta_grid.ravel()                                       # (N_gamma*N_beta,)

    # Assemble parameter matrix for all grid points
    n_points = n_gamma * n_beta
    all_params = np.zeros((n_points, 2 * p))                            # (N_total, 2p)
    all_params[:, 0] = gamma_flat      # first gamma for each point
    all_params[:, p] = beta_flat       # first beta for each point

    # Convert each parameter set into a PennyLane QuantumTape
    tapes = []
    for idx in range(n_points):
        tape = qml.tape.QuantumScript(...)  # build tape from qnode with all_params[idx]
        tapes.append(tape)

    # Execute all tapes in one batch -- exploits simulator vectorisation
    results = qml.execute(tapes, device=qnode.device, gradient_fn=None)
    # results: list of length N_total, each a scalar expectation value

    energies = np.array(results).reshape(n_gamma, n_beta)               # (N_gamma, N_beta)
    approx_ratios = energies / c_max
    return LandscapeData(gamma_vals, beta_vals, energies, approx_ratios)
```

*(Full implementation: `report1_scripts/landscape_sampler.py`)*

### Module 5: `noise_model.py`

```python
# --- noise_model.py pseudocode ---

def create_noisy_device(n_qubits, **noise_kwargs):
    # Device is noise-agnostic; channels are inserted inside the QNode
    return qml.device("default.mixed", wires=n_qubits)


def noisy_qaoa_qnode(cost_hamiltonian, n_qubits, p, noise_params):
    dev = create_noisy_device(n_qubits, **noise_params)
    edges = extract_edges_from_hamiltonian(cost_hamiltonian)

    @qml.qnode(dev, diff_method="best")
    def circuit(params):
        gammas = params[:p]                                   # shape (p,)
        betas = params[p:]                                    # shape (p,)

        # --- Hadamard layer + noise ---
        for q in range(n_qubits):
            qml.Hadamard(wires=q)
            # Insert single-qubit noise after each Hadamard
            qml.DepolarizingChannel(noise_params["depol_rate"], wires=q)
            qml.ThermalRelaxationError(
                noise_params["t1"], noise_params["t2"],
                noise_params["gate_time_1q"], wires=q
            )

        for k in range(p):
            # --- Cost layer + two-qubit noise ---
            for i, j, w in edges:
                qml.IsingZZ(gammas[k] * w, wires=[i, j])
                # Two-qubit noise: ~10x single-qubit depol rate
                for qubit in [i, j]:
                    qml.DepolarizingChannel(
                        noise_params["depol_rate"] * 10, wires=qubit
                    )
                    qml.ThermalRelaxationError(
                        noise_params["t1"], noise_params["t2"],
                        noise_params["gate_time_2q"], wires=qubit
                    )

            # --- Mixer layer + single-qubit noise ---
            for q in range(n_qubits):
                qml.RX(2 * betas[k], wires=q)
                qml.DepolarizingChannel(noise_params["depol_rate"], wires=q)
                qml.ThermalRelaxationError(
                    noise_params["t1"], noise_params["t2"],
                    noise_params["gate_time_1q"], wires=q
                )

        return qml.expval(cost_hamiltonian)
    return circuit


def lindblad_simulation(initial_state, hamiltonian, collapse_ops, tlist,
                        observables):
    # initial_state: qutip.Qobj density matrix of shape (2^N, 2^N)
    # hamiltonian:   qutip.Qobj Hamiltonian of shape (2^N, 2^N)
    # collapse_ops:  list of qutip.Qobj, each shape (2^N, 2^N)
    # tlist:         np.ndarray of time points, shape (T,)

    # Build per-qubit collapse operators as tensor products:
    #   Amplitude damping: L_ad = sqrt(1/T1) * sigma_minus, tensored
    #   Pure dephasing:    L_pd = sqrt(1/T2 - 1/(2*T1)) * Z/2, tensored
    for q in range(n_qubits):
        gamma_1 = 1.0 / T1
        gamma_phi = 1.0 / T2 - 1.0 / (2 * T1)
        # Tensor product: I (x) ... (x) sigma_minus (x) ... (x) I
        L_ad = np.sqrt(gamma_1) * tensor_single_op(qutip.destroy(2), q, n_qubits)
        L_pd = np.sqrt(gamma_phi) * tensor_single_op(qutip.sigmaz() / 2, q, n_qubits)
        collapse_ops.extend([L_ad, L_pd])

    # Solve Lindblad master equation
    result = qutip.mesolve(
        hamiltonian, initial_state, tlist,
        c_ops=collapse_ops,
        e_ops=[observables[label] for label in observables],
    )
    # result.expect[k] is np.ndarray of shape (T,) for k-th observable
    return {label: result.expect[k] for k, label in enumerate(observables)}


def compare_noise_channels(n_qubits, circuit_params, noise_params):
    # 1. PennyLane channel model
    qnode_pl = noisy_qaoa_qnode(...)
    energy_pl = float(qnode_pl(circuit_params))
    # Extract density matrix from PennyLane (via qml.state())
    rho_pl = ...                         # shape (2^N, 2^N) complex

    # 2. QuTiP Lindblad model
    # Initialise |+>^N as a density matrix
    plus = (qutip.basis(2, 0) + qutip.basis(2, 1)).unit()
    rho0 = qutip.tensor([plus * plus.dag()] * n_qubits)    # (2^N, 2^N)
    # Step through each gate as piecewise-constant Hamiltonian evolution
    # ... (build H_gate for each gate, evolve for gate_time) ...
    result = lindblad_simulation(rho0, H_total, collapse_ops, tlist, ...)
    energy_qt = result["H_C"][-1]        # final time-step energy

    # 3. Compare
    fidelity = qutip.fidelity(rho_pl_qobj, rho_qt)
    purity_pl = np.real(np.trace(rho_pl @ rho_pl))
    purity_qt = np.real((rho_qt * rho_qt).tr())
    return {"energy_pennylane": energy_pl, "energy_qutip": energy_qt,
            "fidelity": fidelity, "purity_pl": purity_pl, "purity_qt": purity_qt}
```

*(Full implementation: `report1_scripts/noise_model.py`)*

### Module 6: `hardware_runner.py`

```python
# --- hardware_runner.py pseudocode ---

def create_ibm_device(n_qubits, backend_name, shots, ibm_token):
    import os
    token = ibm_token or os.environ.get("IBMQ_TOKEN")
    # PennyLane-Qiskit plugin handles transpilation internally
    dev = qml.device("qiskit.ibmq", wires=n_qubits,
                      backend=backend_name, shots=shots, ibmqx_token=token)
    return dev


def create_iqm_device(n_qubits, backend_url, shots):
    dev = qml.device("iqm.iqm_device", wires=n_qubits,
                      url=backend_url, shots=shots)
    return dev


def execute_on_hardware(qnode, params, cost_hamiltonian):
    import time
    t0 = time.perf_counter()

    # Strategy 1: direct expectation (preferred)
    energy = float(qnode(params))         # qml.expval(H_C) internally

    # Strategy 2 (alternative): sample-based post-processing
    # samples = qnode_sample(params)      # shape (S, N) binary array
    # for s in range(S):
    #     z = samples[s]                  # shape (N,) binary
    #     spin = 1 - 2 * z               # {0,1} -> {+1,-1}
    #     for i, j, w in edges:
    #         C_s += w/2 * (1 - spin[i] * spin[j])
    # energy = np.mean(C_values)

    wall_time = time.perf_counter() - t0
    return HardwareResult(energy, ..., shots, backend_name, wall_time)


def select_optimal_qubits(backend_name, n_qubits, ibm_token):
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService(channel="ibm_quantum", token=ibm_token)
    backend = service.backend(backend_name)
    props = backend.properties()

    # Build weighted graph over physical qubits
    coupling = backend.configuration().coupling_map   # list of [q_i, q_j] pairs
    G = nx.Graph()
    for qi, qj in coupling:
        cx_error = props.gate_error("cx", [qi, qj])
        G.add_edge(qi, qj, weight=1.0 / cx_error)    # weight = 1 / eps_CX

    for q in G.nodes():
        G.nodes[q]["t1"] = props.t1(q)
        G.nodes[q]["t2"] = props.t2(q)

    # Search for the best size-N connected subgraph via BFS enumeration
    best_score = -np.inf
    best_subgraph = None
    alpha = 0.1
    for seed_node in G.nodes():
        subgraph = bfs_subgraph(G, seed_node, n_qubits)   # first N reachable nodes
        score = (
            sum(G[u][v]["weight"] for u, v in subgraph.edges())
            + alpha * sum(G.nodes[q]["t1"] * G.nodes[q]["t2"]
                          for q in subgraph.nodes())
        )
        if score > best_score:
            best_score = score
            best_subgraph = subgraph
    return list(best_subgraph.nodes())
```

*(Full implementation: `report1_scripts/hardware_runner.py`)*

### Module 7: `error_mitigation.py`

```python
# --- error_mitigation.py pseudocode ---

def global_fold(circuit_fn, scale_factor):
    """Produce a noise-amplified circuit by global unitary folding."""
    n = int((scale_factor - 1) / 2)           # full folding repetitions
    frac = (scale_factor - 1) / 2 - n         # fractional remainder

    # Record original circuit operations from PennyLane tape
    original_tape = qml.tape.QuantumScript(...)
    ops = original_tape.operations              # list of gate operations
    N_gates = len(ops)

    def folded_circuit(params):
        # 1. Apply original circuit U
        for op in ops:
            qml.apply(op)

        # 2. Apply n full folds: (U^dag U) repeated n times
        for _ in range(n):
            # U^dag: apply ops in reverse, each adjointed
            for op in reversed(ops):
                qml.adjoint(op)
            # U: apply ops forward again
            for op in ops:
                qml.apply(op)

        # 3. Partial fold for non-integer scale factors
        if frac > 0:
            n_partial = int(round(frac * N_gates))  # gates to partially fold
            tail_ops = ops[-n_partial:]              # last n_partial gates
            # Adjoint of the tail
            for op in reversed(tail_ops):
                qml.adjoint(op)
            # Re-apply the tail
            for op in tail_ops:
                qml.apply(op)

    return folded_circuit


def zne_extrapolate(qnode, params, scale_factors, extrapolation):
    noisy_values = []                              # shape will be (M,)
    for lam in scale_factors:
        # Build a folded QNode at scale factor lambda
        folded_qnode = build_folded_qnode(qnode, lam)   # uses global_fold
        E_lam = float(folded_qnode(params))
        noisy_values.append(E_lam)

    scale_arr = np.array(scale_factors)            # shape (M,)
    noisy_arr = np.array(noisy_values)             # shape (M,)

    if extrapolation == "linear":
        E0 = linear_extrapolation(scale_arr, noisy_arr)
    elif extrapolation == "richardson":
        E0 = richardson_extrapolation(scale_arr, noisy_arr)

    return E0, {"noisy_values": noisy_arr, "scale_factors": scale_arr}


def richardson_extrapolation(scale_factors, noisy_values):
    # scale_factors: shape (M,)   noisy_values: shape (M,)
    M = len(scale_factors)
    # Build Vandermonde matrix V[i, j] = scale_factors[i]^j
    V = np.vander(scale_factors, N=M, increasing=True)  # shape (M, M)
    # Solve V @ a = E for coefficient vector a
    a = np.linalg.solve(V, noisy_values)                # shape (M,)
    return a[0]                           # a_0 = E(lambda=0)


def linear_extrapolation(scale_factors, noisy_values):
    # Fit E(lambda) = a + b*lambda via least-squares
    # Design matrix: [[1, lambda_1], [1, lambda_2], ...]
    A = np.column_stack([np.ones_like(scale_factors), scale_factors])  # (M, 2)
    coeffs, *_ = np.linalg.lstsq(A, noisy_values, rcond=None)        # (2,)
    return coeffs[0]                      # a = E(0)
```

*(Full implementation: `report1_scripts/error_mitigation.py`)*

### Module 8: `renderer.py`

```python
# --- renderer.py pseudocode ---

def init_opengl_context(width, height, title):
    glfw.init()
    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
    glfw.window_hint(glfw.SAMPLES, 4)              # 4x MSAA
    window = glfw.create_window(width, height, title, None, None)
    glfw.make_context_current(window)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_MULTISAMPLE)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    return window


def create_surface_vbo(gamma_vals, beta_vals, energies):
    n_gamma, n_beta = len(gamma_vals), len(beta_vals)

    # --- Vertex data: 7 floats per vertex (28 bytes) ---
    # [gamma, beta, energy, nx, ny, nz, norm_energy]
    vertices = np.zeros((n_gamma * n_beta, 7), dtype=np.float32)
    for i in range(n_gamma):
        for j in range(n_beta):
            idx = i * n_beta + j
            E = energies[i, j]
            vertices[idx, 0:3] = [gamma_vals[i], beta_vals[j], E]
            # Finite-difference surface normals: n = (-dE/dgamma, -dE/dbeta, 1)
            dEdg = (energies[min(i+1, n_gamma-1), j] - energies[max(i-1, 0), j]) / ...
            dEdb = (energies[i, min(j+1, n_beta-1)] - energies[i, max(j-1, 0)]) / ...
            normal = np.array([-dEdg, -dEdb, 1.0])
            normal /= np.linalg.norm(normal)
            vertices[idx, 3:6] = normal
            vertices[idx, 6] = ...                  # normalised energy in [0, 1]

    # --- Index data: 2 triangles per grid cell ---
    indices = []                                     # uint32
    for i in range(n_gamma - 1):
        for j in range(n_beta - 1):
            tl = i * n_beta + j                      # top-left
            tr = (i + 1) * n_beta + j                # top-right
            bl = i * n_beta + (j + 1)                # bottom-left
            br = (i + 1) * n_beta + (j + 1)          # bottom-right
            indices.extend([tl, bl, tr, bl, br, tr])
    # Total triangles: 2 * (N_gamma - 1) * (N_beta - 1)

    # --- GPU upload ---
    vao = glGenVertexArrays(1)
    vbo = glGenBuffers(1)
    ebo = glGenBuffers(1)
    glBindVertexArray(vao)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_DYNAMIC_DRAW)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, ...)
    # Attribute 0: position (3 floats, offset 0, stride 28)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 28, ctypes.c_void_p(0))
    # Attribute 1: normal (3 floats, offset 12, stride 28)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 28, ctypes.c_void_p(12))
    # Attribute 2: energy (1 float, offset 24, stride 28)
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 28, ctypes.c_void_p(24))
    ...
    return vao, vbo, ebo, len(indices)


def draw_optimizer_trajectory(trajectory, shader_program):
    n = len(trajectory)
    # Each vertex: (gamma, beta, energy, t) -- t in [0,1] for colour
    verts = np.zeros((n, 4), dtype=np.float32)
    for k, (gamma, beta, energy) in enumerate(trajectory):
        verts[k] = [gamma, beta, energy + 0.01, k / max(n - 1, 1)]  # z-offset
    # Upload as GL_STREAM_DRAW, draw as GL_LINE_STRIP
    ...
```

*(Full implementation: `report1_scripts/renderer.py`)*

### Module 9: `main.py`

```python
# --- main.py pseudocode ---

def parse_args():
    parser = argparse.ArgumentParser(description="QAOA Energy Landscape Explorer")
    parser.add_argument("--graph-type", choices=[...], default="cycle")
    parser.add_argument("--n-nodes", type=int, default=6)
    parser.add_argument("--p", type=int, default=1)
    parser.add_argument("--optimizer", choices=["COBYLA","L-BFGS-B","SPSA","Adam"])
    parser.add_argument("--grid-size", type=int, default=50)
    parser.add_argument("--noise", action="store_true")
    parser.add_argument("--hardware", action="store_true")
    parser.add_argument("--backend", type=str, default="ibm_brisbane")
    parser.add_argument("--zne", action="store_true")
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--config", type=str, default="config.yaml")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def run_pipeline(args):
    # 1. Graph generation
    graph = create_graph(args.graph_type, args.n_nodes, seed=args.seed)
    H_C = build_cost_hamiltonian(graph)
    # Brute-force C_max for N <= 20
    if args.n_nodes <= 20:
        c_max = brute_force_max_cut(graph)
    ...

    # 2. Device selection
    if args.hardware:
        dev = create_ibm_device(args.n_nodes, args.backend)
    elif args.noise:
        dev = create_noisy_device(args.n_nodes, ...)
    else:
        dev = qml.device("default.qubit", wires=args.n_nodes)

    qnode = build_qaoa_circuit(H_C, args.n_nodes, args.p, dev)

    # 3. Landscape sampling (p=1 with rendering)
    if args.render and args.p == 1:
        landscape = sample_landscape_parallel(
            qnode, args.p, (0, np.pi), (0, np.pi/2),
            args.grid_size, args.grid_size, c_max)

    # 4. Classical optimisation
    result = optimize_qaoa(qnode, args.p, method=args.optimizer, ...)

    # 5. ZNE error mitigation
    if args.zne and (args.noise or args.hardware):
        E0, details = zne_extrapolate(qnode, result.best_params,
                                       scale_factors=[1, 3, 5])

    # 6. Rendering
    if args.render:
        window = init_opengl_context(1280, 720, "QAOA Landscape")
        vao, vbo, ebo, n_idx = create_surface_vbo(
            landscape.gamma_vals, landscape.beta_vals, landscape.energies)
        # Main render loop
        while not glfw.window_should_close(window):
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            # Draw surface, trajectory, side panel
            ...
            glfw.swap_buffers(window)
            glfw.poll_events()
        glfw.terminate()
```

*(Full implementation: `report1_scripts/main.py`)*

---

## Part 4 — OpenGL Visualization Design

### 1. Full Rendering Pipeline for the 3D $(\beta, \gamma, \langle H_C \rangle)$ Surface

#### Vertex Buffer Layout

Each vertex has 7 float32 values (28 bytes per vertex):

| Offset (bytes) | Attribute | Size | Description |
|:-:|:-:|:-:|:--|
| 0 | Position | 3 floats | $(x, y, z) = (\gamma, \beta, \langle H_C \rangle)$ |
| 12 | Normal | 3 floats | Surface normal $(n_x, n_y, n_z)$ for Phong lighting |
| 24 | Energy | 1 float | Normalised energy $\in [0, 1]$ for colour mapping |

The surface mesh is a regular grid of $(N_\gamma \times N_\beta)$ vertices connected by $(N_\gamma - 1)(N_\beta - 1) \times 2$ triangles (two per grid cell). The index buffer stores triangle indices as uint32.

#### GLSL Shaders

**Vertex Shader (`surface.vert`):**

The vertex shader receives three per-vertex attributes at layout locations 0–2: position $(\gamma, \beta, \langle H_C\rangle)$, the surface normal, and the normalised energy scalar. It transforms the position into world space via the model matrix, computes the correct world-space normal using the inverse-transpose of the model matrix (to handle non-uniform scaling), passes the normalised energy through as a varying, and outputs clip-space position via the MVP chain.

*(Full source: `report1_scripts/shaders/surface.vert`)*

**Fragment Shader (`surface.frag`):**

The fragment shader implements Phong illumination with an energy-based colour map. The normalised energy $e \in [0,1]$ is mapped to a hue via $h = (1 - e) \times 0.66$, producing a blue-to-red gradient (blue = low energy, red = high energy) through an HSV→RGB conversion. Phong shading combines three terms:

- **Ambient** ($0.2 \times$ base colour) — constant fill light.
- **Diffuse** ($\max(\hat{n} \cdot \hat{l},\, 0) \times$ base colour) — Lambertian term.
- **Specular** ($0.5 \times [\max(\hat{v} \cdot \hat{r},\, 0)]^{32}$) — shininess highlights.

The output alpha is $0.9$ (slight transparency) to allow seeing the trajectory through the surface edges.

*(Full source: `report1_scripts/shaders/surface.frag`)*

#### Camera Controls (Arcball Rotation)

The `ArcballCamera` class (Module 8) handles:
- **Left-click drag:** Rotate the view around the scene centre. The yaw/pitch angles update proportionally to mouse displacement.
- **Scroll wheel:** Zoom in/out by adjusting the camera-to-target distance.
- **Middle-click drag:** Pan the target point in the view plane.

GLFW callbacks wire mouse events to the camera. A cursor-position callback computes the delta from the last position and calls `camera.on_mouse_drag(dx, dy)` when the left button is held. A scroll callback forwards the vertical offset to `camera.on_scroll(yoffset)`. Both callbacks are registered via `glfw.set_cursor_pos_callback` and `glfw.set_scroll_callback`.

*(Full implementation: `report1_scripts/renderer.py`)*

#### Live Surface Update During Optimization

As the optimizer progresses, the landscape surface can be updated in two ways:

1. **Full recompute (for changing landscape slice in $p > 1$):** When the fixed parameters change, resample the 2D landscape and call `update_surface_vbo` with new energy data. This uses `glBufferSubData` for an in-place update. The recompute is throttled (e.g., every 10 optimiser steps) to avoid bottlenecking the render loop.

2. **Streaming VBO update (double buffering):** Two VBOs are allocated. While one is bound for rendering, the other is updated on a background thread. After the update completes, the buffer IDs are swapped. This prevents stalls in the render loop during data upload.

### 2. Optimizer Trajectory Overlay

The optimizer trajectory is rendered as a 3D polyline lying on (or slightly above) the energy surface. Each vertex is $(\gamma_k, \beta_k, \langle H_C \rangle_k + \epsilon)$ where $\epsilon \approx 0.01$ prevents z-fighting.

**Trajectory shaders:**

The trajectory vertex shader (`trajectory.vert`) takes position and a progress scalar $t \in [0,1]$ as inputs, applies the combined model–view–projection matrix, and forwards $t$ to the fragment stage. The fragment shader (`trajectory.frag`) interpolates colour from blue $(0.0, 0.3, 1.0)$ at $t = 0$ to red $(1.0, 0.1, 0.0)$ at $t = 1$ using GLSL's `mix` function, with full opacity.

*(Full source: `report1_scripts/shaders/trajectory.vert` and `trajectory.frag`)*

The trajectory VBO is rebuilt every frame (or every N iterations) from the `trajectory_log` list using `GL_STREAM_DRAW`.

### 3. Side Panel

Using Dear ImGui (via `imgui[glfw]` Python bindings) for the info panel. Each frame, a fixed-size window is positioned in the top-left corner and populated with live text fields: circuit depth $p$, current $\langle \hat{H}_C \rangle$, approximation ratio, gradient norm, iteration count, optimiser name, and graph metadata. The ImGui draw data is rendered to the OpenGL framebuffer after the 3D scene.

*(Full implementation: `report1_scripts/renderer.py`)*

---

## Part 5 — Hardware Execution Plan

### 1. IBM Quantum Submission via PennyLane

#### 1.1 Transpilation for a Specific Backend

PennyLane's `qiskit.ibmq` device handles transpilation internally. When the QNode is executed, PennyLane:

1. Converts PennyLane operations to Qiskit gates.
2. Calls `qiskit.transpile()` targeting the backend's native gate set (typically $\{CX, \text{ID}, R_Z, \sqrt{X}, X\}$ for IBM Eagle/Heron processors).
3. Maps virtual qubits to physical qubits based on the coupling map.
4. Applies routing passes (SABRE, stochastic swap) to handle non-adjacent CNOT requirements.

To control transpilation, the PennyLane-Qiskit device accepts a `transpile_options` dictionary specifying `optimization_level` (0–3, where 3 enables the most aggressive gate cancellation and resynthesis passes) and `initial_layout` (the list of physical qubit indices returned by `select_optimal_qubits`). These are forwarded directly to `qiskit.transpile()`.

*(Full implementation: `report1_scripts/hardware_runner.py`)*

#### 1.2 Qubit Selection via Coupling Map Analysis

Use `select_optimal_qubits()` (Module 6). The key heuristic:

1. Load calibration data for all qubits (T1, T2) and all CNOT/ECR gates (error rates).
2. Build a weighted graph over physical qubits.
3. Find the connected subgraph of size $N$ that maximises $\sum_{\text{edges}} 1/\epsilon_{CX} + \alpha \sum_{\text{nodes}} T_1 T_2$.
4. Use this subgraph as the `initial_layout` in transpilation.

For the QAOA circuit on $C_4$ (4 qubits, ring topology), select 4 physical qubits that form a ring in the coupling map with minimal CNOT error.

#### 1.3 Shots: `qml.sample` vs `qml.expval`

**`qml.expval(H_C)` on hardware:**
PennyLane decomposes $H_C = \sum_k c_k P_k$ into Pauli terms, groups commuting terms, and measures each group separately. For MaxCut on $|E|$ edges, $H_C$ has $|E|$ ZZ terms — all mutually commuting (all diagonal), so they can be measured in a single circuit run. This is optimal.

**`qml.sample()` on hardware:**
Returns raw bitstrings, shape `(shots, n_qubits)`. You compute $\langle H_C \rangle$ classically:

$$\langle H_C \rangle \approx \frac{1}{S}\sum_{s=1}^{S} C(\mathbf{z}^{(s)})$$

where $C(\mathbf{z})$ is the MaxCut cost evaluated on bitstring $\mathbf{z}$. This avoids PennyLane's Pauli grouping overhead and gives access to the full distribution for post-processing.

**Shot budget:** For standard deviation $\sigma$ in the energy estimate:

$$\sigma = \frac{\Delta E}{\sqrt{S}}$$

where $\Delta E = E_{\max} - E_{\min}$ is the spectral range of $H_C$. To achieve precision $\epsilon$, need $S \geq (\Delta E / \epsilon)^2$ shots. For MaxCut on 4-10 qubits, $\Delta E \sim 10$, so for $\epsilon = 0.1$: $S \geq 10000$ shots.

#### 1.4 Readout Error Handling

Readout (measurement) errors are the dominant noise source on IBM hardware (typically 1-5% per qubit). Strategies:

1. **Readout error mitigation:** Prepare all $2^N$ computational basis states, measure, and build a calibration matrix $M$ where $M_{ij} = P(\text{measure } i | \text{prepared } j)$. Then invert: $\mathbf{p}_{\text{corrected}} = M^{-1} \mathbf{p}_{\text{raw}}$.

2. **Twirled readout error mitigation (TREX):** Randomly flip measurement outcomes using Pauli-$X$ twirling to convert general readout error into a symmetric channel, then correct analytically.

3. **PennyLane integration:** Use `qml.transforms.mitigate_with_zne` or custom post-processing on the raw counts.

### 2. Error Mitigation Strategy

#### 2.1 ZNE with Linear + Richardson Extrapolation

**Global unitary folding:** To amplify noise by factor $\lambda$, replace circuit $U$ with $U(U^\dagger U)^{(\lambda-1)/2}$:

- $\lambda = 1$: original circuit (depth $d$)
- $\lambda = 3$: circuit followed by $U^\dagger U$ (depth $3d$)
- $\lambda = 5$: circuit followed by $U^\dagger U U^\dagger U$ (depth $5d$)

For each $\lambda$, measure $\langle H_C \rangle_\lambda$.

**Linear extrapolation:** Fit $E(\lambda) = a + b\lambda$ and extrapolate to $\lambda = 0$: $E_0 = a = E(\lambda=1) - b$. Assumes noise effect is linear in noise strength.

**Richardson extrapolation:** Fit $E(\lambda) = \sum_{k=0}^{M-1} a_k \lambda^k$ through $M$ data points and evaluate at $\lambda = 0$. More accurate but more sensitive to statistical noise.

**ZNE for QAOA MaxCut on $C_4$ at $p = 1$:**

With the circuit depth being $|E| + N = 4 + 4 = 8$ gates (4 ZZ-rotations + 4 RX gates):
- $\lambda = 1$: 8 gates
- $\lambda = 3$: 24 gates
- $\lambda = 5$: 40 gates

Use 3 scale factors with Richardson extrapolation by calling `zne_extrapolate(hardware_qnode, optimal_params, scale_factors=[1, 3, 5], extrapolation="richardson")`.

*(Full implementation: `report1_scripts/error_mitigation.py`)*

#### 2.2 Pauli Twirling as Precursor

Pauli twirling converts general coherent errors on 2-qubit gates into stochastic Pauli channels (incoherent noise), which ZNE handles more effectively. For each CNOT (or ECR) gate in the transpiled circuit:

1. Sample random Pauli pair $(P_a, P_b)$ from the twirl group.
2. Replace $\text{CNOT}$ with $P_b \cdot \text{CNOT} \cdot P_a$.
3. Average over many random twirls.

This should be applied before ZNE. PennyLane does not have built-in Pauli twirling, so it must be implemented via tape transforms.

#### 2.3 Realistic Fidelity Expectations

| Configuration | Expected $\langle H_C \rangle$ error | Approx. ratio degradation |
|:--|:-:|:-:|
| $p=1$, 4 qubits, noiseless | 0 | 0 |
| $p=1$, 4 qubits, ibm_brisbane (raw) | ~0.2-0.5 | ~5-15% |
| $p=1$, 4 qubits, ibm_brisbane + ZNE | ~0.05-0.15 | ~2-5% |
| $p=1$, 10 qubits, ibm_brisbane (raw) | ~0.5-1.5 | ~10-30% |
| $p=1$, 10 qubits, ibm_brisbane + ZNE | ~0.2-0.5 | ~5-15% |
| $p=2$, 10 qubits, ibm_brisbane (raw) | ~1.0-3.0 | ~20-50% |
| $p=2$, 10 qubits, ibm_brisbane + ZNE | ~0.3-1.0 | ~10-25% |

These are order-of-magnitude estimates based on typical 2024-2025 IBM Eagle processor performance with $\sim 1\%$ two-qubit gate error rates and $\sim 2\%$ readout error rates.

### 3. Benchmarking Checklist

**Metrics to record for each experiment:**

| Metric | Noiseless Sim | Noisy Sim | Hardware | Hardware+ZNE |
|:--|:-:|:-:|:-:|:-:|
| $\langle H_C \rangle$ | ✓ | ✓ | ✓ | ✓ |
| Approximation ratio $r = \langle H_C \rangle / C_{\max}$ | ✓ | ✓ | ✓ | ✓ |
| Optimal $(\gamma^*, \beta^*)$ | ✓ | ✓ | ✓ | ✓ |
| Optimizer iterations to convergence | ✓ | ✓ | ✓ | ✓ |
| Total circuit evaluations | ✓ | ✓ | ✓ | ✓ |
| Wall-clock time | ✓ | ✓ | ✓ | ✓ |
| Bitstring distribution (top-10 samples) | ✓ | ✓ | ✓ | ✓ |
| Probability of sampling optimal cut | ✓ | ✓ | ✓ | ✓ |

**Plots:**
1. $r$ vs $p$ (circuit depth): noiseless, noisy, hardware, hardware+ZNE on same axes.
2. Energy landscape surface comparison: noiseless vs noisy (as 3D surfaces or heatmaps side by side).
3. Optimizer convergence curves: energy vs iteration for each configuration.
4. Histogram of sampled bitstring costs: noiseless vs hardware.

---

## Part 6 — Research Contribution Angle

### 1. Three Directions Toward a Publishable Paper

**Direction 1: QAOA Landscape Phase Transitions in Frustrated Graphs**

Study the QAOA energy landscape topology as a function of graph frustration index. For anti-ferromagnetic Ising models on non-bipartite graphs (triangular lattices, Petersen graph, random $d$-regular), map the landscape at $p = 1, 2, 3$ and quantify:
- Number of local minima (via landscape sampling + basin detection).
- Hessian spectrum at critical points (positive definite, negative definite, saddle index).
- Correlation between landscape complexity and classical hardness (brute-force time, simulated annealing residual energy).

This extends the 2025 Boulebnane-Montanaro result by moving beyond regular graphs with uniform fields to arbitrary frustrated topologies.

**Direction 2: Noise-Induced Landscape Smoothing and Concentration**

Systematically compare the $(\beta, \gamma)$ landscape under noiseless, depolarising, and thermal relaxation noise at varying strengths. The hypothesis: moderate noise smooths away spurious local minima, potentially improving optimizer convergence (stochastic regularisation). Measure:
- Gradient norm distribution across the landscape as a function of noise rate.
- Number of local minima as noise increases (using persistent homology / topological data analysis on the landscape).
- Optimal parameter stability: how much do $(\gamma^*, \beta^*)$ shift under noise?

This connects to the barren plateau literature and provides practical guidance on when noise helps vs hurts variational optimization.

**Direction 3: Adaptive QAOA Depth Selection via Landscape Diagnostics**

Propose a protocol that examines the $p$-level landscape structure to decide whether increasing $p$ is worthwhile:
- Compute the landscape Hessian at the $p$-level optimum.
- If the smallest Hessian eigenvalue is "sufficiently negative" (surface is steeply curved), the current $p$ is under-parameterised and increasing $p$ will help.
- If the Hessian is nearly flat (plateau), further depth is wasteful.
- Validate on families of graphs with known $p^*$ (the smallest $p$ for which QAOA is provably optimal).

### 2. Most Relevant 2024–2025 Papers

1. **Boulebnane & Montanaro (2025).** "Solving optimization problems with local light cones on fixedangle QAOA." *Quantum*. — Phase transition in QAOA performance with longitudinal field.

2. **Kremenetski, Hogg, & Hadfield (2025).** "Quantum Alternating Operator Ansatz (QAOA+): landscape structure with constraints." *PRX Quantum*. — Extended QAOA landscape analysis beyond MaxCut.

3. **Shaydulin et al. (2024).** "Evidence of scaling advantage for the quantum approximate optimization algorithm on a classically intractable problem." *Science Advances*. — First empirical evidence of QAOA advantage at scale.

4. **Wurtz & Love (2024).** "Counterdiabatically optimized QAOA." *Physical Review A*. — Improved initialisation strategies via counterdiabatic terms.

5. **Cerezo et al. (2024).** "Does provable absence of barren plateaus imply classical simulability?" *Physical Review X*. — Fundamental limits on trainability and simulability.

6. **Pelofske et al. (2024).** "Scaling whole-chip QAOA for higher-order ising spin glass models on heavy-hex graphs." *npj Quantum Information*. — Large-scale QAOA benchmarking on IBM hardware.

7. **Zhou, Wang, Choi et al. (2020).** "Quantum Approximate Optimization Algorithm: Performance, Mechanism, and Implementation on Near-Term Devices." *Physical Review X*. — INTERP and FOURIER strategies (seminal reference).

### 3. Interesting Graph Families and Problem Types

| Graph Family | Why Interesting | Landscape Feature |
|:--|:--|:--|
| **Random 3-regular** | Phase transition at clause density for SAT-like problems | Complex landscape with many near-degenerate minima |
| **Petersen graph** | Highly symmetric, non-planar, frustrated | Discrete symmetry in landscape, symmetry-breaking effects |
| **Triangular lattice** | Geometrically frustrated for AF Ising | Exponentially degenerate ground states, flat landscape regions |
| **Erdős-Rényi $G(n, p_c)$** | Percolation transition at $p_c = 1/n$ | Landscape topology changes at the percolation threshold |
| **Complete bipartite $K_{n,n}$** | MaxCut is trivial but landscape is instructive | Single clean global minimum — contrast with random graphs |
| **Weighted Sherrington-Kirkpatrick** | Complete graph with Gaussian random weights | Full-RSB spin glass: hardest possible QAOA landscape |

The most visually dramatic landscapes come from frustrated systems (triangular, 3-regular non-bipartite, SK model) where near-degenerate minima create complex topography.
