"""Noise simulation: PennyLane channel model and QuTiP Lindblad master equation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pennylane as qml
import qutip

from qaoa_circuit import _extract_edges_from_hamiltonian


@dataclass
class NoisyDeviceConfig:
    device: qml.Device
    depolarizing_rate: float
    thermal_relaxation_t1: float
    thermal_relaxation_t2: float
    gate_time_1q: float
    gate_time_2q: float

def create_noisy_device(
    n_qubits: int,
    depolarizing_rate: float = 0.01,
    thermal_relaxation_t1: float = 50e-6,
    thermal_relaxation_t2: float = 70e-6,
    gate_time_1q: float = 50e-9,
    gate_time_2q: float = 300e-9,
) -> NoisyDeviceConfig:
    """Create a PennyLane mixed-state device for noisy simulation.

    Returns a ``default.mixed`` device.  Actual noise is injected inside
    the QNode via explicit channel gates (``DepolarizingChannel``,
    ``ThermalRelaxationError``), not via device-level noise models.

    Args:
        n_qubits: Number of wires.
        depolarizing_rate: Single-qubit depolarising probability per gate.
        thermal_relaxation_t1: Amplitude-damping T₁ time (seconds).
        thermal_relaxation_t2: Dephasing T₂ time (seconds).
        gate_time_1q: Duration of a single-qubit gate (seconds).
        gate_time_2q: Duration of a two-qubit gate (seconds).

    Returns:
        A ``default.mixed`` PennyLane device.
    """
    return NoisyDeviceConfig(
        device=qml.device("default.mixed", wires=n_qubits),
        depolarizing_rate=depolarizing_rate,
        thermal_relaxation_t1=thermal_relaxation_t1,
        thermal_relaxation_t2=thermal_relaxation_t2,
        gate_time_1q=gate_time_1q,
        gate_time_2q=gate_time_2q,
    )


def noisy_qaoa_qnode(
    cost_hamiltonian: qml.Hamiltonian,
    n_qubits: int,
    p: int,
    noise_params: dict,
) -> qml.QNode:
    """Build a noisy QAOA QNode on a mixed-state simulator.

    The circuit mirrors the ideal QAOA ansatz but inserts noise channels
    after every gate:

    - **After each single-qubit gate** (Hadamard, RX):
        ``qml.DepolarizingChannel(depol_rate)`` and
        ``qml.ThermalRelaxationError(pe=0, t1, t2, gate_time_1q)``.
    - **After each two-qubit gate** (IsingZZ):
        ``qml.DepolarizingChannel(depol_rate × 10)`` on *both* qubits and
        ``qml.ThermalRelaxationError(pe=0, t1, t2, gate_time_2q)``
        on *both* qubits.

    Expected keys in *noise_params*:
        ``depol_rate``, ``t1``, ``t2``, ``gate_time_1q``, ``gate_time_2q``.

    Args:
        cost_hamiltonian: QAOA cost Hamiltonian.
        n_qubits: Number of qubits.
        p: QAOA depth.
        noise_params: Dictionary of noise parameters.

    Returns:
        A ``qml.QNode`` on ``default.mixed`` that returns ``⟨H_C⟩``.
    """
    depol_rate  = noise_params["depol_rate"]
    t1          = noise_params["t1"]
    t2          = noise_params["t2"]
    gt_1q       = noise_params["gate_time_1q"]
    gt_2q       = noise_params["gate_time_2q"]

    dev = qml.device('default.mixed', wires=n_qubits)
    edges = _extract_edges_from_hamiltonian(cost_hamiltonian)

    def _depol_and_thermal(wires, gate_time, rate):
        for w in wires:
            qml.DepolarizingChannel(rate, wires=w)
            qml.ThermalRelaxationError(pe=0, t1=t1, t2=t2, tg=gate_time, wires=w)

    @qml.qnode(dev)
    def circuit(gammas, betas):
        for q in range(n_qubits):
            qml.Hadamard(wires=q)
            _depol_and_thermal([q], gt_1q, depol_rate)

        for k in range(p):
            for i, j, w in edges:
                qml.IsingZZ(2*gammas[k]*w, wires=[i, j])
                _depol_and_thermal([i, j], gt_2q, 10*depol_rate)

            for q in range(n_qubits):
                qml.RX(2*betas[k], wires=q)
                _depol_and_thermal([q], gt_1q, depol_rate)
        return qml.expval(cost_hamiltonian)
    
    return circuit



def lindblad_simulation(
    initial_state: qutip.Qobj,
    hamiltonian: qutip.Qobj,
    collapse_ops: list[qutip.Qobj],
    tlist: np.ndarray,
    observables: list[qutip.Qobj],
) -> dict[str, np.ndarray]:
    """Run a Lindblad master-equation simulation with QuTiP.

    Solves:

        dρ/dt = −i[H, ρ] + Σ_k (L_k ρ L_k† − ½{L_k†L_k, ρ})

    where ``L_k`` are the *collapse_ops* encoding noise channels:
      - Amplitude damping: ``L = √γ₁ · σ⁻``
      - Pure dephasing:    ``L = √γ_φ · σ_z / 2``

    Uses ``qutip.mesolve`` with ``nsteps=10000``.

    Args:
        initial_state: Initial ket or density matrix (``qutip.Qobj``).
        hamiltonian: System Hamiltonian (``qutip.Qobj`` operator).
        collapse_ops: List of Lindblad collapse operators.
        tlist: 1-D array of time points.
        observables: Operators whose expectation values are recorded.

    Returns:
        Dictionary ``{"obs_0": array, "obs_1": array, ...}`` with
        expectation-value time series for each observable.
    """
    for obs in observables:
        if obs.dims != initial_state.dims:
            raise ValueError("Observable dimensions must match initial state.")
    result = qutip.mesolve(
        hamiltonian,
        initial_state,
        tlist,
        c_ops=collapse_ops,
        e_ops=observables,
        nsteps=10000,
    )
    return {f"obs_{i}": np.array(result.expect[i]) for i in range(len(observables))}


def compare_noise_channels(
    n_qubits: int,
    circuit_params: np.ndarray,
    noise_params: dict,
    cost_hamiltonian: qml.Hamiltonian,
) -> dict[str, float]:
    """Compare PennyLane channel model vs QuTiP Lindblad for the same QAOA circuit.

    Steps:
      1. **PennyLane**: run ``noisy_qaoa_qnode`` to get the noisy energy.
      2. **QuTiP**: build the system as a tensor product of *n_qubits* qubits,
         prepare |+⟩^N, construct collapse operators (amplitude damping
         ``√(1/T₁)·σ⁻`` and dephasing ``√γ_φ·σ_z/2`` per qubit, with
         ``γ_φ = 1/T₂ − 1/(2T₁)``), then step through each QAOA gate layer
         using ``qutip.mesolve`` for the appropriate gate duration.
      3. Compute the cost-Hamiltonian expectation on the QuTiP final state.
      4. Optionally compute the state fidelity between the two density matrices.

    Args:
        n_qubits: Number of qubits.
        circuit_params: QAOA parameter vector (shape ``(2·p,)``).
        n_qubits: Number of qubits.
        circuit_params: QAOA parameter vector (shape ``(2·p,)``).
        noise_params: Noise parameter dictionary
            (keys: ``depol_rate``, ``t1``, ``t2``,
            ``gate_time_1q``, ``gate_time_2q``).
        cost_hamiltonian: QAOA cost Hamiltonian.

    Returns:
        Dictionary with keys ``"pennylane_energy"``, ``"qutip_energy"``,
        ``"fidelity_pl_qt"``, ``"purity_pl"``, ``"purity_qt"``.
    """
    # PennyLane simulation
    pl_circuit = noisy_qaoa_qnode(
        cost_hamiltonian=cost_hamiltonian,  # Build from qaoa_circuit.build_cost_hamiltonian
        n_qubits=n_qubits,
        p=len(circuit_params) // 2,
        noise_params=noise_params,
    )
    pennylane_energy = pl_circuit(circuit_params[:len(circuit_params)//2], circuit_params[len(circuit_params)//2:])

    # --- shared derived values ---
    t1         = noise_params["t1"]
    t2         = noise_params["t2"]
    gt_1q      = noise_params["gate_time_1q"]
    gt_2q      = noise_params["gate_time_2q"]
    depol_rate = noise_params["depol_rate"]
    p          = len(circuit_params) // 2
    gammas     = circuit_params[:p]
    betas      = circuit_params[p:]
    edges      = _extract_edges_from_hamiltonian(cost_hamiltonian)

    # --- QuTiP N-qubit operator helpers ---
    si = qutip.qeye(2)
    sx, sz, sm = qutip.sigmax(), qutip.sigmaz(), qutip.sigmam()

    def _local(op: qutip.Qobj, q: int) -> qutip.Qobj:
        lst = [si] * n_qubits
        lst[q] = op
        return qutip.tensor(lst)

    def _two(op1: qutip.Qobj, op2: qutip.Qobj, qi: int, qj: int) -> qutip.Qobj:
        lst = [si] * n_qubits
        lst[qi] = op1
        lst[qj] = op2
        return qutip.tensor(lst)

    # --- initial state |+><+|^N ---
    plus = (qutip.basis(2, 0) + qutip.basis(2, 1)).unit()
    rho  = qutip.ket2dm(qutip.tensor([plus] * n_qubits))

    # --- Lindblad collapse operators per qubit ---
    # gamma_phi = 1/T2 - 1/(2*T1)  (pure dephasing; clamped to 0 for unphysical params)
    gamma_1   = 1.0 / t1
    gamma_phi = max(1.0 / t2 - 1.0 / (2.0 * t1), 0.0)
    c_ops: list[qutip.Qobj] = []
    for q in range(n_qubits):
        c_ops.append(np.sqrt(gamma_1)   * _local(sm, q))        # amplitude damping
        c_ops.append(np.sqrt(gamma_phi) * _local(sz, q) / 2.0)  # pure dephasing

    def _evolve_gate(rho_in: qutip.Qobj, gate_ham: qutip.Qobj, dt: float) -> qutip.Qobj:
        """Evolve density matrix under gate_ham with continuous Lindblad noise for dt."""
        res = qutip.mesolve(
            gate_ham, rho_in, np.array([0.0, dt]), c_ops=c_ops, e_ops=[]
        )
        return res.states[-1]

    # Hadamard layer: e^{-i*pi/2*(X+Z)/sqrt(2)}, generator = pi/(2*gt_1q)*(X+Z)/sqrt(2)
    for q in range(n_qubits):
        H_had = (np.pi / (2.0 * gt_1q)) * (_local(sx, q) + _local(sz, q)) / np.sqrt(2.0)
        rho = _evolve_gate(rho, H_had, gt_1q)

    # QAOA layers
    for k in range(p):
        # cost layer: IsingZZ(2*gamma*w) = e^{-i*gamma*w*ZZ}, generator = (gamma*w/gt_2q)*ZZ
        for i, j, w in edges:
            H_zz = (float(gammas[k]) * w / gt_2q) * _two(sz, sz, i, j)
            rho = _evolve_gate(rho, H_zz, gt_2q)
        # mixer layer: RX(2*beta) = e^{-i*beta*X}, generator = (beta/gt_1q)*X
        for q in range(n_qubits):
            H_rx = (float(betas[k]) / gt_1q) * _local(sx, q)
            rho = _evolve_gate(rho, H_rx, gt_1q)

    # Build QuTiP cost Hamiltonian via qml.matrix (handles arbitrary Pauli structures)
    dims_N = [[2] * n_qubits, [2] * n_qubits]
    H_cost_qt = qutip.Qobj(
        np.zeros((2**n_qubits, 2**n_qubits), dtype=complex), dims=dims_N
    )
    for coeff, op in zip(cost_hamiltonian.coeffs, cost_hamiltonian.ops):
        mat = qml.matrix(op, wire_order=list(range(n_qubits)))
        H_cost_qt += float(np.real(coeff)) * qutip.Qobj(mat, dims=dims_N)

    qutip_energy = float(qutip.expect(H_cost_qt, rho))

    # PennyLane density matrix (same noisy circuit, returns dm for fidelity/purity)
    dev_dm = qml.device("default.mixed", wires=n_qubits)

    def _noise(wires: list[int], gate_time: float, rate: float) -> None:
        for w in wires:
            qml.DepolarizingChannel(rate, wires=w)
            qml.ThermalRelaxationError(pe=0, t1=t1, t2=t2, tg=gate_time, wires=w)

    @qml.qnode(dev_dm)
    def _dm_circuit(gammas_p: np.ndarray, betas_p: np.ndarray) -> np.ndarray:
        for q in range(n_qubits):
            qml.Hadamard(wires=q)
            _noise([q], gt_1q, depol_rate)
        for k in range(p):
            for i, j, w in edges:
                qml.IsingZZ(2.0 * gammas_p[k] * w, wires=[i, j])
                _noise([i, j], gt_2q, 10.0 * depol_rate)
            for q in range(n_qubits):
                qml.RX(2.0 * betas_p[k], wires=q)
                _noise([q], gt_1q, depol_rate)
        return qml.density_matrix(wires=list(range(n_qubits)))

    pl_dm_arr = _dm_circuit(gammas, betas)
    pl_rho = qutip.Qobj(pl_dm_arr.reshape(2**n_qubits, 2**n_qubits), dims=dims_N)

    fidelity_pl_qt = float(qutip.fidelity(pl_rho, rho))
    purity_pl      = float((pl_rho * pl_rho).tr().real)
    purity_qt      = float((rho    * rho   ).tr().real)

    return {
        "pennylane_energy": pennylane_energy,
        "qutip_energy": qutip_energy,
        "fidelity_pl_qt": fidelity_pl_qt,
        "purity_pl": purity_pl,
        "purity_qt": purity_qt,
    }
