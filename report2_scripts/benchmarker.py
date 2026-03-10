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


def run_benchmark(pes_data: object) -> BenchmarkResult:
    """Compute benchmark metrics from a PES scan.

    Starts with methods = ["HF", "VQE"] and energies dict containing
    HF, VQE, and FCI arrays from pes_data. Appends "CISD" and "CCSD(T)"
    if their arrays are not None in pes_data.

    Computes for each method:
      - Correlation energy recovery:
        r = (E_method - E_HF) / (E_FCI - E_HF)
        Use np.where to avoid division by zero when |denom| < 1e-10.
        A value near 1.0 indicates the method captures nearly all
        correlation energy.
      - Absolute error vs FCI: |E_method - E_FCI|.

    Args:
        pes_data: PESData from a scan_pes run, containing energies_hf,
                  energies_vqe, energies_fci, and optionally energies_cisd
                  and energies_ccsd_t arrays.

    Returns:
        BenchmarkResult with per-method energy arrays, correlation recovery
        ratios, and errors vs FCI.
    """
    methods = ["HF", "VQE"]
    energies = {
        "HF": pes_data.energies_hf,
        "VQE": pes_data.energies_vqe,
        "FCI": pes_data.energies_fci
    }
    if pes_data.energies_cisd is not None:
        methods.append("CISD")
        energies["CISD"] = pes_data.energies_cisd
    if pes_data.energies_ccsd_t is not None:
        methods.append("CCSD(T)")
        energies["CCSD(T)"] = pes_data.energies_ccsd_t

    correlation_recovery = {}
    errors_vs_fci = {}
    for method in methods:
        if method == "FCI":
            continue
        e_method = energies[method]
        e_hf = energies["HF"]
        e_fci = energies["FCI"]
        denom = e_fci - e_hf
        r = np.where(np.abs(denom) > 1e-10, (e_method - e_hf) / denom, 0.0)
        correlation_recovery[method] = r
        errors_vs_fci[method] = np.abs(e_method - e_fci)

    return BenchmarkResult(
        bond_lengths=pes_data.bond_lengths,
        methods=methods,
        energies=energies,
        correlation_recovery=correlation_recovery,
        errors_vs_fci=errors_vs_fci
    )


def print_benchmark_table(result: BenchmarkResult,
                          bond_lengths_subset: np.ndarray | None = None) -> str:
    """Format a benchmark comparison table as a string.

    If bond_lengths_subset is None, picks 5 evenly spaced indices via
    np.linspace(0, len-1, 5, dtype=int). Otherwise, finds the closest
    bond length index for each requested value via np.argmin.

    Table columns: R (A) | method1 | method2 | ... | FCI | |VQE-FCI|
    Energies formatted to 6 decimal places, R to 3 decimal places.

    Args:
        result: BenchmarkResult to display.
        bond_lengths_subset: Optional array of specific bond lengths to show.

    Returns:
        Multi-line string with header, separator, and data rows.
    """
    if bond_lengths_subset is None:
        indices = np.linspace(0, len(result.bond_lengths)-1, 5, dtype=int)
    else:
        indices = [np.argmin(np.abs(result.bond_lengths - R)) for R in bond_lengths_subset]

    header = "R (A) | " + " | ".join(f"{m:>10}" for m in result.methods) + " | FCI (ref) | |VQE-FCI|"
    separator = "-" * len(header)
    rows = []
    for idx in indices:
        R = result.bond_lengths[idx]
        row = f"{R:6.3f} | " + " | ".join(f"{result.energies[m][idx]:10.6f}" for m in result.methods) + f" | {result.energies['FCI'][idx]:10.6f} | {result.errors_vs_fci['VQE'][idx]:10.6f}"
        rows.append(row)
    table = "\n".join([header, separator] + rows)
    return table
