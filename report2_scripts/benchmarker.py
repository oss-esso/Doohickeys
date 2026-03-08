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
    ...


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
    ...
