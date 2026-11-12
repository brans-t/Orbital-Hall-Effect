"""Assemble and validate full OAM matrices without silently symmetrizing."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ohk.eigen import inverse_gaps


@dataclass
class AssemblyReport:
    hermiticity_error: float
    trace: complex
    diagonal_filled: bool
    max_abs: float


def diagonal_closed_form(basis):
    inv = inverse_gaps(basis.eps, basis.regularization_eta)
    return 2 * np.sum(np.imag(basis.v[0] * basis.v[1].conj()) * inv, axis=1)


def full_matrix(variant, basis, hermiticity_tol=1e-8):
    result = np.array(variant.matrix(basis), dtype=complex, copy=True)
    if result.shape != (basis.n_bands,) * 2 or not np.isfinite(result).all():
        raise ValueError("OAM variant must return a finite square matrix with num_bands rows")
    if not variant.diagonal:
        np.fill_diagonal(result, diagonal_closed_form(basis))
    error = float(np.linalg.norm(result - result.conj().T) / max(np.linalg.norm(result), 1.0))
    if error > hermiticity_tol:
        raise ValueError(f"Non-Hermitian OAM matrix: residual {error:.3g}")
    return result, AssemblyReport(
        error, complex(np.trace(result)), not variant.diagonal, float(np.max(abs(result)))
    )


def compare(variants, basis):
    return {variant.name: full_matrix(variant, basis)[0] for variant in variants}
