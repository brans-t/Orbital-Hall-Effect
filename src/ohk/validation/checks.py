"""Algebraic residuals, with explicit band-phase alignment for symmetry checks."""

from __future__ import annotations
import numpy as np
from ohk.eigen import EigenBasis
from ohk.oam.assemble import full_matrix


def check_hermiticity(L, reference_scale=None):
    L = np.asarray(L)
    denominator = max(np.linalg.norm(L) if reference_scale is None else reference_scale, 1.0)
    return float(np.linalg.norm(L - L.conj().T) / denominator)


def check_trace(L):
    """Diagnostic only: OAM has no general zero-trace sum rule, even for two bands."""
    return complex(np.trace(L))


def _aligned_residual(model, variant, k, symmetry, sign):
    a, b = EigenBasis.at(model, k), EigenBasis.at(model, -np.asarray(k))
    target = (
        a.U.conj()
        if symmetry == "T"
        else np.asarray(model.spec.metadata.get("inversion_operator", np.eye(a.n_bands))) @ a.U
    )
    # Express the -k operator in the symmetry-transformed +k band basis.
    overlap = b.U.conj().T @ target
    left = overlap.conj().T @ full_matrix(variant, b)[0] @ overlap
    right = full_matrix(variant, a)[0]
    return float(np.max(abs(left - sign * (right.conj() if symmetry == "T" else right))))


def check_time_reversal(model, variant, k, tol=1e-10):
    """Spinless complex-conjugation T in the orbital basis; phases aligned before comparison."""
    return _aligned_residual(model, variant, k, "T", -1)


def check_inversion(model, variant, k, tol=1e-10):
    """Uses spec.metadata['inversion_operator']; defaults to identity."""
    return _aligned_residual(model, variant, k, "I", 1)


def check_gauge_invariance(model, variant, k, rng=None, tol=1e-10):
    """Historical name: checks L' = D^dag L D, not entrywise matrix invariance."""
    rng = np.random.default_rng(0) if rng is None else rng
    basis = EigenBasis.at(model, k)
    phase = np.exp(1j * rng.uniform(-np.pi, np.pi, basis.n_bands))
    rotated = EigenBasis(
        basis.k,
        basis.eps,
        basis.U * phase,
        phase.conj()[None, :, None] * basis.v * phase[None, None, :],
        basis.dH,
    )
    expected = phase.conj()[:, None] * full_matrix(variant, basis)[0] * phase[None, :]
    return float(np.max(abs(full_matrix(variant, rotated)[0] - expected)))


def check_band_diagonal_real(L, tol=1e-10):
    return float(np.max(abs(np.diag(L).imag)))
