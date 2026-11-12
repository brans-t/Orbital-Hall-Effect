"""Kubo OHC, reported numerically in units of e/(2*pi) in 2D.

The raw normalized-BZ average is (2/A) sum_n f_n Im[J_nm v_mn]/gap^2 in units of e.
Multiply by 2*pi exactly once to express it in e/(2*pi). In 3D the unit additionally
contains inverse model length. Broadening is opt-in and applies to both OAM and Kubo gaps.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import warnings
import numpy as np
from ohk.eigen import EigenBasis
from ohk.oam.assemble import full_matrix
from ohk.transport.current import orbital_current
from ohk.transport.grids import BZGrid


@dataclass
class KuboResult:
    sigma: float
    sigma_terms: np.ndarray
    degenerate_fraction: float
    nen: int  # Maximum occupied-band count; use occupied_counts for a metal.
    occupied_counts: np.ndarray | None = None
    min_gap: float = float("inf")
    regularization_eta: float = 0.0
    units: str = "e/(2*pi)"

    def __float__(self):
        return self.sigma


def _response(
    model,
    variant,
    grid,
    current_direction=0,
    field_direction=1,
    degeneracy_tol=1e-8,
    subband_regularization=False,
    eta=None,
):
    if grid.is_path:
        raise ValueError("A high-symmetry path is not a BZ integration grid")
    if not np.allclose(grid.lattice_vectors, model.lattice(), atol=1e-12, rtol=1e-12):
        raise ValueError("grid and model lattice vectors differ")
    d = model.spec.periodic_dim
    if not 0 <= current_direction < d or not 0 <= field_direction < d:
        raise ValueError("current and field directions must be valid Cartesian axes")
    if not np.isfinite(degeneracy_tol) or degeneracy_tol <= 0:
        raise ValueError("degeneracy_tol must be positive and finite")
    broadening = (degeneracy_tol if eta is None else float(eta)) if subband_regularization else 0.0
    if eta is not None and not subband_regularization:
        raise ValueError("eta requires subband_regularization=True")
    if (
        broadening < 0
        or not np.isfinite(broadening)
        or (subband_regularization and broadening == 0)
    ):
        raise ValueError("regularization eta must be positive and finite")
    energies, kernels, gaps = [], [], []
    for k in grid.kpoints:
        basis = EigenBasis.at(model, k)
        gaps.append(basis.min_gap())
        basis.regularization_eta = broadening
        L = full_matrix(variant, basis)[0]
        J = orbital_current(L, basis.v, current_direction)
        gap = basis.eps[:, None] - basis.eps[None, :]
        off = ~np.eye(basis.n_bands, dtype=bool)
        denominator = np.zeros_like(gap)
        denominator[off] = 1 / (gap[off] ** 2 + broadening**2)
        kernels.append(2 * np.imag(J * basis.v[field_direction].T * denominator).sum(axis=1))
        energies.append(basis.eps)
    gaps = np.asarray(gaps)
    fraction = float(np.mean(gaps < degeneracy_tol))
    if fraction > 0.001:
        warnings.warn(
            f"{fraction:.2%} of k points are near-degenerate; check convergence",
            UserWarning,
            stacklevel=3,
        )
    weights = grid.weights if grid.normalized else grid.weights / grid.reciprocal_volume()
    scale = 2 * np.pi / grid.cell_area()
    return (
        np.asarray(energies),
        np.asarray(kernels) * scale,
        weights,
        fraction,
        gaps.min(),
        broadening,
    )


def orbital_hall_conductivity(
    model,
    variant,
    grid,
    mu=0.0,
    current_direction=0,
    field_direction=1,
    degeneracy_tol=1e-8,
    subband_regularization=False,
    eta=None,
):
    if not np.isfinite(mu):
        raise ValueError("mu must be finite")
    eps, kernels, weights, fraction, gap, broadening = _response(
        model,
        variant,
        grid,
        current_direction,
        field_direction,
        degeneracy_tol,
        subband_regularization,
        eta,
    )
    occupation = eps < mu
    terms = (occupation * kernels).sum(axis=1)
    counts = occupation.sum(axis=1)
    units = "e/(2*pi)" if model.spec.periodic_dim == 2 else "e/(2*pi*length)"
    return KuboResult(
        float(weights @ terms),
        terms,
        fraction,
        int(counts.max()),
        counts,
        float(gap),
        broadening,
        units,
    )


def fermi_scan(model, variant, grid, energies: Sequence[float], **kwargs):
    """Cache the band-resolved response once and vary only the occupation."""
    mus = np.asarray(energies, float)
    if mus.ndim != 1 or not np.isfinite(mus).all():
        raise ValueError("energies must be a finite one-dimensional sequence")
    eps, kernels, weights, *_ = _response(model, variant, grid, **kwargs)
    return np.array([weights @ ((eps < mu) * kernels).sum(axis=1) for mu in mus])


def orbital_hall_symbol(model, variant, grid, mu=0.0, degeneracy_tol=1e-8):
    return orbital_hall_conductivity(model, variant, grid, mu, 1, 0, degeneracy_tol)


def convergence_study(model, variant, grid_sizes, mu=0.0, **kwargs):
    """Return mesh values and successive differences; no unjustified extrapolation."""
    sizes = list(grid_sizes)
    if not sizes:
        raise ValueError("grid_sizes must not be empty")
    values = np.array(
        [
            orbital_hall_conductivity(
                model, variant, BZGrid.uniform(model.lattice(), n), mu=mu, **kwargs
            ).sigma
            for n in sizes
        ]
    )
    return {
        "grid_sizes": sizes,
        "sigma": values,
        "successive_changes": np.diff(values),
        "extrapolated": None,
    }
