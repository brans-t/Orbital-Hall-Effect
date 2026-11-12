"""Finite real-space reference and deterministic random multiband models."""

from __future__ import annotations
import numpy as np


def _compare(lattice, variant, momenta=None, band_pairs=None, tol=1e-8):
    from ohk.eigen import EigenBasis
    from ohk.oam.assemble import full_matrix

    result = lattice.compute()
    allowed = lattice.allowed_momenta()
    if momenta is None:
        selected = allowed
    else:
        frac = np.asarray(momenta, float)
        if frac.ndim != 2 or frac.shape[1] != lattice.model.spec.periodic_dim:
            raise ValueError("momenta must be fractional reciprocal vectors")
        index = frac * np.asarray(lattice.size)
        if not np.allclose(index, np.rint(index), atol=1e-10, rtol=0):
            raise ValueError("requested momenta are not allowed by the finite supercell")
        selected = frac @ lattice.model.spec.reciprocal_vectors()
    bulk, finite = [], []
    for k in selected:
        basis = EigenBasis.at(lattice.model, k)
        wave = lattice._plane_waves(k) @ basis.U
        bulk.append(full_matrix(variant, basis)[0])
        finite.append(wave.conj().T @ result.L_site @ wave)
    bulk, finite = np.asarray(bulk), np.asarray(finite)
    if band_pairs is not None:
        pairs = list(band_pairs)
        if not pairs:
            raise ValueError("band_pairs must not be empty")
        bulk = np.array([bulk[:, i, j] for i, j in pairs]).T
        finite = np.array([finite[:, i, j] for i, j in pairs]).T
    diffs = abs(bulk - finite)
    maximum = float(diffs.max())
    return {
        "max_abs_diff": maximum,
        "diffs": diffs,
        "bloch": bulk,
        "lattice": finite,
        "kpoints": selected,
        "passed": maximum <= tol,
    }


def compare_bloch_vs_lattice(model, variant, size, momenta=None, band=0, tol=1e-8):
    from ohk.model.finite import FiniteLattice

    if not 0 <= band < model.num_bands:
        raise ValueError("invalid band")
    return _compare(
        FiniteLattice(model, size),
        variant,
        momenta,
        [(band, j) for j in range(model.num_bands)],
        tol,
    )


def make_random_multiband_model(n_bands=3, n_hops=5, scale=0.35, seed=0):
    from ohk.model.analytic import AnalyticModel, Hopping
    from ohk.model.base import ModelSpec

    if n_bands < 1 or n_hops < 1 or scale <= 0:
        raise ValueError("band/hopping counts and scale must be positive")
    rng = np.random.default_rng(seed)
    directions = [
        (i + 1, 0) if i >= 5 else [(1, 0), (0, 1), (1, 1), (-1, 1), (2, 1)][i]
        for i in range(n_hops)
    ]
    hops = []
    for r in directions:
        t = (
            scale
            * (rng.normal(size=(n_bands, n_bands)) + 1j * rng.normal(size=(n_bands, n_bands)))
            / np.sqrt(n_bands)
        )
        hops.extend(Hopping(i, j, r, t[i, j]) for i in range(n_bands) for j in range(n_bands))
    return AnalyticModel(
        ModelSpec(n_bands, np.eye(2), name=f"random-{n_bands}-band"),
        hops,
        onsite=np.linspace(-1, 1, n_bands),
    )
