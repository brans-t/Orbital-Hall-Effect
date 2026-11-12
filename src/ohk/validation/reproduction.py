"""Batched Bloch and sparse finite-PBC calculations for reproducible figure sweeps.

The Bloch route evaluates hopping Fourier sums and the projected Eq. (16a).
The finite route builds full sparse real-space operators, Fourier-reduces their origin rows
at the allowed discrete momenta, and uses the direct position/velocity definition.
No artificial broadening or adjustable normalization is applied.
"""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix

from ohk.eigen import DegeneracyError
from ohk.transport.grids import BZGrid


@dataclass
class BandResponse:
    kpoints: np.ndarray
    energies: np.ndarray
    oam: np.ndarray
    kernels: np.ndarray  # Per-band response in the package's e/(2*pi) internal units.
    min_gap: float

    def sigma(self, mu=None):
        """mu=None explicitly occupies only band zero (half filling for a two-band model)."""
        if mu is None:
            return float(self.kernels[:, 0].mean())
        if not np.isfinite(mu):
            raise ValueError("mu must be finite")
        return float(np.sum(self.kernels * (self.energies < mu)) / len(self.kpoints))

    def fermi_scan(self, mus):
        """Exact zero-temperature occupation, via a sorted cumulative response."""
        mus = np.asarray(mus, float)
        if mus.ndim != 1 or not np.isfinite(mus).all():
            raise ValueError("mus must be a finite vector")
        order = np.argsort(self.energies.ravel())
        eps = self.energies.ravel()[order]
        terms = self.kernels.ravel()[order] / len(self.kpoints)
        cumulative = np.r_[0.0, np.cumsum(terms)]
        return cumulative[np.searchsorted(eps, mus, side="left")]


def _response(kpoints, h, dh, cell_area, method):
    eps, u = np.linalg.eigh(h)
    dagger = u.conj().swapaxes(-1, -2)
    v = dagger[:, None] @ dh @ u[:, None]
    gap = eps[:, :, None] - eps[:, None, :]
    off = ~np.eye(eps.shape[1], dtype=bool)
    minimum = float(np.min(np.abs(gap[:, off])))
    if minimum < 1e-14:
        raise DegeneracyError("Figure sweep hits an exact band touching")
    inv = np.zeros_like(gap)
    inv[:, off] = 1 / gap[:, off]
    vx, vy = v[:, 0], v[:, 1]
    if method == "projector":
        c = -v * inv[:, None]
        cx, cy = c[:, 0], c[:, 1]
        cxdag, cydag = cx.conj().swapaxes(-1, -2), cy.conj().swapaxes(-1, -2)
        mean = (eps[:, :, None] + eps[:, None, :]) / 2
        projected = 1j * (
            cxdag @ (eps[:, :, None] * cy)
            - mean * (cxdag @ cy)
            - cydag @ (eps[:, :, None] * cx)
            + mean * (cydag @ cx)
        )
        diagonal = np.diagonal(v, axis1=-2, axis2=-1)
        vsum = diagonal[:, :, :, None] + diagonal[:, :, None, :]
        berry = 0.5j * (c[:, 0] * vsum[:, 1] - c[:, 1] * vsum[:, 0])
        angular = projected + berry
    elif method == "direct":
        r = -1j * v * inv[:, None]
        rx, ry = r[:, 0], r[:, 1]
        angular = 0.5 * (rx @ vy - ry @ vx - vx @ ry + vy @ rx)
    else:
        raise ValueError("method must be projector or direct")
    current = 0.5 * (vx @ angular + angular @ vx)
    kernels = (4 * np.pi / cell_area) * np.imag(current * vy.swapaxes(-1, -2) * inv**2).sum(axis=-1)
    return BandResponse(np.asarray(kpoints), eps, angular, kernels, minimum)


def bloch_response(model, kpoints):
    """Batched analytic hopping sums; supports explicit embedded hopping models."""
    kpoints = np.asarray(kpoints, float)
    d, nb = model.spec.periodic_dim, model.num_bands
    if kpoints.ndim != 2 or kpoints.shape[1] != d or not np.isfinite(kpoints).all():
        raise ValueError("kpoints must be a finite Cartesian array")
    h = np.zeros((len(kpoints), nb, nb), complex)
    dh = np.zeros((len(kpoints), d, nb, nb), complex)
    tau = model.orbital_positions
    for hop in model.hoppings:
        displacement = np.asarray(hop.dR) @ model.lattice() + tau[hop.j] - tau[hop.i]
        term = hop.amplitude * np.exp(1j * (kpoints @ displacement))
        h[:, hop.i, hop.j] += term
        dh[:, :, hop.i, hop.j] += 1j * term[:, None] * displacement
    return _response(kpoints, h, dh, abs(np.linalg.det(model.lattice())), "projector")


def finite_pbc_response(model, size=60):
    """Independent sparse N^d torus, exactly reduced into its translation sectors.

    Stores O(N^d * hopping_count), rather than dense O((N^d * num_bands)^2) matrices.
    Fourier reduction is exact at the finite torus's allowed momenta, not interpolation.
    Only H(R) hoppings, their physical displacements and lattice metadata are consumed.
    """
    d, nb = model.spec.periodic_dim, model.num_bands
    sizes = (size,) * d if np.isscalar(size) else tuple(size)
    if len(sizes) != d or any(int(n) != n or n < 1 for n in sizes):
        raise ValueError("size must contain positive integers")
    sizes = tuple(map(int, sizes))
    cells = np.indices(sizes).reshape(d, -1).T
    count, sites = len(cells), len(cells) * nb
    tau = model.orbital_positions
    rows, columns, values, velocity_values = [], [], [], [[] for _ in range(d)]
    for hop in model.hoppings:
        target = (cells + np.asarray(hop.dR, int)) % sizes
        target_index = np.ravel_multi_index(target.T, sizes)
        rows.append(np.arange(count) * nb + hop.i)
        columns.append(target_index * nb + hop.j)
        amplitudes = np.full(count, hop.amplitude, complex)
        values.append(amplitudes)
        delta = np.asarray(hop.dR) @ model.lattice() + tau[hop.j] - tau[hop.i]
        for j in range(d):
            velocity_values[j].append(1j * delta[j] * amplitudes)
    row, col = np.concatenate(rows), np.concatenate(columns)

    def reduce_operator(elements):
        operator = coo_matrix((np.concatenate(elements), (row, col)), shape=(sites, sites)).tocsr()
        difference = operator - operator.conj().T
        if difference.nnz and np.max(abs(difference.data)) > 1e-12:
            raise ValueError("Non-Hermitian finite real-space operator")
        origin = operator[:nb].toarray().reshape(nb, count, nb).transpose(1, 0, 2)
        cell_operator = origin.reshape(*sizes, nb, nb)
        return (np.fft.ifftn(cell_operator, axes=tuple(range(d))) * count).reshape(count, nb, nb)

    h = reduce_operator(values)
    dh = np.stack([reduce_operator(elements) for elements in velocity_values], axis=1)
    kpoints = (cells / np.asarray(sizes)) @ model.spec.reciprocal_vectors()
    embedding = np.exp(1j * (kpoints @ tau.T))
    phase = embedding.conj()[:, :, None] * embedding[:, None, :]
    h *= phase
    dh *= phase[:, None]
    return _response(kpoints, h, dh, abs(np.linalg.det(model.lattice())), "direct")


def uniform_response(model, n=240):
    return bloch_response(model, BZGrid.uniform(model.lattice(), n).kpoints)
