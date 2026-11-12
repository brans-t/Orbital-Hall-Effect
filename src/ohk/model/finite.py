"""Independent finite-PBC oracle from real-space hopping matrices.

Velocity uses the UNWRAPPED displacement of each hopping, including boundary winding. A naive
commutator with sawtooth site coordinates is not a periodic bulk velocity. Position diagonals
are removed in the ENERGY basis, not the site basis (where doing so would erase position).

Translations resolve the momentum sectors before diagonalization, avoiding arbitrary mixing
of equal-energy states at different momenta. No Bloch H, dH or production OAM formula is used.
"""

from __future__ import annotations
from dataclasses import dataclass
from itertools import product
import numpy as np


def oam_from_velocity(energies, velocity):
    """Direct symmetrized position/velocity definition, within one momentum sector."""
    gap = energies[:, None] - energies[None, :]
    mask = ~np.eye(len(energies), dtype=bool)
    if np.any(abs(gap[mask]) < 1e-14):
        from ohk.eigen import DegeneracyError

        raise DegeneracyError("Degenerate finite-lattice momentum sector")
    position = np.zeros_like(velocity, dtype=complex)
    position[:, mask] = -1j * velocity[:, mask] / gap[mask]
    x, y = position[:2]
    vx, vy = velocity[:2]
    return 0.5 * (x @ vy - y @ vx - vx @ y + vy @ x)


@dataclass
class FiniteLatticeResult:
    energies: np.ndarray
    eigenvectors: np.ndarray
    L: np.ndarray
    positions: np.ndarray
    L_site: np.ndarray
    momenta: np.ndarray | None = None

    def at_momentum_indices(self, i, j):
        """2x2 block for two energy-sorted state indices (not momentum indices)."""
        return self.L[np.ix_([i, j], [i, j])]


class FiniteLattice:
    def __init__(self, model, size):
        if not hasattr(model, "hoppings"):
            raise TypeError("FiniteLattice requires an explicit real-space hopping list")
        self.model = model
        d = model.spec.periodic_dim
        size = (size,) * d if np.isscalar(size) else tuple(size)
        if len(size) != d or any(int(n) != n or n < 1 for n in size):
            raise ValueError("size must contain positive integers for every periodic dimension")
        self.size = tuple(map(int, size))
        self.cells = np.array(list(product(*[range(n) for n in self.size])), dtype=int)
        self.n_cells = len(self.cells)
        self.n_sites = self.n_cells * model.num_bands
        self._operators = None

    def _build_operators(self):
        if self._operators is not None:
            return self._operators
        model, n = self.model, self.n_sites
        h = np.zeros((n, n), complex)
        v = np.zeros((model.spec.periodic_dim, n, n), complex)
        tau = getattr(
            model, "orbital_positions", np.zeros((model.num_bands, model.spec.periodic_dim))
        )
        for index, cell in enumerate(self.cells):
            for hop in model.hoppings:
                target = (cell + np.array(hop.dR, int)) % self.size
                target_cell = np.ravel_multi_index(tuple(target), self.size)
                i, j = index * model.num_bands + hop.i, target_cell * model.num_bands + hop.j
                displacement = np.array(hop.dR) @ model.lattice() + tau[hop.j] - tau[hop.i]
                h[i, j] += hop.amplitude
                v[:, i, j] += 1j * displacement * hop.amplitude
        for value in (h, v):
            if not np.allclose(value, value.conj().swapaxes(-1, -2), atol=1e-12, rtol=0):
                raise ValueError("Finite-lattice hopping matrices are not Hermitian")
        self._operators = h, v
        return h, v

    def build_hamiltonian(self):
        return self._build_operators()[0].copy()

    def site_positions(self):
        tau = getattr(
            self.model,
            "orbital_positions",
            np.zeros((self.model.num_bands, self.model.spec.periodic_dim)),
        )
        return (self.cells[:, None, :] @ self.model.lattice() + tau[None, :, :]).reshape(
            self.n_sites, self.model.spec.periodic_dim
        )

    def allowed_momenta(self):
        return (self.cells / np.array(self.size)) @ self.model.spec.reciprocal_vectors()

    def _plane_waves(self, k):
        phase = np.exp(1j * (self.site_positions() @ k)) / np.sqrt(self.n_cells)
        q = np.zeros((self.n_sites, self.model.num_bands), complex)
        q[np.arange(self.n_sites), np.arange(self.n_sites) % self.model.num_bands] = phase
        return q

    @staticmethod
    def velocity(H, positions):
        """Open-boundary -i[r,H]; periodic calculations must use hopping windings instead."""
        positions = np.asarray(positions)
        return tuple(
            -1j * (positions[:, j, None] - positions[None, :, j]) * H
            for j in range(positions.shape[1])
        )

    @staticmethod
    def oam_operator(H, positions, drop_diagonal=True):
        """Open-boundary OAM with zero energy-basis position diagonal."""
        if not drop_diagonal:
            raise ValueError("The energy-basis position diagonal must be removed")
        energies, u = np.linalg.eigh(H)
        velocity = np.array(FiniteLattice.velocity(H, positions))
        L = oam_from_velocity(energies, u.conj().T @ velocity @ u)
        return u @ L @ u.conj().T

    def compute(self):
        h, v = self._build_operators()
        energies, states, blocks = [], [], []
        n = self.model.num_bands
        momenta = self.allowed_momenta()
        for k in momenta:
            q = self._plane_waves(k)
            eps, u = np.linalg.eigh(q.conj().T @ h @ q)
            wave = q @ u
            velocity = wave.conj().T @ v @ wave
            blocks.append(oam_from_velocity(eps, velocity))
            energies.extend(eps)
            states.append(wave)
        states = np.concatenate(states, axis=1)
        L = np.zeros((self.n_sites,) * 2, complex)
        for index, block in enumerate(blocks):
            L[index * n : (index + 1) * n, index * n : (index + 1) * n] = block
        L_site = states @ L @ states.conj().T
        order = np.argsort(energies, kind="stable")
        return FiniteLatticeResult(
            np.array(energies)[order],
            states[:, order],
            L[np.ix_(order, order)],
            self.site_positions(),
            L_site,
            momenta,
        )

    def compare_with_bloch(self, model, band_pairs, tol=1e-8):
        if model is not self.model:
            raise ValueError("Comparison must use the model that built this finite lattice")
        from ohk.oam.covariant import CovariantProjector
        from ohk.validation.oracle import _compare

        return _compare(self, CovariantProjector(), band_pairs=band_pairs, tol=tol)
