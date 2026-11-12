"""Haldane benchmark with the complete d0 term and physical orbital embedding."""

from __future__ import annotations
import numpy as np
from ohk.model.analytic import AnalyticModel, Hopping
from ohk.model.base import ModelSpec


class HaldaneModel(AnalyticModel):
    def __init__(self, t1=1.0, t2=0.15, M=0.4, phi=0.0, nearest_neighbour_distance=1.0):
        if not np.isfinite([t1, t2, M, phi, nearest_neighbour_distance]).all():
            raise ValueError("Haldane parameters must be finite")
        if nearest_neighbour_distance <= 0:
            raise ValueError("nearest_neighbour_distance must be positive")
        self.t1, self.t2, self.M, self.phi = t1, t2, M, phi
        self.a = nearest_neighbour_distance
        spec = ModelSpec(
            2,
            self.lattice_vectors(self.a),
            num_orbitals_per_site=(1, 1),
            name="Haldane",
            metadata=dict(t1=t1, t2=t2, M=M, phi=phi),
        )
        hops = [Hopping(0, 1, r, t1) for r in [(0, 0), (0, -1), (-1, 0)]]
        for r in [(1, -1), (-1, 0), (0, 1)]:
            hops.extend(
                [Hopping(0, 0, r, t2 * np.exp(1j * phi)), Hopping(1, 1, r, t2 * np.exp(-1j * phi))]
            )
        super().__init__(
            spec,
            hops,
            onsite=np.array([M, -M]),
            orbital_positions=np.array([[0.0, 0.0], [self.a, 0.0]]),
        )

    @staticmethod
    def delta_vectors(a=1.0):
        return a * np.array([[1.0, 0.0], [-0.5, np.sqrt(3) / 2], [-0.5, -np.sqrt(3) / 2]])

    @staticmethod
    def lattice_vectors(a=1.0):
        d = HaldaneModel.delta_vectors(a)
        return np.array([d[0] - d[2], d[0] - d[1]])

    @staticmethod
    def next_nearest_vectors(a=1.0):
        a1, a2 = HaldaneModel.lattice_vectors(a)
        return np.array([a1 - a2, -a1, a2])

    def f(self, k):
        return np.exp(1j * (self.delta_vectors(self.a) @ self._momentum(k))).sum()

    def d0(self, k):
        return float(
            2
            * self.t2
            * np.cos(self.phi)
            * np.cos(self.next_nearest_vectors(self.a) @ self._momentum(k)).sum()
        )

    def dz(self, k):
        return float(
            self.M
            - 2
            * self.t2
            * np.sin(self.phi)
            * np.sin(self.next_nearest_vectors(self.a) @ self._momentum(k)).sum()
        )

    def chern_number(self, n_k=64, band=0):
        """FHS Chern number in Cartesian orientation and a periodic orbital gauge."""
        if not isinstance(n_k, (int, np.integer)) or n_k < 3 or band not in (0, 1):
            raise ValueError("n_k must be >=3 and band must be 0 or 1")
        reciprocal = self.spec.reciprocal_vectors()
        states = np.empty((n_k, n_k, 2), complex)
        for i in range(n_k):
            for j in range(n_k):
                k = np.array([i, j]) @ reciprocal / n_k
                eps, u = np.linalg.eigh(self.H(k))
                if eps[1] - eps[0] < 1e-10:
                    raise ValueError("Chern number undefined at a band touching")
                states[i, j] = np.exp(1j * (self.orbital_positions @ k)) * u[:, band]
        links = []
        for axis in (0, 1):
            link = np.sum(states.conj() * np.roll(states, -1, axis=axis), axis=-1)
            if np.any(abs(link) < 1e-12):
                raise ValueError("Singular FHS links; increase n_k")
            links.append(link / abs(link))
        x, y = links
        flux = np.angle(x * np.roll(y, -1, axis=0) * np.roll(x.conj(), -1, axis=1) * y.conj()).sum()
        # Overlap links have the opposite sign to A=i<u|du>.
        return int(np.rint(-flux / (2 * np.pi) * np.sign(np.linalg.det(self.lattice()))))

    def phase_boundary(self):
        if self.t2 == 0 or abs(self.M) > 3 * np.sqrt(3) * abs(self.t2):
            return float("nan")
        return float(np.arcsin(abs(self.M) / (3 * np.sqrt(3) * abs(self.t2))))
