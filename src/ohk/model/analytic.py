"""Hopping model: H_ij(k) = sum_R t_ij(R) exp(+ik.(R+tau_j-tau_i)).

Only missing conjugate partners are added. Orbital positions are Cartesian and zero by default.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ohk.model.base import TightBindingModel, ModelSpec


@dataclass(frozen=True)
class Hopping:
    i: int
    j: int
    dR: tuple[int, ...]
    amplitude: complex

    def __post_init__(self):
        if any(not isinstance(i, (int, np.integer)) or i < 0 for i in (self.i, self.j)):
            raise ValueError("orbital indices must be non-negative integers")
        if len(self.dR) not in (2, 3) or any(int(x) != x for x in self.dR):
            raise ValueError("dR must contain two or three integer components")
        if not np.isfinite(self.amplitude):
            raise ValueError("amplitude must be finite")


class AnalyticModel(TightBindingModel):
    def __init__(
        self, spec: ModelSpec, hoppings, onsite=None, assume_hermitian=False, orbital_positions=None
    ):
        self.spec = spec
        shape = (spec.num_bands, spec.periodic_dim)
        self.orbital_positions = (
            np.zeros(shape)
            if orbital_positions is None
            else np.array(orbital_positions, dtype=float, copy=True)
        )
        if self.orbital_positions.shape != shape or not np.isfinite(self.orbital_positions).all():
            raise ValueError("orbital_positions has an invalid shape or non-finite values")
        terms = {}
        for hop in hoppings:
            if max(hop.i, hop.j) >= spec.num_bands or len(hop.dR) != spec.periodic_dim:
                raise ValueError("hopping incompatible with ModelSpec")
            key = (hop.i, hop.j, tuple(hop.dR))
            terms[key] = terms.get(key, 0j) + hop.amplitude
        if not assume_hermitian:
            for (i, j, r), amplitude in list(terms.items()):
                partner = (j, i, tuple(-x for x in r))
                if partner not in terms:
                    terms[partner] = np.conj(amplitude)
                elif not np.isclose(terms[partner], np.conj(amplitude), atol=1e-12, rtol=1e-12):
                    raise ValueError("inconsistent conjugate hopping partners")
        if onsite is not None:
            onsite = np.asarray(onsite)
            if onsite.shape != (spec.num_bands,) or not np.isfinite(onsite).all():
                raise ValueError("onsite must contain num_bands finite real energies")
            if np.any(np.abs(onsite.imag) > 1e-14):
                raise ValueError("onsite energies must be real")
            for i, energy in enumerate(onsite):
                key = (i, i, (0,) * spec.periodic_dim)
                terms[key] = terms.get(key, 0j) + energy
        self.hoppings = tuple(Hopping(i, j, r, complex(t)) for (i, j, r), t in terms.items())
        self._i = np.array([h.i for h in self.hoppings], dtype=int)
        self._j = np.array([h.j for h in self.hoppings], dtype=int)
        self._t = np.array([h.amplitude for h in self.hoppings], dtype=complex)
        cells = np.array([h.dR for h in self.hoppings], dtype=float).reshape(-1, spec.periodic_dim)
        self._displacements = (
            cells @ self.lattice()
            + self.orbital_positions[self._j]
            - self.orbital_positions[self._i]
        )

    def H(self, k):
        values = self._t * np.exp(1j * (self._displacements @ self._momentum(k)))
        h = np.zeros((self.num_bands,) * 2, dtype=complex)
        np.add.at(h, (self._i, self._j), values)
        return h

    def dH(self, k):
        values = self._t * np.exp(1j * (self._displacements @ self._momentum(k)))
        dh = np.zeros((self.spec.periodic_dim, self.num_bands, self.num_bands), dtype=complex)
        for direction in range(self.spec.periodic_dim):
            np.add.at(
                dh[direction], (self._i, self._j), 1j * self._displacements[:, direction] * values
            )
        return dh

    def hermitize(self):
        return AnalyticModel(self.spec, self.hoppings, orbital_positions=self.orbital_positions)
