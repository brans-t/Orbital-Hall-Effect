"""L0 contract: Cartesian momenta and derivatives; primitive vectors are rows."""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import warnings
import numpy as np


@dataclass(frozen=True)
class ModelSpec:
    num_bands: int
    lattice_vectors: np.ndarray
    periodic_dim: int = 2
    num_orbitals_per_site: tuple[int, ...] = ()
    name: str = ""
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.num_bands, (int, np.integer)) or self.num_bands < 1:
            raise ValueError("num_bands must be a positive integer")
        a = np.array(self.lattice_vectors, dtype=float, copy=True)
        if self.periodic_dim not in (2, 3) or a.shape != (self.periodic_dim,) * 2:
            raise ValueError("lattice_vectors must have shape (periodic_dim, periodic_dim)")
        if not np.isfinite(a).all() or abs(np.linalg.det(a)) < 1e-14:
            raise ValueError("lattice_vectors must be finite and linearly independent")
        a.setflags(write=False)
        object.__setattr__(self, "lattice_vectors", a)

    def reciprocal_vectors(self):
        return 2 * np.pi * np.linalg.inv(self.lattice_vectors).T


class TightBindingModel(ABC):
    spec: ModelSpec

    @abstractmethod
    def H(self, k):
        """Hermitian Hamiltonian, shape (num_bands, num_bands)."""

    @property
    def num_bands(self):
        return self.spec.num_bands

    def lattice(self):
        return self.spec.lattice_vectors.copy()

    def _momentum(self, k):
        k = np.asarray(k, dtype=float)
        if k.shape != (self.spec.periodic_dim,) or not np.isfinite(k).all():
            raise ValueError("k must be a finite Cartesian vector of dimension periodic_dim")
        return k

    def dH(self, k):
        if not getattr(self, "_warned_fd", False):
            warnings.warn(
                "Using finite-difference dH; prefer an analytic derivative",
                UserWarning,
                stacklevel=2,
            )
            self._warned_fd = True
        return self.dH_finite_difference(k)

    def dH_finite_difference(self, k, step=1e-5):
        k = self._momentum(k)
        if not np.isfinite(step) or step <= 0:
            raise ValueError("step must be positive and finite")
        offsets = np.eye(len(k)) * step
        return np.array([(self.H(k + q) - self.H(k - q)) / (2 * step) for q in offsets])

    def energies(self, k):
        return np.linalg.eigvalsh(self.H(k))

    def is_hermitian(self, k, tol=1e-10):
        h = self.H(k)
        return bool(np.allclose(h, h.conj().T, atol=tol, rtol=0))

    def validate(self, k=None):
        k = self._momentum(np.zeros(self.spec.periodic_dim) if k is None else k)
        h, dh = np.asarray(self.H(k)), np.asarray(self.dH(k))
        shape = (self.num_bands,) * 2
        if h.shape != shape or dh.shape != (self.spec.periodic_dim, *shape):
            raise ValueError("H or dH has an invalid shape")
        for value in (h, dh):
            if not np.isfinite(value).all():
                raise ValueError("H and dH must be finite")
            if not np.allclose(value, value.conj().swapaxes(-1, -2), atol=1e-10, rtol=0):
                raise ValueError("H and dH must be Hermitian")
