"""Cartesian BZ quadrature. Normalized weights sum to one; paths are not quadrature grids."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class BZGrid:
    kpoints: np.ndarray
    weights: np.ndarray
    normalized: bool
    lattice_vectors: np.ndarray
    labels: tuple[str, ...] = ()
    is_path: bool = False

    def __post_init__(self):
        from ohk.model.base import ModelSpec

        self.lattice_vectors = np.array(self.lattice_vectors, float, copy=True)
        d = len(self.lattice_vectors)
        ModelSpec(1, self.lattice_vectors, periodic_dim=d)
        self.kpoints = np.array(self.kpoints, float, copy=True)
        self.weights = np.array(self.weights, float, copy=True)
        if self.kpoints.ndim != 2 or self.kpoints.shape[1] != d or not len(self.kpoints):
            raise ValueError("kpoints must have shape (n_k>0, dimension)")
        if self.weights.shape != (len(self.kpoints),):
            raise ValueError("weights must have length n_k")
        if not np.isfinite(self.kpoints).all() or not np.isfinite(self.weights).all():
            raise ValueError("kpoints and weights must be finite")
        if np.any(self.weights < 0):
            raise ValueError("weights must be non-negative")
        target = 1.0 if self.normalized else self.reciprocal_volume()
        if not np.isclose(self.weights.sum(), target, rtol=1e-10, atol=1e-12):
            raise ValueError("weights must sum to 1 (normalized) or the Cartesian BZ volume")

    def __len__(self):
        return len(self.kpoints)

    @classmethod
    def uniform(cls, lattice_vectors, n, shift=0.0):
        a = np.asarray(lattice_vectors, float)
        d = len(a)
        counts = np.repeat(n, d) if np.isscalar(n) else np.asarray(n)
        if (
            counts.shape != (d,)
            or not np.isfinite(counts).all()
            or np.any(counts < 1)
            or np.any(counts != np.rint(counts))
        ):
            raise ValueError("n must contain positive integer mesh sizes")
        shift = np.broadcast_to(np.asarray(shift, float), (d,))
        if not np.isfinite(shift).all():
            raise ValueError("shift must be finite")
        frac = np.stack(
            np.meshgrid(
                *[(np.arange(int(q)) + s) / q for q, s in zip(counts, shift)], indexing="ij"
            ),
            axis=-1,
        )
        kpoints = frac.reshape(-1, d) @ (2 * np.pi * np.linalg.inv(a).T)
        return cls(kpoints, np.full(len(kpoints), 1 / len(kpoints)), True, a)

    @classmethod
    def path(cls, lattice_vectors, vertices, labels, points_per_segment=50):
        a, vertices = np.asarray(lattice_vectors, float), np.asarray(vertices, float)
        if vertices.ndim != 2 or vertices.shape[1] != len(a) or len(vertices) < 2:
            raise ValueError("vertices must contain at least two dimension-matched points")
        if (
            len(labels) != len(vertices)
            or not isinstance(points_per_segment, int)
            or points_per_segment < 2
        ):
            raise ValueError("one label per vertex and points_per_segment>=2 are required")
        frac = np.concatenate(
            [
                np.linspace(x, y, points_per_segment, endpoint=False)
                for x, y in zip(vertices[:-1], vertices[1:])
            ]
            + [vertices[-1:]]
        )
        k = frac @ (2 * np.pi * np.linalg.inv(a).T)
        return cls(k, np.full(len(k), 1 / len(k)), True, a, tuple(labels), True)

    def reciprocal_vectors(self):
        return 2 * np.pi * np.linalg.inv(self.lattice_vectors).T

    def check_reciprocal(self, tol=1e-12):
        error = float(
            np.max(
                abs(
                    self.lattice_vectors @ self.reciprocal_vectors().T / (2 * np.pi)
                    - np.eye(len(self.lattice_vectors))
                )
            )
        )
        if error > tol:
            raise ValueError(f"Reciprocal basis residual {error}")
        return error

    def reciprocal_volume(self):
        return (2 * np.pi) ** len(self.lattice_vectors) / self.cell_area()

    def cell_area(self):
        return float(abs(np.linalg.det(self.lattice_vectors)))
