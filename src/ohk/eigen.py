"""L1 eigenbasis: O -> D^dag O D under band phases; observables are invariant."""

from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
import numpy as np


class DegeneracyError(ValueError):
    """The isolated-band expression is undefined at an exact band degeneracy."""


def inverse_gaps(eps, eta=0.0):
    eps = np.asarray(eps, dtype=float)
    if eps.ndim != 1 or not np.isfinite(eps).all() or not np.isfinite(eta) or eta < 0:
        raise ValueError("eps must be a finite vector and eta must be non-negative")
    gap = eps[:, None] - eps[None, :]
    off = ~np.eye(len(eps), dtype=bool)
    if eta == 0 and np.any(np.abs(gap[off]) < 1e-14):
        raise DegeneracyError("Exact band degeneracy: isolated-band OAM is undefined")
    result = np.zeros_like(gap)
    if eta:
        result[off] = gap[off] / (gap[off] ** 2 + eta**2)
    else:
        result[off] = 1 / gap[off]
    return result


def velocity_matrix(H, dH, U):
    h, dh, u = np.asarray(H), np.asarray(dH), np.asarray(U)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or u.shape != h.shape:
        raise ValueError("H and U must be matching square matrices")
    if dh.ndim != 3 or dh.shape[1:] != h.shape:
        raise ValueError("dH must have shape (dimension, num_bands, num_bands)")
    return u.conj().T @ dh @ u


def derivative_coefficients(eps, v, eta=0.0):
    return -np.asarray(v) * inverse_gaps(eps, eta)[None, :, :]


@dataclass
class EigenBasis:
    k: np.ndarray
    eps: np.ndarray
    U: np.ndarray
    v: np.ndarray
    dH: np.ndarray | None = None
    regularization_eta: float = 0.0

    @property
    def n_bands(self):
        return len(self.eps)

    @classmethod
    def at(cls, model, k):
        k = model._momentum(k)
        model.validate(k)
        h, dh = model.H(k), model.dH(k)
        eps, u = np.linalg.eigh(h)
        return cls(k.copy(), eps, u, velocity_matrix(h, dh, u), dh)

    @cached_property
    def C(self):
        return derivative_coefficients(self.eps, self.v, self.regularization_eta)

    def derivative_coefficients(self):
        return self.C

    def berry_connection(self, derivative_vectors=None):
        """Diagonal i<u|du>, requiring explicit derivatives in a specified smooth gauge."""
        if derivative_vectors is None:
            raise ValueError("Berry connection requires explicit gauge-fixed derivative_vectors")
        du = np.asarray(derivative_vectors)
        if du.shape != self.v.shape:
            raise ValueError("derivative_vectors must have shape (dimension, num_bands, num_bands)")
        return 1j * np.einsum("an,jan->jn", self.U.conj(), du)

    def min_gap(self):
        return float(np.diff(self.eps).min()) if self.n_bands > 1 else float("inf")

    def degenerate_mask(self, tol=1e-8):
        if tol <= 0:
            raise ValueError("tol must be positive")
        mask = np.zeros(self.n_bands, bool)
        close = np.diff(self.eps) < tol
        mask[:-1] |= close
        mask[1:] |= close
        return mask

    def check(self, tol=1e-10):
        diagnostics = {
            "unitarity": float(np.max(abs(self.U.conj().T @ self.U - np.eye(self.n_bands)))),
            "velocity_hermiticity": float(np.max(abs(self.v - self.v.conj().swapaxes(-1, -2)))),
            "velocity_diagonal_imag": float(
                np.max(abs(np.diagonal(self.v, axis1=-2, axis2=-1).imag))
            ),
        }
        if not np.isfinite(self.C).all() or any(x > tol for x in diagnostics.values()):
            raise ValueError(f"Eigenbasis invariant failure: {diagnostics}")
        return diagnostics
