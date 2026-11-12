"""Eq. (16a), using zero diagonal derivative coefficients (parallel-transport gauge)."""

from __future__ import annotations
import numpy as np
from ohk.oam.base import OAMVariant


class CovariantProjector(OAMVariant):
    name = "covariant-16a"
    reference = "Gautam & Satpathy, PRB 114, 214403 (2026), Eq. (16a)"
    gauge_invariant = True  # Historical flag: gauge-covariant matrix, invariant observables.
    diagonal = False

    def __init__(self, include_diagonal=False):
        self.include_diagonal = include_diagonal
        self.diagonal = bool(include_diagonal)

    def matrix(self, basis):
        result = self._projected_part(basis) + self._berry_part(basis)
        if not self.include_diagonal:
            np.fill_diagonal(result, 0)
        return result

    @staticmethod
    def _projected_part(basis):
        cx, cy = basis.C[:2]
        mean = (basis.eps[:, None] + basis.eps[None, :]) / 2
        xy = cx.conj().T @ (basis.eps[:, None] * cy) - mean * (cx.conj().T @ cy)
        yx = cy.conj().T @ (basis.eps[:, None] * cx) - mean * (cy.conj().T @ cx)
        return 1j * (xy - yx)

    @staticmethod
    def _berry_part(basis):
        a = 1j * basis.C
        velocity_diag = np.diagonal(basis.v, axis1=-2, axis2=-1)
        vsum = velocity_diag[:, :, None] + velocity_diag[:, None, :]
        return 0.5 * (a[0] * vsum[1] - a[1] * vsum[0])
