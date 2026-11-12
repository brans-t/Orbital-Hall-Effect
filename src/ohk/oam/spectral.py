"""Eq. (16b), with denominators eps_n-eps_m and therefore a minus sum sign.

The paper prints a PLUS with the OPPOSITE denominators eps_m-eps_n. Both are equivalent.
sign=+1 below reproduces the erroneous transcription in the original project design notes,
not the published paper. It exists only as a negative control.
"""

from __future__ import annotations
import numpy as np
from ohk.eigen import inverse_gaps
from ohk.oam.base import OAMVariant


class SpectralSum(OAMVariant):
    name = "spectral-16b"
    reference = "Gautam & Satpathy, PRB 114, 214403 (2026), Eq. (16b)"
    gauge_invariant = True
    diagonal = True

    def __init__(self, sign=-1):
        if sign not in (-1, 1):
            raise ValueError("sign must be -1 (paper) or +1 (incorrect-transcription control)")
        self.sign = sign
        if sign == 1:
            self.name = "spectral-wrong-sign-control"

    def matrix(self, basis):
        result = self.off_diagonal(basis, self.sign)
        np.fill_diagonal(result, self.diagonal_elements(basis))
        return result

    @staticmethod
    def off_diagonal(basis, sign=-1):
        if sign not in (-1, 1):
            raise ValueError("sign must be -1 or +1")
        inv = inverse_gaps(basis.eps, basis.regularization_eta)
        vx, vy = basis.v[:2]
        result = np.zeros((basis.n_bands,) * 2, complex)
        for n in range(basis.n_bands):
            for p in range(basis.n_bands):
                if n == p:
                    continue
                first = ((vx[n, n] + vx[p, p]) * vy[n, p] - (vy[n, n] + vy[p, p]) * vx[n, p]) * inv[
                    n, p
                ]
                total = 0j
                for m in range(basis.n_bands):
                    if m != n and m != p:
                        total += (inv[n, m] + inv[p, m]) * (
                            vx[n, m] * vy[m, p] - vy[n, m] * vx[m, p]
                        )
                result[n, p] = 0.5j * (first + sign * total)
        return result

    @staticmethod
    def diagonal_elements(basis):
        inv = inverse_gaps(basis.eps, basis.regularization_eta)
        return 2 * np.sum(np.imag(basis.v[0] * basis.v[1].conj()) * inv, axis=1)
