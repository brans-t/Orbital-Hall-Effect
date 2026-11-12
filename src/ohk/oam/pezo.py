"""Truncated OAM: retain only the intermediate-band sum in Eq. (16b).

This is the FIRST term as ordered in the paper, not the diagonal-velocity term mistakenly
identified in the original design. Thus off-diagonal elements vanish for every two-band model.
"""

from __future__ import annotations
import numpy as np
from ohk.oam.base import OAMVariant
from ohk.eigen import inverse_gaps


class PezoTruncated(OAMVariant):
    name = "pezo-truncated"
    reference = "Pezo, Ovalle, Manchon, PRB 106, 104414 (2022); truncated comparison"
    gauge_invariant = True
    diagonal = False

    def __init__(self, include_diagonal=False):
        self.include_diagonal = include_diagonal
        self.diagonal = bool(include_diagonal)

    def matrix(self, basis):
        inv = inverse_gaps(basis.eps, basis.regularization_eta)
        vx, vy = basis.v[:2]
        result = np.zeros((basis.n_bands,) * 2, complex)
        for n in range(basis.n_bands):
            for p in range(basis.n_bands):
                if n == p:
                    continue
                for m in range(basis.n_bands):
                    if m != n and m != p:
                        result[n, p] -= (
                            0.5j
                            * (inv[n, m] + inv[p, m])
                            * (vx[n, m] * vy[m, p] - vy[n, m] * vx[m, p])
                        )
        if self.include_diagonal:
            from ohk.oam.assemble import diagonal_closed_form

            np.fill_diagonal(result, diagonal_closed_form(basis))
        return result

    @staticmethod
    def is_blind_to_sum_term(n_bands):
        return n_bands == 2
