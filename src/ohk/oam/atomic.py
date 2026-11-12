"""Atomic-centered angular momentum in explicitly specified local multiplets."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ohk.oam.base import OAMVariant


@dataclass(frozen=True)
class ACASpec:
    """orbital_l contains ONE ell per local multiplet; a multiplet occupies 2l+1 bands.

    For a p multiplet use (1,), not (1,1,1). The basis is |ell,m>, m=-ell,...,+ell.
    local_matrices may instead specify custom blocks in the actual orbital basis.
    """

    orbital_l: tuple[int, ...]
    quantized_axis: str = "z"
    local_matrices: tuple[np.ndarray, ...] | None = None

    def __post_init__(self):
        if self.quantized_axis not in ("x", "y", "z"):
            raise ValueError("quantized_axis must be 'x', 'y' or 'z'")
        if not self.orbital_l or any(int(ell) != ell or ell < 0 for ell in self.orbital_l):
            raise ValueError("orbital_l must contain non-negative integers")


class AtomCentered(OAMVariant):
    name = "aca"
    reference = "atom-centered approximation"
    gauge_invariant = True
    diagonal = True

    def __init__(self, spec: ACASpec):
        self.spec = spec

    def matrix(self, basis):
        if self.spec.local_matrices is not None:
            blocks = [np.asarray(b, complex) for b in self.spec.local_matrices]
        else:
            blocks = []
            for ell in self.spec.orbital_l:
                m = np.arange(-ell, ell + 1)
                plus = np.diag(np.sqrt(ell * (ell + 1) - m[:-1] * (m[:-1] + 1)), -1)
                blocks.append(
                    {
                        "z": self.local_lz(ell),
                        "x": (plus + plus.T) / 2,
                        "y": (plus - plus.T) / (2j),
                    }[self.spec.quantized_axis]
                )
        local = np.zeros((basis.n_bands,) * 2, complex)
        offset = 0
        for block in blocks:
            if block.ndim != 2 or block.shape[0] != block.shape[1]:
                raise ValueError("local angular-momentum blocks must be square")
            n = len(block)
            if not np.isfinite(block).all() or not np.allclose(block, block.conj().T):
                raise ValueError("local angular-momentum blocks must be finite and Hermitian")
            if offset + n > basis.n_bands:
                raise ValueError("local multiplets exceed num_bands")
            local[offset : offset + n, offset : offset + n] = block
            offset += n
        if offset != basis.n_bands:
            raise ValueError("local multiplet dimensions must sum to num_bands")
        return basis.U.conj().T @ local @ basis.U

    @staticmethod
    def local_lz(ell):
        if not isinstance(ell, (int, np.integer)) or ell < 0:
            raise ValueError("ell must be a non-negative integer")
        return np.diag(np.arange(-ell, ell + 1, dtype=float))
