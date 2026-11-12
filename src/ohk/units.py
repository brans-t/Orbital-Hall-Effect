"""Output-only physical conversion for Hamiltonian energies in eV and lengths in Angstrom.

Internal m_e=hbar=e=1 does not make an arbitrary eV/Angstrom model atomic-unit calibrated.
L/hbar acquires the factor m_e*E0*a0^2/hbar^2. OHC has units of charge in 2D,
not Siemens: it describes angular-momentum current, not electrical current.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

HBAR = 1.054_571_817e-34
M_E = 9.109_383_7015e-31
ELECTRON_CHARGE = 1.602_176_634e-19
BOHR_MAGNETON = 9.274_010_0783e-24
hbar = m_e = e = 1.0


@dataclass(frozen=True)
class UnitSystem:
    hbar: float = HBAR
    energy: float = 1.0  # eV
    length: float = 1.0  # Angstrom

    def __post_init__(self):
        if (
            not np.isfinite([self.hbar, self.energy, self.length]).all()
            or min(self.hbar, self.energy, self.length) <= 0
        ):
            raise ValueError("unit scales must be positive and finite")

    @property
    def oam_scale(self):
        return M_E * self.energy * ELECTRON_CHARGE * (self.length * 1e-10) ** 2 / self.hbar**2

    def oam_to_hbar(self, L_internal):
        """Dimensionless physical L/hbar, including the declared model scales."""
        return np.asarray(L_internal) * self.oam_scale

    def oam_to_si(self, L_internal):
        """Angular momentum in J*s."""
        return self.oam_to_hbar(L_internal) * self.hbar

    def conductivity_to_si(self, sigma_internal, dimension=2):
        """OHC in C (2D) or C/m (3D), from coefficients in e/(2*pi[*length])."""
        if dimension not in (2, 3):
            raise ValueError("dimension must be 2 or 3")
        return (
            np.asarray(sigma_internal)
            * self.oam_scale
            * ELECTRON_CHARGE
            / (2 * np.pi)
            / (self.length * 1e-10) ** (dimension - 2)
        )
