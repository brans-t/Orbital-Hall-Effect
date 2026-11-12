"""Orbital Hall kit: tight-binding OAM, transport and independent numerical validation.

Internally hbar=m_e=e=1; k is Cartesian and velocity indices are [direction, band, band].
Off-diagonal matrices are gauge covariant; diagonal expectations and Kubo contractions are
gauge invariant. See docs/design/05-implementation-notes.md for audited conventions.
"""

from ohk.units import HBAR, M_E, ELECTRON_CHARGE
from ohk.eigen import EigenBasis, DegeneracyError
from ohk.model import HaldaneModel, AnalyticModel, Hopping, ModelSpec, Wannier90Model
from ohk.oam import CovariantProjector, SpectralSum, full_matrix
from ohk.transport import BZGrid, orbital_hall_conductivity, fermi_scan

__version__ = "0.1.0"
__all__ = [
    "__version__",
    "HBAR",
    "M_E",
    "ELECTRON_CHARGE",
    "EigenBasis",
    "DegeneracyError",
    "HaldaneModel",
    "AnalyticModel",
    "Hopping",
    "ModelSpec",
    "Wannier90Model",
    "CovariantProjector",
    "SpectralSum",
    "full_matrix",
    "BZGrid",
    "orbital_hall_conductivity",
    "fermi_scan",
]
