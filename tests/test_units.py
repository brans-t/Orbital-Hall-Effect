"""Dimensional restoration of arbitrary eV/Angstrom Hamiltonians."""

import numpy as np
import pytest
from ohk.units import UnitSystem, HBAR, M_E, ELECTRON_CHARGE, hbar, m_e, e


def test_internal_constants_are_unity():
    assert hbar == m_e == e == 1


def test_oam_restores_hbar_scaling():
    # Calibrate E0*a0^2 = hbar^2/m so one internal unit equals one hbar.
    energy = HBAR**2 / (M_E * ELECTRON_CHARGE * 1e-20)
    units = UnitSystem(energy=energy)
    assert units.oam_to_hbar(1.0) == pytest.approx(1)
    assert units.oam_to_si(1.0) == pytest.approx(HBAR, rel=1e-12, abs=0)
    # Arbitrary eV/Angstrom models need the additional physical scale.
    assert UnitSystem(energy=2, length=3).oam_scale == pytest.approx(18 * UnitSystem().oam_scale)


def test_conductivity_unit_conversion():
    u = UnitSystem()
    expected = u.oam_scale * ELECTRON_CHARGE / (2 * np.pi)
    assert u.conductivity_to_si(1.0) == pytest.approx(expected, rel=1e-12, abs=0)
    assert u.conductivity_to_si(1.0, dimension=3) == pytest.approx(
        expected / 1e-10, rel=1e-12, abs=0
    )


def test_invalid_unit_scale():
    with pytest.raises(ValueError):
        UnitSystem(length=0)
