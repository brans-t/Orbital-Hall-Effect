"""Finite PBC spectra and OAM, preserving hopping winding and momentum sectors."""

import numpy as np
import pytest
from ohk.model.finite import FiniteLattice
from ohk.oam import CovariantProjector
from ohk.validation.oracle import compare_bloch_vs_lattice


@pytest.mark.slow
def test_oracle_matches_bloch_haldane(haldane_topological):
    result = compare_bloch_vs_lattice(haldane_topological, CovariantProjector(), 12)
    assert result["max_abs_diff"] < 1e-10


def test_oracle_matches_multiband_and_rectangular_cell(random_three_band):
    lattice = FiniteLattice(random_three_band, (3, 4))
    result = lattice.compare_with_bloch(random_three_band, [(0, 0), (0, 1), (1, 2)])
    assert result["max_abs_diff"] < 1e-10


def test_oracle_rejects_diagonal_position_elements():
    with pytest.raises(ValueError, match="energy-basis"):
        FiniteLattice.oam_operator(np.diag([0.0, 1.0]), np.eye(2), drop_diagonal=False)


def test_oracle_hamiltonian_is_hermitian(haldane_topological):
    lattice = FiniteLattice(haldane_topological, 4)
    h = lattice.build_hamiltonian()
    np.testing.assert_allclose(h, h.conj().T, atol=1e-15)
    _, velocity = lattice._build_operators()
    np.testing.assert_allclose(velocity, velocity.conj().swapaxes(-1, -2), atol=1e-15)


def test_supercell_energies_match_bloch(haldane_topological):
    lattice = FiniteLattice(haldane_topological, 4)
    expected = sorted(e for k in lattice.allowed_momenta() for e in haldane_topological.energies(k))
    actual = np.linalg.eigvalsh(lattice.build_hamiltonian())
    np.testing.assert_allclose(actual, expected, atol=1e-12)
    result = lattice.compute()
    np.testing.assert_allclose(result.energies, expected, atol=1e-12)
    h = lattice.build_hamiltonian()
    np.testing.assert_allclose(
        h @ result.eigenvectors, result.eigenvectors * result.energies, atol=1e-12
    )


def test_oracle_only_accepts_allowed_momenta(haldane_topological):
    with pytest.raises(ValueError, match="not allowed"):
        compare_bloch_vs_lattice(
            haldane_topological, CovariantProjector(), 3, momenta=np.array([[0.12, 0.7]])
        )


def test_finite_reference_does_not_call_bloch_formulas(haldane_topological, monkeypatch):
    lattice = FiniteLattice(haldane_topological, 3)

    def forbidden(*args, **kwargs):
        raise AssertionError("The finite reference must not call a Bloch formula")

    monkeypatch.setattr(haldane_topological, "H", forbidden)
    monkeypatch.setattr(haldane_topological, "dH", forbidden)
    monkeypatch.setattr(CovariantProjector, "matrix", forbidden)
    result = lattice.compute()
    assert np.isfinite(result.L).all()
    assert np.max(abs(result.L)) > 0.01


def test_three_dimensional_model_matches_real_space():
    from ohk.model import AnalyticModel, Hopping, ModelSpec

    model = AnalyticModel(
        ModelSpec(3, np.diag([1.2, 1.3, 1.4]), periodic_dim=3),
        [
            Hopping(0, 1, (1, 0, 1), 0.2 + 0.1j),
            Hopping(1, 2, (0, 1, 0), -0.3j),
            Hopping(0, 2, (1, 1, 0), 0.1),
        ],
        onsite=np.array([-1.0, 0.3, 1.1]),
    )
    lattice = FiniteLattice(model, (2, 3, 2))
    assert lattice.compare_with_bloch(model, [(0, 0), (0, 1), (1, 2)])["passed"]
