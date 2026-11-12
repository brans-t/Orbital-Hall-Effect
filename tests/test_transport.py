"""Kubo normalization, occupation, symmetry, broadening and convergence."""

import numpy as np
import pytest
from ohk.eigen import EigenBasis, DegeneracyError
from ohk.model import HaldaneModel
from ohk.oam import CovariantProjector, SpectralSum
from ohk.oam.assemble import full_matrix
from ohk.transport.current import orbital_current, orbital_current_tensor
from ohk.transport.grids import BZGrid
from ohk.transport.kubo import (
    orbital_hall_conductivity,
    orbital_hall_symbol,
    fermi_scan,
    convergence_study,
)


def test_ohc_sign_matches_paper(haldane_topological):
    m = haldane_topological
    r = orbital_hall_conductivity(m, CovariantProjector(), BZGrid.uniform(m.lattice(), 30))
    assert r.sigma < 0
    assert r.degenerate_fraction == 0
    assert r.min_gap > 0.3


def test_two_band_sum_sign_cannot_flip_ohc(haldane_topological):
    m = haldane_topological
    grid = BZGrid.uniform(m.lattice(), 18)
    a = orbital_hall_conductivity(m, SpectralSum(-1), grid)
    b = orbital_hall_conductivity(m, SpectralSum(1), grid)
    assert a.sigma == pytest.approx(b.sigma, abs=1e-13)


def test_ohc_sixty_mesh_reference(haldane_topological):
    m = haldane_topological
    r = orbital_hall_conductivity(m, CovariantProjector(), BZGrid.uniform(m.lattice(), 60))
    # Reproducible mesh value in e/(2*pi); raw e units would be -0.1689667...
    assert r.sigma == pytest.approx(-1.0616491612557013, abs=1e-10)


def test_ohc_grows_near_phase_boundary(haldane_gap_closing_phi):
    values = []
    for offset in [0.1, 0.01]:
        m = HaldaneModel(phi=haldane_gap_closing_phi + offset)
        # Follow the gap energy, which is NOT zero because d0 shifts the Dirac cones.
        mu = -3 * m.t2 * np.cos(m.phi)
        values.append(
            abs(
                orbital_hall_conductivity(
                    m, CovariantProjector(), BZGrid.uniform(m.lattice(), 30), mu=mu
                ).sigma
            )
        )
    assert values[1] > 10 * values[0]


def test_fermi_scan_plateau():
    m = HaldaneModel(phi=np.pi / 2)
    scan = fermi_scan(
        m, CovariantProjector(), BZGrid.uniform(m.lattice(), 24), [-20, -0.1, 0.0, 0.1, 20]
    )
    assert abs(scan[0]) < 1e-12 and abs(scan[-1]) < 1e-12
    np.testing.assert_allclose(scan[1:4], scan[2], atol=1e-12)
    assert scan[2] < -0.1


def test_orbital_current_is_hermitian(random_three_band):
    b = EigenBasis.at(random_three_band, [0.3, 0.7])
    L = full_matrix(CovariantProjector(), b)[0]
    J = orbital_current_tensor(L, b.v)
    np.testing.assert_allclose(J, J.conj().swapaxes(-1, -2), atol=1e-12)
    assert np.linalg.norm(b.v[0] @ L - L @ b.v[0]) > 0.001


def test_kubo_matches_literal_double_sum(random_three_band):
    m = random_three_band
    grid = BZGrid.uniform(m.lattice(), (5, 4), shift=0.5)
    total = 0.0
    for k, weight in zip(grid.kpoints, grid.weights):
        b = EigenBasis.at(m, k)
        L = full_matrix(CovariantProjector(), b)[0]
        J = orbital_current(L, b.v)
        for n in range(m.num_bands):
            for p in range(m.num_bands):
                if n != p and b.eps[n] < 0:
                    total += (
                        weight * 2 * np.imag(J[n, p] * b.v[1, p, n]) / (b.eps[n] - b.eps[p]) ** 2
                    )
    expected = total * 2 * np.pi / grid.cell_area()
    result = orbital_hall_conductivity(m, CovariantProjector(), grid)
    assert result.sigma == pytest.approx(expected, abs=1e-12)


def test_fermi_scan_reuses_diagonalization(haldane_topological, monkeypatch):
    m = haldane_topological
    grid = BZGrid.uniform(m.lattice(), 5)
    original, calls = np.linalg.eigh, []

    def counted(h):
        calls.append(1)
        return original(h)

    monkeypatch.setattr(np.linalg, "eigh", counted)
    mus = [-0.5, 0.0, 0.5]
    scan = fermi_scan(m, CovariantProjector(), grid, mus)
    assert len(calls) == len(grid)
    expected = [orbital_hall_conductivity(m, CovariantProjector(), grid, mu=x).sigma for x in mus]
    np.testing.assert_allclose(scan, expected, atol=1e-12)


def test_unnormalized_grid_is_equivalent(haldane_topological):
    m = haldane_topological
    grid = BZGrid.uniform(m.lattice(), 8)
    raw = BZGrid(grid.kpoints, grid.weights * grid.reciprocal_volume(), False, m.lattice())
    a = orbital_hall_conductivity(m, CovariantProjector(), grid).sigma
    b = orbital_hall_conductivity(m, CovariantProjector(), raw).sigma
    assert a == pytest.approx(b, abs=1e-12)


def test_exact_degeneracy_requires_explicit_regularization(haldane_gap_closing_phi):
    m = HaldaneModel(phi=haldane_gap_closing_phi)
    grid = BZGrid.uniform(m.lattice(), 3)
    with pytest.raises(DegeneracyError):
        orbital_hall_conductivity(m, CovariantProjector(), grid)
    with pytest.warns(UserWarning, match="near-degenerate"):
        r = orbital_hall_conductivity(
            m, CovariantProjector(), grid, subband_regularization=True, eta=0.001
        )
    assert np.isfinite(r.sigma)
    assert r.degenerate_fraction > 0 and r.regularization_eta == 0.001


def test_axis_swap(haldane_topological):
    m = haldane_topological
    grid = BZGrid.uniform(m.lattice(), 24)
    a = orbital_hall_conductivity(m, CovariantProjector(), grid).sigma
    b = orbital_hall_symbol(m, CovariantProjector(), grid).sigma
    assert a == pytest.approx(-b, abs=1e-12)


def test_convergence_study_returns_measured_changes(haldane_topological):
    r = convergence_study(haldane_topological, CovariantProjector(), [6, 12, 18])
    assert len(r["sigma"]) == 3
    np.testing.assert_allclose(r["successive_changes"], np.diff(r["sigma"]))
    assert r["extrapolated"] is None


def test_path_is_not_an_integration_grid(haldane_topological):
    m = haldane_topological
    grid = BZGrid.path(m.lattice(), [[0, 0], [0.5, 0.5]], ["G", "M"], 3)
    with pytest.raises(ValueError, match="not a BZ"):
        orbital_hall_conductivity(m, CovariantProjector(), grid)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"mu": float("nan")},
        {"eta": 0.1},
        {"subband_regularization": True, "eta": -1},
        {"current_direction": 3},
        {"degeneracy_tol": 0},
    ],
)
def test_invalid_transport_inputs(haldane_topological, kwargs):
    m = haldane_topological
    with pytest.raises(ValueError):
        orbital_hall_conductivity(m, CovariantProjector(), BZGrid.uniform(m.lattice(), 3), **kwargs)


@pytest.mark.parametrize("n", [0, 2.5, [2, 0], [2, 3, 4]])
def test_invalid_grid_sizes(haldane_topological, n):
    with pytest.raises(ValueError):
        BZGrid.uniform(haldane_topological.lattice(), n)


def test_invalid_weights():
    with pytest.raises(ValueError, match="sum"):
        BZGrid(np.zeros((2, 2)), np.ones(2), True, np.eye(2))
