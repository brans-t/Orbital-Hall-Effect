"""Large-sweep acceleration must agree with the established public calculations."""

import numpy as np
import pytest

from ohk import EigenBasis, CovariantProjector, full_matrix, orbital_hall_conductivity, BZGrid
from ohk.model.finite import FiniteLattice
from ohk.validation.reproduction import bloch_response, finite_pbc_response


@pytest.mark.parametrize(
    "fixture", ["haldane_topological", "random_three_band", "random_four_band"]
)
def test_batched_projector_matches_public_api(fixture, request):
    model = request.getfixturevalue(fixture)
    grid = BZGrid.uniform(model.lattice(), (5, 4))
    response = bloch_response(model, grid.kpoints)
    for index, k in enumerate(grid.kpoints):
        basis = EigenBasis.at(model, k)
        np.testing.assert_allclose(response.energies[index], basis.eps, atol=1e-12)
        np.testing.assert_allclose(
            response.oam[index], full_matrix(CovariantProjector(), basis)[0], atol=1e-12
        )
    expected = orbital_hall_conductivity(model, CovariantProjector(), grid).sigma
    assert response.sigma(mu=0) == pytest.approx(expected, abs=1e-12)
    np.testing.assert_allclose(
        response.fermi_scan([-20, -0.2, 0, 0.2, 20]),
        [response.sigma(mu=mu) for mu in [-20, -0.2, 0, 0.2, 20]],
        atol=1e-12,
    )


def test_sparse_torus_matches_dense_oracle(haldane_topological):
    model = haldane_topological
    small = FiniteLattice(model, 3).compute()
    response = finite_pbc_response(model, 3)
    np.testing.assert_allclose(np.sort(response.energies.ravel()), small.energies, atol=1e-12)
    bloch = bloch_response(model, response.kpoints)
    np.testing.assert_allclose(
        response.oam.diagonal(axis1=-2, axis2=-1),
        bloch.oam.diagonal(axis1=-2, axis2=-1),
        atol=1e-12,
    )
    np.testing.assert_allclose(abs(response.oam), abs(bloch.oam), atol=1e-12)
    np.testing.assert_allclose(response.kernels, bloch.kernels, atol=1e-11)


def test_sparse_torus_is_independent_of_bloch_model_methods(haldane_topological, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Finite PBC must use real-space operators")

    monkeypatch.setattr(haldane_topological, "H", forbidden)
    monkeypatch.setattr(haldane_topological, "dH", forbidden)
    monkeypatch.setattr(CovariantProjector, "matrix", forbidden)
    response = finite_pbc_response(haldane_topological, 6)
    assert response.sigma(mu=0) < 0
