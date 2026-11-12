"""Physical symmetry checks with band-phase alignment."""

import numpy as np
from ohk.model import HaldaneModel
from ohk.eigen import EigenBasis
from ohk.oam import CovariantProjector
from ohk.oam.assemble import full_matrix
from ohk.validation.checks import check_time_reversal, check_gauge_invariance, check_inversion


def test_haldane_phi_zero_has_time_reversal(haldane_trivial, sampler):
    for k in sampler():
        assert check_time_reversal(haldane_trivial, CovariantProjector(), k) < 1e-12


def test_haldane_phi_quarter_breaks_time_reversal(haldane_topological, sampler):
    assert (
        max(check_time_reversal(haldane_topological, CovariantProjector(), k) for k in sampler())
        > 0.01
    )


def test_diagonal_odd_under_k_flip_in_trivial_phase(haldane_trivial):
    k = np.array([0.7, 1.1])
    a = full_matrix(CovariantProjector(), EigenBasis.at(haldane_trivial, k))[0]
    b = full_matrix(CovariantProjector(), EigenBasis.at(haldane_trivial, -k))[0]
    np.testing.assert_allclose(a.diagonal(), -b.diagonal(), atol=1e-12)


def test_l12_vanishes_at_phi_half_pi(sampler):
    # d0=0 at pi/2, but the model is NOT degenerate there for the default mass.
    model = HaldaneModel(phi=np.pi / 2)
    for k in sampler():
        b = EigenBasis.at(model, k)
        assert abs(full_matrix(CovariantProjector(), b)[0][0, 1]) < 1e-12


def test_l11_equals_l22_for_two_band_model(haldane_topological, sampler):
    for k in sampler():
        L = full_matrix(CovariantProjector(), EigenBasis.at(haldane_topological, k))[0]
        np.testing.assert_allclose(L[0, 0], L[1, 1], atol=1e-12)


def test_hermiticity(random_three_band, sampler):
    for k in sampler():
        L = full_matrix(CovariantProjector(), EigenBasis.at(random_three_band, k))[0]
        np.testing.assert_allclose(L, L.conj().T, atol=1e-12)


def test_gauge_covariance(random_four_band, sampler):
    for k in sampler():
        assert check_gauge_invariance(random_four_band, CovariantProjector(), k) < 1e-12


def test_inversion_with_explicit_orbital_operator(sampler):
    model = HaldaneModel(M=0, phi=np.pi / 4)
    model.spec.metadata["inversion_operator"] = np.array([[0, 1], [1, 0]])
    for k in sampler():
        assert check_inversion(model, CovariantProjector(), k) < 1e-12
