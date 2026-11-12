"""Eigenbasis contracts and independent derivatives."""

import numpy as np
import pytest
from ohk.eigen import EigenBasis, DegeneracyError
from ohk.model import HaldaneModel, AnalyticModel, Hopping, ModelSpec, TightBindingModel


def test_velocity_diagonal_is_real(haldane_topological, sampler):
    for k in sampler():
        b = EigenBasis.at(haldane_topological, k)
        assert np.max(abs(np.diagonal(b.v, axis1=-2, axis2=-1).imag)) < 1e-14


def test_velocity_diagonal_equals_band_gradient(haldane_topological, sampler):
    m = haldane_topological
    for k in sampler():
        b = EigenBasis.at(m, k)
        for j, delta in enumerate(np.eye(2) * 1e-5):
            gradient = (m.energies(k + delta) - m.energies(k - delta)) / 2e-5
            np.testing.assert_allclose(b.v[j].diagonal(), gradient, atol=1e-8)


def test_eigenbasis_shape_contract(haldane_topological):
    b = EigenBasis.at(haldane_topological, [0.0, 0.0])
    assert b.eps.shape == (2,)
    assert b.v.shape == b.C.shape == b.dH.shape == (2, 2, 2)
    assert b.C is b.derivative_coefficients()


def test_derivative_coefficients_have_zero_diagonal(random_three_band):
    b = EigenBasis.at(random_three_band, [0.3, -0.7])
    np.testing.assert_array_equal(np.diagonal(b.C, axis1=-2, axis2=-1), 0.0)


def test_derivative_coefficients_reproduce_spectral_expansion(random_three_band):
    m, k, step = random_three_band, np.array([0.3, -0.7]), 1e-6
    b = EigenBasis.at(m, k)
    for j, offset in enumerate(np.eye(2) * step):
        vectors = []
        for point in (k + offset, k - offset):
            _, u = np.linalg.eigh(m.H(point))
            overlap = np.diag(b.U.conj().T @ u)
            vectors.append(u * (overlap.conj() / abs(overlap)))
        du = (vectors[0] - vectors[1]) / (2 * step)
        np.testing.assert_allclose(b.U.conj().T @ du, b.C[j], atol=3e-9, rtol=0)


def test_non_hermitian_hamiltonian_is_rejected():
    m = AnalyticModel(ModelSpec(2, np.eye(2)), [Hopping(0, 1, (0, 0), 1)], assume_hermitian=True)
    with pytest.raises(ValueError, match="Hermitian"):
        EigenBasis.at(m, [0, 0])


def test_min_gap_and_degeneracy_mask(haldane_gap_closing_phi):
    m = HaldaneModel(phi=haldane_gap_closing_phi)
    candidates = [
        np.array(u) @ m.spec.reciprocal_vectors() for u in [(1 / 3, 2 / 3), (2 / 3, 1 / 3)]
    ]
    bases = [EigenBasis.at(m, k) for k in candidates]
    b = min(bases, key=lambda b: b.min_gap())
    assert b.min_gap() < 1e-12
    assert b.degenerate_mask().all()
    with pytest.raises(DegeneracyError):
        b.derivative_coefficients()


def test_eigenbasis_check_passes(random_four_band, sampler):
    for k in sampler():
        assert max(EigenBasis.at(random_four_band, k).check().values()) < 1e-12


def test_berry_connection_needs_gauge_derivatives(haldane_topological):
    b = EigenBasis.at(haldane_topological, [0.3, 0.7])
    with pytest.raises(ValueError, match="explicit"):
        b.berry_connection()
    du = b.U @ b.C
    np.testing.assert_allclose(b.berry_connection(du), 0, atol=1e-15)


def test_black_box_fallback_warns_once():
    class BlackBox(TightBindingModel):
        spec = ModelSpec(2, np.eye(2))

        def H(self, k):
            return np.diag([np.cos(k[0]), 2 + np.sin(k[1])])

    m = BlackBox()
    with pytest.warns(UserWarning, match="finite-difference"):
        np.testing.assert_allclose(
            m.dH([0.1, 0.2])[:, [0, 1], [0, 1]], [[-np.sin(0.1), 0], [0, np.cos(0.2)]], atol=1e-9
        )
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        m.dH([0.1, 0.2])
    assert not caught
