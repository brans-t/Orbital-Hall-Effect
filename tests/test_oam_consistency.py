"""Independent elementwise definition versus projector and spectral sum; >=3 bands mandatory."""

import numpy as np
import pytest
from ohk.eigen import EigenBasis
from ohk.oam import CovariantProjector, SpectralSum
from ohk.oam.assemble import full_matrix, diagonal_closed_form

TOL = 1e-12


def direct_definition(basis, model=None):
    n = basis.n_bands
    r = np.zeros_like(basis.v)
    for i in range(n):
        for j in range(n):
            if i != j:
                r[:, i, j] = -1j * basis.v[:, i, j] / (basis.eps[i] - basis.eps[j])
    result = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(n):
            for m in range(n):
                result[i, j] += 0.5 * (
                    r[0, i, m] * basis.v[1, m, j]
                    - r[1, i, m] * basis.v[0, m, j]
                    - basis.v[0, i, m] * r[1, m, j]
                    + basis.v[1, i, m] * r[0, m, j]
                )
    return result


@pytest.mark.parametrize(
    "model_fixture", ["haldane_topological", "random_three_band", "random_four_band"]
)
def test_16a_matches_direct_definition(model_fixture, request, sampler):
    model = request.getfixturevalue(model_fixture)
    for k in sampler(7):
        b = EigenBasis.at(model, k)
        np.testing.assert_allclose(
            full_matrix(CovariantProjector(), b)[0], direct_definition(b), atol=TOL, rtol=0
        )
        np.testing.assert_allclose(
            CovariantProjector(True).matrix(b), direct_definition(b), atol=TOL, rtol=0
        )


@pytest.mark.parametrize(
    "model_fixture", ["haldane_topological", "random_three_band", "random_four_band"]
)
def test_16b_matches_direct_definition(model_fixture, request, sampler):
    model = request.getfixturevalue(model_fixture)
    for k in sampler(7):
        b = EigenBasis.at(model, k)
        np.testing.assert_allclose(SpectralSum().matrix(b), direct_definition(b), atol=TOL, rtol=0)


def test_wrong_transcription_sign_is_detected_above_two_bands(random_three_band, sampler):
    residuals = []
    for k in sampler():
        b = EigenBasis.at(random_three_band, k)
        residuals.append(np.max(abs(SpectralSum(1).matrix(b) - direct_definition(b))))
    assert max(residuals) > 0.01


def test_sum_sign_is_invisible_for_two_bands(haldane_topological, sampler):
    for k in sampler():
        b = EigenBasis.at(haldane_topological, k)
        np.testing.assert_allclose(SpectralSum(1).matrix(b), SpectralSum(-1).matrix(b), atol=TOL)


def test_diagonal_closed_form_matches_definition(random_three_band, sampler):
    for k in sampler():
        b = EigenBasis.at(random_three_band, k)
        np.testing.assert_allclose(
            diagonal_closed_form(b), direct_definition(b).diagonal(), atol=TOL, rtol=0
        )


def test_diagonal_is_real(random_three_band, sampler):
    for k in sampler():
        L = SpectralSum().matrix(EigenBasis.at(random_three_band, k))
        assert np.max(abs(L.diagonal().imag)) < TOL


def test_published_denominator_convention(random_four_band):
    b = EigenBasis.at(random_four_band, [0.4, -0.7])
    expected = np.diag(diagonal_closed_form(b)).astype(complex)
    # Paper: PLUS sum, denominators eps_m-eps_n.
    for n in range(b.n_bands):
        for p in range(b.n_bands):
            if n == p:
                continue
            vx, vy = b.v[:2]
            term = ((vx[n, n] + vx[p, p]) * vy[n, p] - (vy[n, n] + vy[p, p]) * vx[n, p]) / (
                b.eps[n] - b.eps[p]
            )
            for m in range(b.n_bands):
                if m != n and m != p:
                    term += (1 / (b.eps[m] - b.eps[n]) + 1 / (b.eps[m] - b.eps[p])) * (
                        vx[n, m] * vy[m, p] - vy[n, m] * vx[m, p]
                    )
            expected[n, p] = 0.5j * term
    np.testing.assert_allclose(expected, SpectralSum().matrix(b), atol=TOL)
