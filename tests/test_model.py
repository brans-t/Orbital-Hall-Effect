"""Model construction, analytic derivatives and native Wannier90 I/O."""

import numpy as np
import pytest
from ohk.model import HaldaneModel, AnalyticModel, Hopping, ModelSpec
from ohk.model.wannier90 import Wannier90Model


def test_haldane_gamma_energies(haldane_topological):
    m = haldane_topological
    # d0 shifts BOTH energies; the old design incorrectly wrote a symmetric spectrum.
    expected = 6 * m.t2 * np.cos(m.phi) + np.array([-1, 1]) * np.sqrt(9 * m.t1**2 + m.M**2)
    np.testing.assert_allclose(m.energies([0, 0]), expected, atol=1e-13)


def test_haldane_traceless_would_be_wrong(haldane_topological):
    assert abs(np.trace(haldane_topological.H([0, 0]))) > 1


def test_haldane_hamiltonian_hermitian(haldane_topological, sampler):
    for k in sampler():
        m = haldane_topological
        m.validate(k)
        expected = np.array(
            [[m.d0(k) + m.dz(k), m.t1 * m.f(k)], [m.t1 * m.f(k).conjugate(), m.d0(k) - m.dz(k)]]
        )
        np.testing.assert_allclose(m.H(k), expected, atol=2e-15)


def test_haldane_analytic_derivative_matches_finite_difference(haldane_topological, sampler):
    for k in sampler():
        np.testing.assert_allclose(
            haldane_topological.dH(k),
            haldane_topological.dH_finite_difference(k),
            atol=2e-9,
            rtol=0,
        )


def test_lattice_vectors_give_correct_reciprocal(haldane_topological):
    m = haldane_topological
    np.testing.assert_allclose(
        m.lattice() @ m.spec.reciprocal_vectors().T, 2 * np.pi * np.eye(2), atol=2e-15
    )


def test_chern_number_phase_boundary(haldane_gap_closing_phi):
    for phi, expected in [
        (0.0, 0),
        (0.15 * np.pi, 0),
        (0.25 * np.pi, -1),
        (0.5 * np.pi, -1),
        (-0.25 * np.pi, 1),
    ]:
        m = HaldaneModel(phi=phi)
        assert m.chern_number(30) == expected
        assert m.chern_number(30, band=1) == -expected
    assert HaldaneModel().phase_boundary() == pytest.approx(haldane_gap_closing_phi)


def test_wannier90_round_trip(tmp_path):
    # Non-unit weights, complex hopping, skewed 3D lattice.
    a = np.array([[2.0, 0.2, 0.1], [0.0, 3.0, 0.4], [0.3, 0.0, 4.0]])
    rng = np.random.default_rng(4)
    t = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    on = np.array([[0.7, 0.3j], [-0.3j, -0.2]])
    spec = ModelSpec(2, a, periodic_dim=3)
    model = Wannier90Model(
        spec, np.array([on, 2 * t, 2 * t.conj().T]), np.array([[0, 0, 0], a[0], -a[0]]), [1, 2, 2]
    )
    path = tmp_path / "toy_hr.dat"
    model.to_hr_dat(path)
    loaded = Wannier90Model.from_hr_dat(path, a)
    for k in [[0.1, 0.2, -0.3], [0.7, -0.2, 0.4]]:
        expected = on + t * np.exp(1j * (a[0] @ k)) + t.conj().T * np.exp(-1j * (a[0] @ k))
        np.testing.assert_allclose(loaded.H(k), expected, atol=1e-14)
        np.testing.assert_allclose(loaded.dH(k), model.dH(k), atol=1e-14)
        np.testing.assert_allclose(loaded.dH(k), loaded.dH_finite_difference(k), atol=1e-8)
    np.testing.assert_array_equal(loaded.degeneracies, [1, 2, 2])


def test_wannier90_requires_lattice_vectors(tmp_path):
    with pytest.raises(ValueError, match="lattice_vectors"):
        Wannier90Model.from_hr_dat(tmp_path / "missing.dat")


@pytest.mark.parametrize("content", ["comment\n2\n1\n", "comment\n2\n1\n1\n0 0 0 1 1 0 0\n"])
def test_wannier90_rejects_truncated_file(tmp_path, content):
    p = tmp_path / "bad_hr.dat"
    p.write_text(content)
    with pytest.raises(ValueError):
        Wannier90Model.from_hr_dat(p, np.eye(2))


def test_hermitize_is_idempotent():
    spec = ModelSpec(2, np.eye(2))
    hop = Hopping(0, 1, (1, 0), 0.7 + 0.2j)
    m = AnalyticModel(spec, [hop])
    assert len(m.hoppings) == 2
    np.testing.assert_allclose(m.H([0.4, 0.3]), m.hermitize().H([0.4, 0.3]))
    complete = AnalyticModel(spec, m.hoppings)
    np.testing.assert_allclose(complete.H([0.4, 0.3]), m.H([0.4, 0.3]))
    with pytest.raises(ValueError, match="inconsistent"):
        AnalyticModel(spec, [hop, Hopping(1, 0, (-1, 0), 1.0)])


def test_empty_hopping_list():
    m = AnalyticModel(ModelSpec(2, np.eye(2)), [], onsite=np.array([1.0, 2.0]))
    np.testing.assert_allclose(m.H([1, 1]), np.diag([1, 2]))
    np.testing.assert_allclose(m.dH([1, 1]), 0.0)


def test_invalid_model_input():
    with pytest.raises(ValueError):
        ModelSpec(2, np.zeros((2, 2)))
    with pytest.raises(ValueError):
        Hopping(0, 1, (0.5, 0), 1)
    with pytest.raises(ValueError):
        HaldaneModel(nearest_neighbour_distance=0)
    with pytest.raises(ValueError):
        HaldaneModel().H([0, 0, 0])
