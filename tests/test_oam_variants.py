"""OAM approximation semantics and assembler failure detection."""

import numpy as np
import pytest
from ohk.eigen import EigenBasis
from ohk.oam import CovariantProjector, PezoTruncated, AtomCentered, OAMVariant
from ohk.oam.atomic import ACASpec
from ohk.oam.assemble import full_matrix, compare


def test_pezo_misses_two_band_off_diagonal(haldane_topological):
    b = EigenBasis.at(haldane_topological, [0.7, 1.1])
    full = full_matrix(CovariantProjector(), b)[0]
    trunc = full_matrix(PezoTruncated(), b)[0]
    assert abs(full[0, 1]) > 0.01
    assert trunc[0, 1] == 0
    np.testing.assert_allclose(trunc.diagonal(), full.diagonal(), atol=1e-12)


def test_aca_s_orbitals_are_zero(haldane_topological):
    b = EigenBasis.at(haldane_topological, [0.7, 1.1])
    L, report = full_matrix(AtomCentered(ACASpec((0, 0))), b)
    np.testing.assert_array_equal(L, 0)
    assert not report.diagonal_filled


def test_atomic_multiplet_algebra():
    from ohk.model import AnalyticModel, ModelSpec

    m = AnalyticModel(ModelSpec(3, np.eye(2)), [], onsite=np.array([-1, 0, 1]))
    b = EigenBasis.at(m, [0, 0])
    x, y, z = [AtomCentered(ACASpec((1,), axis)).matrix(b) for axis in ("x", "y", "z")]
    np.testing.assert_allclose(x @ y - y @ x, 1j * z, atol=1e-14)
    np.testing.assert_allclose(x @ x + y @ y + z @ z, 2 * np.eye(3), atol=1e-14)


def test_aca_requires_complete_multiplets(haldane_topological):
    b = EigenBasis.at(haldane_topological, [0, 0])
    with pytest.raises(ValueError, match="exceed"):
        AtomCentered(ACASpec((1,))).matrix(b)


def test_assembler_rejects_nonhermitian_variant(haldane_topological):
    class Bad(OAMVariant):
        name = "bad"

        def matrix(self, basis):
            return np.array([[0, 1], [0, 0]], complex)

    with pytest.raises(ValueError, match="Non-Hermitian"):
        full_matrix(Bad(), EigenBasis.at(haldane_topological, [0.3, 0.7]))


def test_variant_flag_cannot_be_shadowed():
    with pytest.raises(TypeError, match="must be a bool"):

        class Bad(OAMVariant):
            def diagonal(self):
                return 0

            def matrix(self, basis):
                return 0


def test_comparison_shares_basis(haldane_topological):
    b = EigenBasis.at(haldane_topological, [0.3, 0.7])
    result = compare([CovariantProjector(), PezoTruncated()], b)
    assert len(result) == 2
