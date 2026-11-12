"""Shared test fixtures.

Design note
-----------
Every fixture here is deliberately small and fast enough to run on each commit.  The expensive
cross-checks (finite-lattice oracle at ``N=24``, high ``k``-count conductivity sweeps) are marked
``slow`` and included by default; use ``-m "not slow"`` to omit them.
"""

from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture(scope="session")
def haldane_t2():
    return 0.15


@pytest.fixture(scope="session")
def haldane_m():
    return 0.4


@pytest.fixture(scope="session")
def haldane_gap_closing_phi():
    """``phi_c = arcsin(M / (3 sqrt(3) t2))`` in radians, the topological transition."""
    return np.arcsin(0.4 / (3.0 * np.sqrt(3.0) * 0.15))


@pytest.fixture
def haldane_trivial():
    """Haldane model with time-reversal symmetry (``phi = 0``)."""
    from ohk.model.haldane import HaldaneModel

    return HaldaneModel(t1=1.0, t2=0.15, M=0.4, phi=0.0)


@pytest.fixture
def haldane_topological():
    """Haldane model in the Chern-insulating phase (``phi = pi/4``)."""
    from ohk.model.haldane import HaldaneModel

    return HaldaneModel(t1=1.0, t2=0.15, M=0.4, phi=0.25 * np.pi)


@pytest.fixture
def random_three_band():
    """Random 3-band model — the minimum needed to exercise the Eq. (16b) sum term.

    Do not replace this with a two-band model: the summation over ``m != n, n'`` is empty for two
    bands, which makes the printed-sign error in Eq. (16b) undetectable
    (docs/design/04-validation-report.md §3.2).
    """
    from ohk.validation.oracle import make_random_multiband_model

    return make_random_multiband_model(n_bands=3, seed=7)


@pytest.fixture
def random_four_band():
    from ohk.validation.oracle import make_random_multiband_model

    return make_random_multiband_model(n_bands=4, seed=11)


@pytest.fixture
def sampler():
    """Deterministic momentum sampler over a few fractional coordinates."""

    def _sample(n: int = 5, seed: int = 0) -> np.ndarray:
        rng = np.random.default_rng(seed)
        return rng.uniform(-np.pi, np.pi, size=(n, 2))

    return _sample
