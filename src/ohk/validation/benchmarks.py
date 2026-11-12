"""Executable Haldane benchmarks. Supersedes the inconsistent design-phase numbers.

The sign=-1 spectral expression matches the PUBLISHED denominator convention. A sign=+1
control does not change any two-band result. Numerical OHC includes the e/(2*pi) conversion.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class HaldaneBenchmark:
    t1: float = 1.0
    t2: float = 0.15
    M: float = 0.4
    phi_c_over_pi: float = float(np.arcsin(0.4 / (3 * np.sqrt(3) * 0.15)) / np.pi)
    chern_by_phi: dict = field(default_factory=lambda: {0.0: 0, 0.15: 0, 0.25: -1, 0.5: -1})
    min_bands_for_sum_term: int = 3
    sigma_oh_polarity_corrected_sign: int = -1


HALDANE_BENCHMARK = HaldaneBenchmark()


def run_haldane_benchmark(variant=None, n_k=60):
    from ohk.model.haldane import HaldaneModel
    from ohk.oam.covariant import CovariantProjector
    from ohk.oam.spectral import SpectralSum
    from ohk.oam.assemble import full_matrix
    from ohk.oam.pezo import PezoTruncated
    from ohk.transport.grids import BZGrid
    from ohk.transport.kubo import orbital_hall_conductivity
    from ohk.eigen import EigenBasis
    from ohk.validation.checks import check_gauge_invariance, check_time_reversal
    from ohk.validation.oracle import make_random_multiband_model, compare_bloch_vs_lattice
    from ohk.model.finite import oam_from_velocity

    variant = CovariantProjector() if variant is None else variant
    m = HaldaneModel(phi=np.pi / 4)
    checks = {}
    k = np.array([0.71, 1.13])
    for n in (3, 4):
        random = make_random_multiband_model(n_bands=n, seed=n)
        basis = EigenBasis.at(random, k)
        direct = oam_from_velocity(basis.eps, basis.v)
        checks[f"projector_vs_definition_{n}band"] = float(
            np.max(abs(full_matrix(CovariantProjector(), basis)[0] - direct))
        )
        checks[f"spectral_vs_definition_{n}band"] = float(
            np.max(abs(SpectralSum().matrix(basis) - direct))
        )
        checks[f"gauge_covariance_{n}band"] = check_gauge_invariance(random, variant, k)
    basis = EigenBasis.at(m, k)
    checks["projector_vs_spectral"] = float(
        np.max(abs(full_matrix(CovariantProjector(), basis)[0] - SpectralSum().matrix(basis)))
    )
    checks["two_band_sign_control_equivalence"] = float(
        np.max(abs(SpectralSum(1).matrix(basis) - SpectralSum(-1).matrix(basis)))
    )
    checks["pezo_two_band_off_diagonal"] = float(abs(PezoTruncated().matrix(basis)[0, 1]))
    checks["time_reversal"] = check_time_reversal(HaldaneModel(phi=0), variant, k)
    checks["finite_lattice"] = compare_bloch_vs_lattice(m, CovariantProjector(), 4)["max_abs_diff"]
    chern = {
        phi: HaldaneModel(phi=phi * np.pi).chern_number(max(24, n_k))
        for phi in HALDANE_BENCHMARK.chern_by_phi
    }
    checks["chern_phase_values"] = float(
        max(abs(chern[p] - c) for p, c in HALDANE_BENCHMARK.chern_by_phi.items())
    )
    sigma = orbital_hall_conductivity(m, variant, BZGrid.uniform(m.lattice(), n_k)).sigma
    # Polarity is a benchmark only for the full OAM definitions.
    if variant.name in ("covariant-16a", "spectral-16b"):
        checks["ohc_polarity"] = 0.0 if sigma < 0 else 1.0
    return {
        "sigma_oh": sigma,
        "sigma_units": "e/(2*pi)",
        "chern_by_phi": chern,
        "gamma_energies": m.energies([0, 0]),
        "phi_c_over_pi": m.phase_boundary() / np.pi,
        "checks": checks,
        "passed": all(np.isfinite(x) and x < 1e-8 for x in checks.values()),
    }
