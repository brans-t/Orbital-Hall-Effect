# Orbital Hall kit (ohk)

A Python scientific computing library following the original L0-L5 architecture. It computes orbital angular momentum matrices and orbital Hall conductivity from tight-binding or Wannier90 Hamiltonians, with an independent finite-periodic-lattice reference calculation.

Based on Gautam & Satpathy, [Phys. Rev. B **114**, 214403 (2026)](https://doi.org/10.1103/pd3s-q2z7). Corrections to the original design's denominator transcription, phase conventions, periodic boundaries, and units are documented in the [implementation and validation notes](docs/design/05-implementation-notes.md).

## Implementation status

| Stage | Implemented |
|---|---|
| P0 | General hopping models, complete Haldane model, analytic derivatives, eigenbasis, FHS Chern number |
| P1 | Projected OAM expression and independent finite-periodic-lattice reference |
| P2 | Spectral sum, Pezo truncation, ACA, and three-/four-band cross-checks |
| P3 | Kubo/BZ integration, Fermi-energy scans, mesh convergence sequences, explicit broadening |
| P4 | JSON reports, CSV path scans, and three command-line entry points |
| P5 | Native Wannier90 I/O, degeneracy weights, 3D models, and an end-to-end smoke test with real Fe data |

The real-material test validates file loading and the calculation pipeline. Material predictions still require checks of mesh convergence, broadening, and the position-operator approximation.

## Installation and validation

Python 3.10 or later:

~~~bash
python -m pip install -e ".[dev]"
python -m pytest
ohk-validate --n-k 30 --output validation.json
~~~

The current Windows workspace also has a project environment ready to use:

~~~powershell
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\ohk-validate.exe --n-k 30
~~~

## Python API

~~~python
import numpy as np
from ohk import (
    HaldaneModel, EigenBasis, CovariantProjector,
    BZGrid, full_matrix, orbital_hall_conductivity, fermi_scan,
)

model = HaldaneModel(t1=1.0, t2=0.15, M=0.4, phi=np.pi / 4)
variant = CovariantProjector()
basis = EigenBasis.at(model, [0.7, 1.1])
L, diagnostics = full_matrix(variant, basis)

grid = BZGrid.uniform(model.lattice(), n=60)
result = orbital_hall_conductivity(model, variant, grid, mu=0.0)
print(result.sigma, result.units)
# -1.0616491612557013, e/(2*pi)

scan = fermi_scan(model, variant, grid, [-0.5, 0.0, 0.5])
~~~

By default, `variant.matrix(basis)` supplies only the off-diagonal part for the projected and truncated variants. Use `full_matrix` to fill the diagonal and check Hermiticity, or use `CovariantProjector(include_diagonal=True)`.

`k` is Cartesian, lattice vectors are rows, and the first index of `v[j,n,m]` is the direction. The Haldane model includes intracell orbital positions. Positions in `AnalyticModel` default to zero and can be supplied through `orbital_positions`.

## Command line

~~~bash
ohc-run --phi 0.7853981633974483 --n-k 60 --mu 0 --output ohc.json
ohc-run --n-k 30 --fermi -0.5 0 0.5 --convergence 20 40 60
oam-scan --points-per-segment 80 --output oam.csv
ohk-validate --n-k 30 --output validation.json
python examples/haldane.py --n-k 30 --output output/haldane.json
~~~

`--phi` is in radians. `--variant` accepts `covariant`, `spectral`, `pezo`, or `aca`. The default scan follows Γ-K-K′. Use `--vertices "0,0;0.333333,0.666667"` to specify fractional reciprocal coordinates. The CSV contains band energies and the real and imaginary parts of every OAM matrix element.

## Wannier90 and the material example

An `_hr.dat` file does not contain the lattice. Supply a separate text matrix with lattice vectors as rows:

~~~bash
ohc-run --hr material_hr.dat --lattice lattice.txt --n-k 12 --mu 9.1839 --output material.json
~~~

The public WannierTools Fe example can be downloaded at a pinned revision. The script verifies SHA-256 checksums for the archive and data file and stores them only in the unversioned `.local/fe/` directory:

~~~bash
python scripts/fetch_fe_example.py
ohc-run --hr .local/fe/wannier90_hr.dat --lattice .local/fe/lattice.txt --n-k 3 --mu 9.1839 --eta 0.01 --output .local/fe/smoke.json
~~~

This command is a coarse-mesh smoke test. `--eta` explicitly enables regularization. The Hamiltonian-only interface uses the periodic cell gauge and a diagonal-position approximation. It does not import `_r.dat`, the full anomalous-position matrix, or WS-distance corrections, and should not be equated with a complete first-principles OHC calculation. The Python API also provides a limited adapter for spinless, fully periodic legacy PythTB models.

## Numerics and units

- Internally, `hbar=m_e=e=1`. Off-diagonal OAM transforms **covariantly** under eigenstate phase changes; diagonal expectations and Kubo contractions are invariant.
- Conductivity coefficients are reported in `e/(2*pi)` in 2D and `e/(2*pi*length)` in 3D. Their SI units are C and C/m, respectively, rather than the S used for electrical conductance.
- Restoring `L/hbar` for an eV/Å model also requires the factor `m_e E0 a0²/hbar²`. Specify the scales with `ohk.units.UnitSystem`.
- Exact degeneracies raise `DegeneracyError` by default. Near-degenerate values are retained and flagged. Broadening is an explicit approximation; there is no default truncation.
- `convergence_study` returns measured sequences and successive differences. It does not automatically claim convergence or extrapolate from a few meshes.
- The complete Haldane model's `d0` term shifts the gap center. `mu=0` does not guarantee a half-filled insulating state. Studies of transition singularities must choose a chemical potential inside the gap.

## Reproducing the paper's key figures

~~~powershell
.\.venv\Scripts\python.exe scripts/reproduce_paper.py --output output/paper-reproduction
~~~

This generates the OAM path in Figure 1, the flux sweep in Figure 2, and the Fermi-energy scan in Figure 3, together with PNG/SVG files, raw CSV data, and convergence tables. By default, it compares a 240×240 Bloch mesh with an independent 60×60 finite periodic lattice. See the [reproduction notes](docs/reproduction.md) for the flux convention, the eV/Å scale assumption, and convergence limits near the transition. In a new environment, first run `python -m pip install -e ".[plot]"`.

## Structure and tests

The [original architecture](docs/design/03-architecture.md) retains the L0-L5 layering. Where historical notes conflict with the implementation, the [implementation corrections](docs/design/05-implementation-notes.md) take precedence.

Tests cover three independent formula paths, three-/four-band negative controls, phase covariance and symmetries, finite lattices, the Chern transition, an independent Kubo double sum, diagonalization reuse in scans, Wannier90 round trips and invalid inputs, units, reports, and the CLI. CI is configured for Linux/Windows and Python 3.10/3.13. External material downloads are not required for offline tests.

MIT; see [LICENSE](LICENSE).
