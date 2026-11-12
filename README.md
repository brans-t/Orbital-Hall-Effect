# Orbital Hall kit (ohk)

Compute orbital angular momentum (OAM) and orbital Hall conductivity (OHC) from
tight-binding or Wannier90 Hamiltonians. Run a bundled example, edit a JSON input,
and find numerical tables and figures in one output directory.

## 1. Install the environment

Requirements: **Python 3.10 or later**. NumPy and SciPy are required; the `plot`
extra installs Matplotlib for PDF/PNG figures. No GPU, external service or model
download is needed for the first example.

Run these commands from the repository root.

### Windows / PowerShell

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[plot]"
.\.venv\Scripts\python.exe examples/first_run.py
```

If `.venv` already exists, use its interpreter and skip the first command. These
commands do not require activating PowerShell scripts.

### Linux / macOS

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[plot]"
.venv/bin/python examples/first_run.py
```

For table-only use, install with `python -m pip install -e .` and select
`--formats csv`. For development and all tests, install `".[dev,plot]"` instead.
Use your environment's Python in the commands below. The installed commands
`ohc-run`, `oam-scan` and `ohk-validate` are in `.venv/Scripts/` on Windows and
`.venv/bin/` on Linux/macOS; activating the environment makes them available by name.

## 2. First run: gapped graphene

```powershell
.\.venv\Scripts\python.exe examples/first_run.py
```

The example reads [inputs/graphene.json](inputs/graphene.json). It uses the Haldane
implementation with nearest-neighbor hopping only (`t2=0`), two orbitals per unit
cell, hopping `2.8 eV`, onsite energies `+/-0.25 eV`, and nearest-neighbor distance
`1.42 Angstrom`. The band gap is `0.5 eV`.

It calculates the conductivity at `mu=0` and a Fermi-energy scan from `-1.5` to
`1.5` at 121 points, on a `120 x 120` momentum mesh. No network access is needed.
At the supplied mesh, the internal conductivity at `mu=0` is approximately
`-47.18951595`. This is a reproducible example value, not an automatic convergence
certificate or a physical-unit conversion.

Open **`output/first-run/conductivity.pdf`** to see the curve. The terminal prints
absolute paths to every file after a successful run.

```text
inputs/
  graphene.json                 Editable, self-contained example input
  wannier90.example.json        Template for your material data
  wannier90/                    Your hr.dat and lattice files (Git-ignored)
examples/
  first_run.py                  One-command beginner example
output/                         Generated results (Git-ignored)
  first-run/
    conductivity.csv            Fermi energy and calculated OHC
    conductivity.pdf            Vector figure of those same values
    input.json                  Resolved inputs; can be run again
    report.json                 Results, diagnostics and output-file list
```

The PDF shows the sampled calculation, with no smoothing. Stair steps away from
the gap can reflect the finite momentum mesh. A finer mesh is needed before
interpreting detailed structure.

## 3. Choose outputs and their location

Default formats are **CSV and PDF**. PNG is available for a quick image preview.

```powershell
# CSV and PDF (default)
.\.venv\Scripts\python.exe examples/first_run.py

# CSV only; does not require Matplotlib
.\.venv\Scripts\python.exe examples/first_run.py --formats csv

# PDF only
.\.venv\Scripts\python.exe examples/first_run.py --formats pdf

# CSV, PDF and PNG, in a custom directory
.\.venv\Scripts\python.exe examples/first_run.py --formats csv pdf png --output-dir output/my-run

# Smaller mesh for a faster functionality check, not a converged calculation
.\.venv\Scripts\python.exe examples/first_run.py --n-k 24 --output-dir output/quick-check
```

`report.json` and the resolved `input.json` are always saved by the conductivity
workflow, even when CSV is disabled. The report retains the numerical Fermi scan,
grid size, operator variant, broadening, minimum sampled band separation and
near-degenerate fraction. Its `metadata.outputs` lists the files from that run.

Running again in the same directory replaces files of the same names. Older
formats from earlier runs are not deleted; use the printed list/report to identify
current outputs, or choose a new directory for each study. Source model files under
`inputs/` are kept separate from these generated files.

## 4. Edit a model input

With the environment active, the first-run script is equivalent to:

```bash
ohc-run --config inputs/graphene.json
```

Configs are flat JSON objects. These are the supported settings:

| JSON key | Meaning | Default without a config |
|---|---|---|
| `model` | `haldane` or `wannier90` | `haldane` |
| `t1`, `t2` | Nearest- and next-nearest-neighbor hopping | `1.0`, `0.15` |
| `mass` | Opposite onsite energies `+mass` and `-mass` | `0.4` |
| `phi` | Next-nearest-neighbor hopping phase, in radians | `pi/4` |
| `nearest_neighbour_distance` | Haldane nearest-neighbor distance | `1.0` |
| `hr`, `lattice` | Wannier90 data and lattice-matrix paths | Required for Wannier90 |
| `variant` | `covariant`, `spectral`, `pezo`, or `aca` | `covariant` |
| `n_k` | Grid points per periodic direction | `40` |
| `mu` | Fermi energy for the reported single-point OHC | `0.0` |
| `fermi_range` | Inclusive scan: `[start, stop, count]` | No scan |
| `fermi` | Explicit Fermi-energy array; alternative to `fermi_range` | No scan |
| `convergence` | Additional mesh sizes, e.g. `[24, 48, 96]` | No mesh study |
| `eta` | Positive opt-in broadening for OAM and Kubo denominators | Unbroadened |
| `formats` | Array containing `csv`, `pdf`, and/or `png` | `["csv", "pdf"]` |
| `output_dir` | Destination for all files | `output/ohc/` |

Unknown keys, invalid numeric types and invalid scan ranges are rejected. Paths
inside a config are relative to **that JSON file**. Command-line paths are relative
to **the current working directory**. CLI flags override JSON values; specifying
`--fermi` replaces a configured `fermi_range` and vice versa.

```bash
ohc-run --config inputs/graphene.json --n-k 180 --output-dir output/graphene-180
ohc-run --config inputs/graphene.json --fermi-range -2 2 161 --formats pdf
ohc-run --config inputs/graphene.json --convergence 24 48 96 --output-dir output/mesh-study
ohc-run --config output/first-run/input.json --output-dir output/repeated-run
```

A mesh study also writes `convergence.csv` and/or `convergence.pdf`. It reports
the measured sequence without declaring convergence or extrapolating it.

The CLI has built-in Haldane/graphene and Wannier90 inputs. General hopping models
are available through the Python API; arbitrary hopping lists are not a JSON CLI
input format in this release.

## 5. Supply a Wannier90 material model

An `_hr.dat` contains the real-space Hamiltonian matrices. The lattice is supplied
separately as a whitespace-delimited `2 x 2` or `3 x 3` matrix, with primitive
vectors as rows in the model's length unit. For example, a 2D lattice file is:

```text
2.0 0.0
1.0 1.7320508075688772
```

Place your material files under `inputs/wannier90/`, copy/edit
[inputs/wannier90.example.json](inputs/wannier90.example.json), and set the correct
`mu`, lattice, mesh and any explicit broadening. Then run:

```bash
ohc-run --config inputs/wannier90.example.json
```

The template intentionally contains no material Hamiltonian. It is not runnable
until you supply those files. Alternatively, use explicit paths:

```bash
ohc-run --hr inputs/wannier90/wannier90_hr.dat --lattice inputs/wannier90/lattice.txt --mu 0 --n-k 12 --output-dir output/material --formats csv pdf
```

For an optional public Fe functionality check:

```bash
python scripts/fetch_fe_example.py
ohc-run --hr .local/fe/wannier90_hr.dat --lattice .local/fe/lattice.txt --mu 9.1839 --n-k 3 --eta 0.01 --output-dir output/fe-smoke --formats csv
```

The downloader verifies pinned checksums and uses `.local/fe/`. This is a coarse
smoke test; it does **not** give a converged Fe conductivity.

## 6. Other commands and the Python API

```bash
# Momentum-path bands and OAM: output/oam/oam.csv, oam.pdf, report.json
oam-scan --points-per-segment 80 --output-dir output/oam --formats csv pdf

# Custom path: vertices are fractional reciprocal coordinates
oam-scan --vertices "0,0;0.333333,0.666667" --formats pdf

# Independent numerical checks: output/validation/report.json
ohk-validate --n-k 30
```

The default OAM path is Gamma-K-Kprime in 2D, or Gamma-X in 3D. The CSV contains
Cartesian momenta, path distance, every band energy, and real/imaginary parts of
every OAM matrix element. The PDF shows all bands, all diagonal OAM elements, and
the Frobenius norm of the off-diagonal OAM matrix. Its JSON also retains the full
table when CSV is disabled.

Legacy commands `ohc-run --output result.json` and `oam-scan --output path.csv`
remain supported as single-file exports. Use `--output-dir` with `--formats` for
the new workflow. Parent output directories are created automatically.

For a custom model, use `AnalyticModel(ModelSpec(...), hoppings, onsite=...,
orbital_positions=...)` with `Hopping` entries. A minimal calculation is:

```python
from ohk import HaldaneModel, CovariantProjector, BZGrid, orbital_hall_conductivity

model = HaldaneModel(t1=1.0, t2=0.0, M=0.4)
grid = BZGrid.uniform(model.lattice(), n=60)
result = orbital_hall_conductivity(model, CovariantProjector(), grid, mu=0.0)
print(result.sigma, result.units)
```

The Python API returns arrays/results without automatically writing files. The
CLI/example handles exports. Use `EigenBasis.at` and `full_matrix` for full OAM
matrices; calling a variant's `matrix` directly can return only its off-diagonal
part. `k` is Cartesian, lattice vectors are rows, and `v[j,n,m]` uses the direction
as its first index. General model orbital positions default to zero; supply the
intended embedding explicitly. A limited spinless, fully periodic legacy PythTB
adapter is also available.

## 7. Units, scope and troubleshooting

- **Internal units:** `hbar=m_e=e=1`. Output columns named
  `sigma_oh_internal` are numerical coefficients in internal `e/(2*pi)` in 2D
  or `e/(2*pi*length)` in 3D. The example CSV is not already SI-calibrated.
- **Physical conversion:** for energies in eV and distances in Angstrom, use
  `ohk.units.UnitSystem`. Its `oam_scale` converts internal OAM to `L/hbar` and
  also multiplies the OHC coefficient. At unit scales it is about `0.131234212`.
  `conductivity_to_si(..., dimension=2 or 3)` returns C or C/m, not Siemens.
  There is no automatic extra spin-degeneracy factor.
- **Exact degeneracy:** default isolated-band expressions raise
  `DegeneracyError`. `--eta` enables an explicit approximation; examine eta and
  mesh dependence rather than treating finite output as validation.
- **Wannier scope:** H(R)-only input uses a diagonal-position approximation. It
  does not load full `_r.dat` anomalous-position matrices or WS-distance
  corrections. Finite numbers are not a complete first-principles validation.
- **Different OAM definitions:** `pezo`, `aca` and `covariant` need not predict
  the same observable. CLI ACA is only defined for the built-in s-orbital model;
  use an explicit `ACASpec` in Python for other orbital multiplets.
- **No PDF produced:** install `".[plot]"` in the same environment used to run
  the command, or request `--formats csv`. The program checks this before computing.
- **Missing input file:** paths in JSON are relative to that JSON's directory.
  The Wannier template needs your own Hamiltonian and lattice files.
- **Wrong output directory:** inspect the absolute paths printed at the end or
  `metadata.outputs` in `report.json`.
- **Large changes with mesh size:** run `--convergence`; a near-degeneracy
  fraction of zero does not establish numerical convergence.

## Development, architecture and paper reproduction

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest
python -m ruff check src tests examples scripts
python -m ruff format --check src tests examples scripts
```

```text
src/ohk/model/        L0: Hamiltonians, hopping models, Wannier90, finite lattices
src/ohk/eigen.py      L1: Eigenbasis and velocity matrices
src/ohk/oam/          L2: OAM definitions and matrix assembly
src/ohk/transport/    L3: BZ integration, currents and Kubo response
src/ohk/validation/  L4: Independent numerical reference calculations
src/ohk/cli.py        L5: Commands connecting inputs, computations and exports
src/ohk/config.py     JSON validation and reproducible resolved inputs
src/ohk/outputs.py    CSV / PDF / PNG presentation of computed data
tests/               Offline numerical and workflow regression tests
scripts/             Optional material download and paper reproduction
docs/                Theory, architecture, implementation and reproduction notes
```

The original L0-L5 numerical architecture is retained. Historical design notes
are kept for provenance; the [implementation notes](docs/design/05-implementation-notes.md)
take precedence where they differ. The redundant JSON-only example was replaced
by `examples/first_run.py`; the unused `io` installation extra was removed.

Based on Gautam & Satpathy,
[Phys. Rev. B 114, 214403 (2026)](https://doi.org/10.1103/pd3s-q2z7).
To regenerate the three original Haldane figure studies:

```bash
python scripts/reproduce_paper.py --output output/paper-reproduction
```

This separate research script exports PNG/SVG figures, CSV data and audit reports.
Its flux conventions, physical-scale assumptions and convergence limits are
documented in [docs/reproduction.md](docs/reproduction.md). General-purpose PDF
exports use the commands described above.

The figure-export workflow was informed by Scientific Agent Skills:
Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026),
*Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*,
[arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065).

MIT; see [LICENSE](LICENSE).
