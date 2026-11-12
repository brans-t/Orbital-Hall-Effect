# 05 · Implementation corrections and validation record

The implementation retains the L0-L5 layers, module boundaries, and pluggable OAM interface from document 03, completing the runnable P0-P5 pipeline. Documents 01, 03, and 04 remain as design history. The following conclusions supersede their original claims after checking the paper and performing independent calculations.

## 1. Denominators and signs in Eq. (16b)

On page 3 of the repository PDF and in [the paper's Eq. (16b)](https://arxiv.org/html/2611.04188v1), the sum has a **positive** sign with denominators **eps_m-eps_n**. The old notes transcribed the denominators as **eps_n-eps_m** while retaining the positive sign, introducing a transcription error rather than identifying a printing error in the paper.

The code uses eps_n-eps_m denominators with a negative sum sign, equivalent to the published expression. `SpectralSum(sign=+1)` is retained only as a negative control for the incorrect transcription.

For two bands, sum(m != n,n') is empty. Flipping its sign therefore changes neither the off-diagonal matrix elements nor the conductivity. Tests detect the transcription error with three-/four-band models and check equivalence for two bands.

## 2. Phase covariance and the Berry connection

Under U -> U D, both v and L satisfy O -> D† O D. Off-diagonal complex phases change; diagonal expectations, off-diagonal absolute values, and Kubo contractions are invariant. Symmetry tests align the band bases before comparing them.

Single-k-point spectral data cannot determine the diagonal Berry connection. The code computes it only when explicit eigenstate derivatives in a specified gauge are provided. The main route uses derivative coefficients with zero diagonal, and the omitted gauge terms cancel in pairs. The old explanation that the cross product of A_n and a real v_nn must vanish is incorrect.

## 3. Pezo, ACA, and the trace

The first term of the published Eq. (16b) is the intermediate-band sum. The Pezo comparison retains this sum and omits the diagonal-velocity term, making off-diagonal elements zero in any two-band model. The original skeleton incorrectly identified the diagonal-velocity term as the first term.

`ACASpec((1,))` denotes one three-dimensional p multiplet; `ACASpec((0,0))` denotes two s orbitals. ACA supplies the entire matrix, so the assembler does not overwrite its zeros with the complete OAM diagonal.

There is no general zero-trace rule for OAM: for example, L11=L22 in a two-band model is usually nonzero.

## 4. Independent finite-periodic-lattice reference

`FiniteLattice` constructs H_site and velocity from the explicit hopping list alone. It does not call the model's H(k), dH(k), or the production OAM formula.

Periodic velocity uses the actual hopping displacement, including boundary winding:

~~~text
v_ij = +i * (R_hop + tau_j - tau_i) * t_ij
~~~

Substituting wrapped site coordinates directly into -i[r,H] produces boundary errors. The static `velocity(H,positions)` helper therefore applies only to open boundaries.

Finite-lattice translation symmetry defines momentum sectors that are diagonalized independently. The real-space velocity is projected into each sector, off-diagonal position is reconstructed as r_nm=-i*v_nm/(eps_n-eps_m), and symmetrized r×v is evaluated directly. The position diagonal is set to zero in the **energy eigenbasis**; removing every position diagonal in the site basis would erase the entire position operator.

Equal-energy states at different momenta are identified by translation quantum numbers. Comparisons project into the same Bloch states rather than guessing correspondence from a global energy sort. This dense implementation is intended for small reference systems.

## 5. Haldane, Chern numbers, and occupations

The complete Γ-point energies are:

~~~text
eps_± = 6*t2*cos(phi) ± sqrt(9*t1² + M²)
~~~

At the default parameters and phi=pi/4, they are [-2.39015309, 3.66294529], not the old description's ±3.0265. Analytic derivatives are cross-checked against central differences.

FHS uses the periodic cell gauge and corrects for the original lattice vectors' negative orientation. The principal arg and Im(log) are equivalent on the same branch. With the default positive-flux convention, the lowest-band Chern number changes from 0 to -1.

The d0 term shifts the Dirac-point energy to -3*t2*cos(phi), so mu=0 is not a universal half-filled, in-gap choice. The transition-divergence test follows this gap center when choosing the chemical potential.

## 6. Kubo normalization and units

For a normalized mesh with sum(w)=1:

~~~text
sigma_raw = (2/A_cell) * sum_k w_k sum_n f_n Im sum_m J_nm v_mn / gap_nm²
~~~

This raw coefficient is in units of e. Expressing it in e/(2*pi) requires another factor of 2*pi. The old notes labeled the raw value of approximately -0.17 as e/(2*pi), omitting this normalization factor.

For a model with an energy scale of E0 eV and a length scale of a0 Å:

~~~text
L/hbar = L_internal * m_e*(E0*e)*(a0*1e-10)²/hbar²
~~~

`UnitSystem` implements this scale. The SI unit of orbital conductivity is C in 2D and C/m in 3D, not the S used for electrical conductance.

Fermi scans cache band-resolved responses, diagonalizing once per k point. Results retain pointwise occupation counts rather than describing a metal with a single integer.

Exact degeneracies raise by default; near-degenerate values are retained and their fraction reported. Enabling eta applies explicit broadening to OAM and Kubo expressions. This is an approximation: different variants need not agree at finite eta, and both eta -> 0 and mesh convergence require checks. The convergence helper returns measured sequences and differences without an unjustified extrapolation.

## 7. Validation evidence

Measured on Windows / Python 3.13 / NumPy 2.3.5:

| Check | Result |
|---|---|
| Automated tests | 90 passed, 90% overall code coverage |
| Ruff and formatting checks | All passed |
| Three formula routes on three-/four-band models | Benchmark residual < 5e-16 |
| Finite lattice versus Bloch | Benchmark residual about 1.33e-15; 12×12 test < 1e-10 |
| Chern numbers at phi/pi=0,.15,.25,.5 | 0,0,-1,-1 |
| Haldane, mu=0, 60×60 mesh, phi=pi/4 | -1.0616491612557013 e/(2*pi) |
| One-command validation at n_k=60 | All checks passed; maximum residual about 1.33e-15 |

These are executable validation results and fixed-mesh regression data, not a claim that arbitrary material predictions are converged. The full test outcome is determined by the actual pytest, Ruff, and CLI output.

## 8. Wannier90 material test and scope

The parser supports 2D/3D lattices, Fortran D exponents, degeneracy weights, block completeness checks, and round-trip output. It does not double conjugate hoppings by default. An hr.dat file does not contain lattice vectors, and their absence raises an explicit error.

The Fe test uses [this pinned WannierTools revision](https://github.com/quanshengwu/wannier_tools/tree/362171ce35450118828c6960a882e3a6d08c73f5/examples/Fe) and is reproduced by `scripts/fetch_fe_example.py`. Third-party data are not versioned.

~~~text
commit: 362171ce35450118828c6960a882e3a6d08c73f5
archive SHA256: 388792dc706600f1b63e57e1c0a8ba06eb32b60d43e24a3c7138d548faee269f
hr.dat SHA256: 80bc1e558e0095532613e2b829ccded2b3e00652abd1f2da15cf47588a2bcd12
num_bands: 16
nrpts: 1367
grid: 3×3×3
mu: 9.1839 eV
eta: 0.01 eV
min_gap: 0.0013646982373032301 eV
sigma: -27.947882312295473 e/(2*pi*model_length)
~~~

This result demonstrates that a real file passes through parsing, eigenbasis construction, OAM, Kubo integration, and JSON output. It is not a converged Fe OHC prediction. The hr-only implementation uses the periodic cell gauge and diagonal-position approximation; it does not implement full Wannier position matrices, anomalous position, or WS-distance corrections. The legacy PythTB adapter is limited to spinless, fully periodic models. Full position matrices, sparse large-system methods, and adaptive integration remain future extensions.
