# 03 · Library architecture

> **Historical design record:** This document preserves the original design. Claims about the Eq. (16b) sign, phases, Pezo truncation, the periodic position operator, units, and some benchmarks were corrected during implementation. [05-implementation-notes.md](05-implementation-notes.md) and the executable tests take precedence.

Goal: a reusable Python library accepting tight-binding or Wannier90 models and producing gauge-invariant OAM matrix elements and orbital Hall conductivity, with a finite-lattice real-space reference for validation.

The original design's first principle is that **$L_{nn'}$ represents a family of definitions rather than one universal formula**. The paper's Eq. (16) is one member; the original claim of a printing error was later corrected in the implementation notes. An $L$ definition must therefore be pluggable rather than hard-coded.

---

## 1. Layers

```
+------------------------------------------------------------------+
| L5 CLI/scripts     ohc-run, oam-scan, validate                   |  <- user entry
+------------------------------------------------------------------+
| L4 Reports         Report: tensors, figures, checks              |
+------------------------------------------------------------------+
| L3 Transport       Kubo: conductivity, Fermi scan                |
+------------------------------------------------------------------+
| L2 OAM core        OAMVariant family and assembler               |  <- library core
+------------------------------------------------------------------+
| L1 Band quantities EigenBasis: eps, v, dH, C_j                   |
+------------------------------------------------------------------+
| L0 Models          Model -> H(k), dH(k)                          |
+------------------------------------------------------------------+
```

Dependencies point strictly downward: L2 uses the public L0/L1 interfaces; L3 consumes $L$ from L2 and $v$ from L1. Each layer can be replaced independently, for example by a WannierBerri L0 backend.

**L1 is the central data contract:** calculations need $\varepsilon$, $v^j_{nm}$, $\partial_{k_j}\mathcal H$, or equivalent $C_j$ at a given $\mathbf k$. All higher layers use this object.

---

## 2. Modules and responsibilities

### L0 `ohk.model`

| Class/function | Responsibility |
|---|---|
| `TightBindingModel` (abstract base) | Required interface: `H(k) -> (N_b,N_b) ndarray`, `dH(k) -> (2,N_b,N_b)`, `lattice()` returning real-space vectors, and `num_bands` |
| `AnalyticModel` | Constructs analytic `H` and `dH` from `(i, j, R, amplitude)` hoppings for hand-built models such as Haldane |
| `Wannier90Model` | Reads `_hr.dat`, or delegates to PythTB `W90`, and supplies analytic `dH` |
| `FiniteLattice` | Builds real-space $\mathcal H$, position, and velocity on an $N\times N$ periodic lattice; **for L5 validation only** |
| `HaldaneModel` | Built-in paper benchmark with the complete $d_0(\mathbf k)$ term; see theory pitfall 3 |

Use analytic `dH` whenever possible. The base class provides `dH_finite_difference` as a documented fallback only when analytic derivatives cannot be supplied.

### L1 `ohk.eigen`

| Class/function | Responsibility |
|---|---|
| `EigenBasis` (dataclass) | One `at(k)` result: `eps`, `U`, `v` ($v^j_{nm}$), and `dH`, with lazy `C_j` caching |
| `derivative_coefficients(eps, v)` | Builds $C_j$; see theory §4.3 |
| `velocity_matrix(H, dH, U)` | $v^j = U^\dagger \partial_{k_j}H\,U$ |
| `berry_connection(eps, v)` | Optional computation of $\mathcal A_n$ for variants that need it; unnecessary on the main route |

**Original gauge contract:** use `U` only to construct quantities described here as gauge invariant: $v$, $L$, and $\sigma$. Downstream code must not rely directly on eigenvector phases. `berry_connection` is explicitly phase sensitive. The original test proposed gauge changes with unchanged $L$; the implementation notes clarify the covariance law.

### L2 `ohk.oam`: the core

This layer is intended to distinguish the library from existing tools.

```
oam/
├── base.py           OAMVariant abstract interface
├── covariant.py      Eq. (16a), projected form; main route
├── spectral.py       Eq. (16b), spectral form with corrected transcription
├── atomic.py         ACA, atom-centered approximation; negative control
├── pezo.py           Pezo-Manchon truncation; paper Ref. 26
└── assemble.py       Complete L matrix: diagonal and off-diagonal parts
```

Abstract `OAMVariant` interface:

```python
class OAMVariant(ABC):
    name: str
    @abstractmethod
    def matrix(self, basis: EigenBasis) -> np.ndarray:
        """Return complex L[n,n'] of shape (N_b, N_b), with hbar=m_e=1."""
```

The four variants differ only in their `matrix` implementation:

| Variant | Expression | Purpose |
|---|---|---|
| `CovariantProjector` | Eq. (16a), matrix form from §4.3 | **Main route**, originally described as gauge invariant; $O(N_b^3)$ |
| `SpectralSum` | Eq. (16b), proposed corrected sign and closed-form diagonal | Cross-check; faster for small $N_b$ |
| `AtomCentered` | $\sum_i \delta_{\text{cell}}\,\mathbf L_i$ | Negative control: identically zero for s orbitals; the paper's key §IV comparison |
| `PezoTruncated` | Retains only the first term of Eq. (16b) | Reproduces the two-band $L_{12}\equiv0$ limitation discussed in the paper |

**Why implement inadequate variants:** the paper compares differing literature results. Implementing them enables quantitative comparisons on the same model, for example $L_{12}=0$ for PezoTruncated on Haldane versus about 0.04 on the main route. This provides research value and a strong regression check.

**Assembly:** `assemble.full_matrix()` supplies the closed-form diagonal and variant off-diagonal elements, then checks Hermiticity. It raises when the norm of `L - L†` exceeds tolerance.

### L3 `ohk.transport`

| Class/function | Responsibility |
|---|---|
| `orbital_current(L, v, axis)` | Symmetrized $J_x = \frac12(v_xL_z+L_zv_x)$ |
| `KuboConductivity` | Consumes `L_variant`, velocity, and occupations; returns $\sigma^{\rm OH}_{xy}$ |
| `fermi_scan(L_variant, energies)` | Scans $E_F$ and returns an array |
| `sum_rule_check` | Original proposal: check $\mathrm{tr}\,L=0$, Hermiticity, and symmetries |

Separate integration from mesh generation. `BZGrid` supplies $\mathbf k$ points and weights, including uniform Monkhorst-Pack, Haldane geometry, and custom meshes. This supports reproducing Figures 2-3.

### L4 `ohk.report`

Package tensors, checks, and optional figures in `Report`. Separate computation from presentation: `Report` consumes pure data without importing the optional matplotlib dependency, allowing CI to validate calculations without plotting.

### L5 `ohk.cli`

Three commands: `ohc-run` computes a model's $\sigma$, `oam-scan` exports $L_{nn'}$ along a k path, and `validate` runs consistency checks. CLI code is a thin wrapper over L0-L4.

---

## 3. Data flow for an OHC calculation

```
model (L0)
   │  for each k in BZGrid
   ▼
EigenBasis.at(k)  ──► eps, U, v, dH        (L1)
   │
   ▼
variant.matrix(basis) ──► L (N_b×N_b)      (L2)
   │
   ▼
orbital_current(L, v) ──► J_x              (L3)
   │
   ▼
Kubo sum over occupied n, all m  ──► σ(k)  (L3)
   │
   ▼
Σ_k weight·σ(k) / A  ──► σ_OH             (L3)
   │
   ▼
Report(σ, checks, metadata)               (L4)
```

All variants at a k point reuse one `EigenBasis`. Only one diagonalization is required, making side-by-side comparison inexpensive.

---

## 4. Tensor and unit conventions

Perform **no unit conversion inside the computation**. Set $\hbar=m_e=e=1$ and use lattice lengths. Repeated $m_e/\hbar$ factors are easy to mix or omit; conversion to SI with $e=\hbar=1$ is performed only in the output layer.

| Quantity | Internal convention | Output conversion |
|---|---|---|
| $L$ | In $\hbar$ units, `hbar=1` | Multiply by $\hbar$ |
| $v$ | $[\text{energy}]\cdot[\text{length}]$, `hbar=1` | Multiply by $1/\hbar$ |
| $\varepsilon$ | Natural model units | Multiply by the energy unit |
| $\sigma^{\rm OH}$ | $e/2\pi$ units in 2D, with cell area $A$ | Multiply by $e/2\pi$ or $e/(2\pi a)$ |

Index conventions: `L[n, n']`, `v[j, n, m]` with $j\in\{x,y\}$, and `dH[j, n, m]`. **Use `[direction, band, band]` for $v$ throughout**, matching `EigenBasis` and preventing confusion between $v_{nm}$ and $v_{mn}$.

---

## 5. Numerical strategy

**Diagonalization:** `numpy.linalg.eigh` once per k point is practical for $N_b\le 200$. Larger systems may use `scipy.linalg.eigh` with `subset_by_index` for a window near occupied bands.

**Degeneracy handling:** denominators $\varepsilon_n-\varepsilon_m$ in $L_{nn'}$ are singular at band degeneracies. The original strategy is:

1. Compute the minimum gap $\delta_{\min}$ at each k point.
2. Mark a point `degenerate` if $\delta_{\min}<\tau_{\text{deg}}$, with default $10^{-8}$.
3. Report the degenerate-point fraction and warn above the default 0.1% threshold.
4. Optionally regularize within a degenerate subspace by diagonalizing a combination of $\mathbf v$ and truncating inverse eigenvalue differences. **This is not the default** because it changes the quantity's definition; report the issue by default.

**Near a transition:** the $\sigma$ divergence at a Chern change in Figure 2 is physical. Flag points near a topological transition when phi approaches $\phi_c$ rather than flattening peaks.

**k-mesh convergence:** the original design proposes `convergence_study(grid_sizes)` returning $\sigma$ versus mesh size and an extrapolation. It describes the Figure 2 line as the $k$ continuum limit and circles as a $60\times60$ lattice. Record their difference quantitatively.

**BZ integration:** use equal-weight Monkhorst-Pack by default. For Haldane's honeycomb lattice, construct reciprocal vectors as $2\pi A^{-T}$ from real-space vectors. Test $a_i\cdot b_j=2\pi\delta_{ij}$ explicitly.

---

## 6. Extension points

| Extension | Interface | Future application |
|---|---|---|
| Model format | Subclass `TightBindingModel` | TBmodels HDF5, pyqula, noncollinear spin |
| New $L$ definition | Subclass `OAMVariant` | Sastges-Baek and Go-Lee expressions including anomalous position |
| Response function | Add a `KuboX` subclass in L3 | Orbital Edelstein effect, orbital magnetization, spin Hall comparisons |
| Integrator | Replace `BZGrid` | Adaptive meshes, tetrahedra, WannierBerri backend |
| Backend | Implement `TightBindingModel` for a WannierBerri `System` | Large first-principles systems |

These interfaces distinguish a reusable engine from a one-off script.

---

## 7. Package layout

```
orbital-hall/
├── pyproject.toml
├── README.md
├── docs/design/                    01-theory · 02-tooling-survey · 03-architecture · 04-validation
├── src/ohk/
│   ├── __init__.py                 __version__, public exports
│   ├── units.py                    hbar, m_e, conversion factors
│   ├── model/
│   │   ├── base.py                 TightBindingModel, ModelSpec
│   │   ├── analytic.py             AnalyticModel: hoppings -> H, dH
│   │   ├── haldane.py              HaldaneModel benchmark
│   │   ├── wannier90.py            Wannier90Model (_hr.dat)
│   │   └── finite.py               FiniteLattice oracle
│   ├── eigen.py                    EigenBasis, derivative_coefficients
│   ├── oam/
│   │   ├── base.py                 OAMVariant
│   │   ├── covariant.py            Eq.(16a)
│   │   ├── spectral.py             Eq. (16b), corrected transcription
│   │   ├── atomic.py               ACA negative control
│   │   ├── pezo.py                 Pezo truncation
│   │   └── assemble.py             Diagonal/off-diagonal assembly
│   ├── transport/
│   │   ├── current.py              Orbital-current operator
│   │   ├── kubo.py                 KuboConductivity, fermi_scan
│   │   └── grids.py                BZGrid
│   ├── validation/
│   │   ├── checks.py               Hermiticity, trace, symmetries
│   │   ├── oracle.py               Finite-lattice comparison
│   │   └── benchmarks.py           Haldane benchmark values
│   ├── report.py                   Report
│   └── cli.py
└── tests/
    ├── conftest.py                 Haldane / random-model fixtures
    ├── test_model.py
    ├── test_eigen.py
    ├── test_oam_consistency.py     Three-path comparison; core test
    ├── test_symmetry.py            Time-reversal / inversion constraints
    ├── test_oracle.py              Finite lattice versus Bloch
    ├── test_transport.py           OHC sign and scale
    └── test_units.py
```

Naming: `ohk`, for orbital Hall kit, is short and distinct. Modules are organized around physical objects rather than mathematical operations so users can find code through physical intuition.

---

## 8. Test strategy as part of the architecture

Visual agreement with a paper figure does not establish correctness. Tests rely on **agreement between independent computational paths**:

| Test | Scope | Original expectation |
|---|---|---|
| `test_oam_consistency` | Same k point: `CovariantProjector`, `SpectralSum`, and pointwise definition | Difference $<10^{-12}$ |
| `test_oracle` | Finite-lattice real-space $L$ versus Bloch $L$ | Difference $<10^{-10}$ at allowed momenta |
| `test_symmetry` | $L_{11}(k)=-L_{11}(-k)$ at $\phi=0$; $|L_{12}|$ vanishes at $\phi=\pi/2$ | $<10^{-12}$ |
| `test_degenerate_bands` | Two bands cannot test the $L_{12}$ sum term; use a random three-band model | Three-band residual $<10^{-12}$ |
| `test_transport` | Calibrate the sign of $\sigma^{\rm OH}$ | Historical claim: proposed correction negative, printed form positive |
| `test_chern` | Haldane Chern number changes at $\phi_c$ | $0\to-1$ |
| `test_units` | Scale of $L$ after restoring $\hbar$ | Comparable to $\hbar$ |

**Essential regression:** random three-band tests in `test_degenerate_bands` are required because two-band models are blind to the Eq. (16b) sum. Testing Haldane alone misses the sign difference described in the original design. Keep this test in CI; the implementation notes correct the original attribution of that difference to the paper's printing.

---

## 9. Staged implementation plan (original design only)

| Stage | Scope | Acceptance |
|---|---|---|
| P0 | L0 models, L1 `EigenBasis`, and tests | Haldane bands and Chern number agree with the paper |
| P1 | L2 `CovariantProjector` and finite-lattice oracle | Three-path agreement $<10^{-12}$ |
| P2 | Remaining L2 variants and comparison | PezoTruncated gives $L_{12}=0$ on Haldane |
| P3 | L3 Kubo and BZGrid | $\sigma^{\rm OH}(\phi)$ reproduces Figure 2's shape and polarity |
| P4 | L4 Report and L5 CLI | `validate` runs the consistency suite |
| P5 | Wannier90 integration | End-to-end calculation with a real-material `_hr.dat` |

The original delivery was a skeleton: interfaces, dataclasses, abstract bases, and module docstrings with `NotImplementedError` bodies. P0 would fill these bodies in a subsequent stage without changing the architecture.
