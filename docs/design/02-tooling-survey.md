# 02 · Survey of existing computational tools

The design-stage conclusion was that **no existing tool directly supplies and validates the complete set of $L_{nn'}(\mathbf k)$ matrix elements defined by the paper's Eq. (16)**. The proposed approach is to use mature tight-binding/Wannier tools for $\mathcal H(\mathbf k)$ and $\partial_{k_j}\mathcal H(\mathbf k)$, and implement $L$ and the OHC layer locally. This is the original survey, not a current software-capability audit.

## 1. Ecosystem overview

The Wannier community maintains the self-describing `wannier-developers/wannier-ecosystem-registry`, accompanied by review arXiv:2312.10769, categorizing Wannier-related software by capability. The survey found Kubo AHE/SHE tools under linear response, but **no registry entry explicitly listing orbital Hall effects or OAM matrix elements**. This was a central finding of the original survey.

## 2. Individual assessments

### 2.1 Wannier90

Purpose: construct maximally localized Wannier functions from first principles, producing `_hr.dat` (real-space $H(\mathbf R)$), `_mmn`, `_eig`, `_chk`, `_uHu`, `_uIu`, and related files.

Project value: **an important input source, rather than the OAM engine**. The original design treats $H(\mathbf R)$ in `_hr.dat` as sufficient input for $\mathcal H(\mathbf k)$ and its analytic derivatives. Berry postprocessing through `wannier90.x -berry` or `postw90` provides AHE, SHE, and Berry curvature.

Main limitation identified by the survey: the official Berry module **does not provide the required OHE current vertex**. Its built-in vertices are spin or charge currents; the orbital current $J_x=\frac12(v_xL_z+L_zv_x)$ must be added. `_uHu` stores $\langle u_{mk}|\hat H|u_{n,k+b}\rangle$ for orbital magnetic-moment calculations, which is not identical to the paper's $L$.

Practical note: `uHu_formatted = .true.` requests text `_uHu` output for checking its layout. The survey found no complete column-by-column format specification and recommended checking the generating source code.

### 2.2 WannierBerri

Purpose: high-performance Wannier interpolation and BZ integration, with Berry-curvature, AHC, SHC, and orbital-magnetization calculators.

Project value: **a reference implementation and potential backend**.

- `wannierberri.utils.mmn2uHu` generates `_uHu` from `_mmn` and `_eig`, documented as matrix elements needed for orbital magnetic moments. Its $\langle u|\mathbf r\times\mathbf v|u\rangle$ modern-theory route is related to the gauge-covariant viewpoint discussed in §2.2 of the theory notes.
- The `spin_Hall` calculator accepts `spin_current_type` (`"ryoo"` or `"qiao"`) and supports symmetrization. The intended OHE pattern is **a general calculator with an injected, symmetrized current vertex**.
- It separates formulas (`formula/covariant.py`) from calculators (`calculators/`) and abstracts the `spin_Hall` current operator $J$. This separation directly informed the proposed architecture.

**Existing OHC support:** The survey did not find an orbital Hall calculator listed in the WannierBerri documentation. It identified an orbital-moment route for $L$ and a Kubo route for SHC, but no documented layer joining them into the OHC current vertex. The proposed contribution is **that connecting layer**, rather than another interpolation engine.

### 2.3 WannierTools

Purpose: topological properties of Wannier90 tight-binding models, including surface states, Wilson loops, and AHC/SHC.

The survey records that **GitHub release notes for v2.7.2, around June 2025, state "Added orbital Hall effects calculations"**. It identified this as an existing public OHE capability. Inputs are `wannier90_hr.dat` and `wt.in`, using NAMELIST and INPUT_CARDS formats.

Project value and positioning: **a comparison target as well as a related tool**. Its existence supports demand for OHE calculations in Wannier workflows. The original assessment questioned whether its Fortran and fixed-file workflow was suitable as the sole independent check of the paper's formulas, describing the relevant implementation as unavailable for inspection. The proposed role of this project is **a Python reference for Eq. (16) with a finite-lattice oracle**. WannierTools may provide a third-party cross-check when its vertex definition matches the paper.

### 2.4 PythTB

Purpose: Python tight-binding model construction and topological/quantum-geometric analysis.

- The survey describes `tb_model`/`TBModel`/`W90` interfaces, with `W90("path").model` reading Wannier90 output, `W90._read_hr_dat` parsing `_hr.dat`, and `min_hopping_norm` pruning small hoppings.
- Interfaces such as `_model.get_hoppings()` expose hopping lists, while `lattice()`/`sites()` provide positions.

Limitation identified by the survey: **no built-in OHE/OAM capability**, or general Kubo transport. Its value is **convenient model definition**, especially hand-built Haldane models and user-supplied models. The `W90` reader can inform, or potentially supply, the project's `_hr.dat` reader.

### 2.5 TBmodels

Purpose: build tight-binding models from Wannier90 output and store them in HDF5.

Value: `Model` reads `*_hr.dat` with `write_hr = .true.` enabled in Wannier90 and provides `eval`, `bands`, and `wannier_centers`. It is a lightweight file-to-model adapter. The original design preferred **PythTB's `W90` route** over adding both adapters, because PythTB also provides the hand-built-model API needed for the Haldane benchmark.

### 2.6 Other related codes

| Code | Relationship to OHE | Assessment |
|---|---|---|
| **Rustb** (docs.rs/Rustb) | Documentation states that it calculates orbital Hall conductivity and orbital angular momentum | OHC support, but Rust and limited documentation leave input conventions unclear; a potential third-party check |
| **pyqula** (joselado) | Python TB library with Kubo/Hall linear response | Kubo implementation may be informative; no OAM vertex identified in the survey |
| **TBPLaS 2.0** | Large-system Kubo-Bastin Hall transport | Real-space method complements the project's k-space route, especially for strong disorder |
| **TightBinding++** | Kubo-Greenwood conductivity | Another methodological reference |
| **Go & Lee, PRB 109, 174435 (2024)** | Wannier-interpolated OHC emphasizing orbital dependence of anomalous position | **No public repository found in the survey**; its conclusion that orbital-dependent position corrections matter supports gauge-covariant treatment |

### 2.7 Parallel expressions in the literature

The survey identified related expressions for $L$ that **differ from one another**, matching the conflicting-results discussion in the paper's introduction:

- **Pezo, Ovalle, Manchon, PRB 106, 104414 (2022)** [paper Ref. 26]: retains only the first term of Eq. (16b) for off-diagonal elements. The paper states that **this gives identically zero $L_{12}$ for every two-band system**, while direct Haldane calculations give $L_{12}\ne0$.
- **Cysne, Souza, Rappoport, PRR 8, 033086 (2026)** [paper Ref. 27]: the paper reports that the final expression agrees with its Eq. (16).
- **Sastges & Baek, arXiv:2604.08280**: another Wannier-basis OAM construction based on modern orbital-magnetization theory.

Architecture implication: **the choice of formula must be explicit, rather than hard-coded**. The library should support several $L$ definitions, initially named paper-16a, paper-16b-corrected, Pezo-style-truncated, and Rappoport-style, and compare them using the same $\mathbf v$ and $\varepsilon$. This separates a reusable library from a one-off script. See the architecture document's proposed `oam_variants` component.

## 3. Capability gaps

| Capability | Existing option | Gap identified by the survey |
|---|---|---|
| $\mathcal H(\mathbf k)$ and analytic $\partial_k\mathcal H$ | Wannier90 `_hr.dat`, PythTB, TBmodels | None |
| Wannier interpolation, BZ integration, k-point parallelism | WannierBerri | None; backend reuse is possible |
| Orbital magnetic moment $\langle u|\mathbf r\times\mathbf v|u\rangle$ | WannierBerri `mmn2uHu`/orbital-magnetization route | Band-diagonal contribution only |
| **Complete OAM matrix, including off-diagonal elements** | **No general implementation identified** | **Core contribution of this library** |
| **Orbital-current vertex $J_x=\frac12(v_xL_z+L_zv_x)$** | **None identified** | **Core contribution of this library** |
| **OHC Kubo integration** | WannierTools, described in the survey as a closed implementation, and Rustb | An independent reference implementation |
| **Finite-lattice real-space oracle** | None identified | Project validation layer |
| Side-by-side comparison of $L$ definitions | None identified | A distinguishing design feature |

## 4. Proposed technology choices

Start with `Python 3.10+`, `numpy`, and `scipy`. Keep `PythTB` optional for `_hr.dat` input and hand-built models, avoiding a dependency on its data structures. Implement initial BZ integration locally for Haldane's nonstandard geometry and custom meshes, then consider WannierBerri for larger systems. No C++/Fortran kernel is proposed initially: NumPy $O(N_b^3)$ products are adequate for $N_b\lesssim100$, and correctness comes first.

Use `pytest` for tests and optionally `h5py` for caching large-model $\mathcal H(\mathbf k)$ samples.
