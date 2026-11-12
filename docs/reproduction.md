# Reproducing the paper's key figures

The target is Gautam & Satpathy, *Phys. Rev. B* **114**, 214403 (2026),
[DOI](https://doi.org/10.1103/pd3s-q2z7), and the public
[arXiv text](https://arxiv.org/html/2611.04188v1). The code computes the theoretical model directly; it does not trace, digitize, or fit the authors' curves.

## Running the calculation

From the project root:

~~~powershell
.\.venv\Scripts\python.exe scripts/reproduce_paper.py --output output/paper-reproduction
~~~

In a new environment, first run `python -m pip install -e ".[plot]"`. The default calculation takes a few minutes, mostly for the full phase sweep. No online download of the original figures is needed; the source PDF is already in the repository. Outputs include 300 DPI PNG, SVG, CSV, convergence tables, and a JSON manifest.

## Model and required conventions

- Parameters are `t1=1`, `t2=0.15`, and `M=0.4`. Figure 1 uses `phi/pi=0,0.15,0.25,0.5`.
- Positive flux in the paper gives a lowest-band Chern number of `+1`; positive flux in the library's original convention gives `-1`. The script explicitly uses `library_phi=-paper_phi` without changing the model's existing convention.
- The path is K-Γ-K′. In this project's coordinates, K has fractional reciprocal coordinates `(1/3,2/3)`, and K′ is its negative.
- Figure 2 occupies the lowest band at every k point: half filling. It is not a scan at fixed `mu=0`. The `d0` term shifts the gap center; the corresponding central chemical potential, `-3*t2*cos(phi)`, is also saved in the CSV.
- Figure 3 scans the Fermi energy at `phi=pi/2`, using the strict zero-temperature occupation `eps<E_F`.
- The energy unit is explicitly assumed to be **1 eV**, and the nearest-neighbor distance **1 Å**. The physical scale factor restoring `L/hbar` and the conductivity coefficient is `m_e*E0*a0²/hbar² = 0.13123421196457827`. This natural interpretation agrees with both Figures 1 and 3, but **neither the main text nor the public TeX source explicitly confirms these two scales**. Absolute-amplitude reproduction is therefore conditional on this assumption; it is not presented as an input explicitly supplied by the authors. `--energy-ev` and `--length-angstrom` change this assumption. Raw internal-unit data are retained in the CSV files.
- Figures 2 and 3 plot **negative conductivity**, following the paper's vertical axes. There is no additional fitted coefficient, smoothing, broadening, or artificial suppression of peaks.

## Two independent calculation routes

Lines use analytic hopping Fourier sums, batched eigendecomposition, and the projected form of Eq. (16a). The batched route is cross-checked against the established single-point public API on three- and four-band models.

Circles come from an actual **60×60 finite periodic lattice**. Sparse real-space Hamiltonian and velocity matrices are built with boundary windings, then exactly decomposed into translation-momentum sectors by a discrete Fourier transform. OAM is obtained from off-diagonal position matrix elements and the symmetrized `r×v` definition. This route does not call the Bloch model's `H(k)`, `dH(k)`, or projected OAM formula. It agrees with the original dense oracle on small lattices. The FFT is an exact symmetry reduction of the finite lattice.

The default Bloch route uses a **240×240** mesh; the phase plot also shows **120×120** results. Both are finite-mesh approximations, not claimed analytic continuum limits.

## Outputs and interpretation

| File | Contents |
|---|---|
| `fig1_oam.png/.svg` | Off-diagonal and diagonal OAM at four fluxes, with an additional full-range diagonal panel |
| `fig2_phase.png/.svg` | Half-filled conductivity versus flux, with the paper's display window beside a full logarithmic view |
| `fig3_fermi.png/.svg` | Fermi-energy scan and insulating plateau |
| `fig1_bloch.csv`, `fig1_finite.csv` | Path data and independent finite-lattice checks |
| `fig2_phase.csv` | Raw conductivity and minimum gap for every phase and mesh |
| `fig3_fermi.csv` | Fermi-energy scan in both internal and physical output units |
| `convergence.csv` | Measured sequences for meshes of 30, 60, 120, 240, and 480 |
| `manifest.json` | Parameters, conventions, software versions, source/PDF hashes, and calculation diagnostics |

The paper's Figure 1 display window clips valley extrema. For example, at `phi=0.15*pi`, the minimum calculated `L11/hbar` is approximately `-6.3983`; the third panel shows the full range.

Exact gap-closing points in Figure 2 are saved as `NaN`. Curves break at these points rather than connecting across undefined values. Circles in the left panel are the predefined subset with `Δ(phi/pi)=0.05`. The right panel shows every finite-lattice phase sample, including the particularly large mesh effects near the transition. All samples are in the CSV.

## Measured results and limitations

| Check | Result |
|---|---|
| Figure 1: maximum difference in OAM absolute values between the 60×60 lattice and Bloch route, internal units | `1.56e-13` |
| Figure 1: `max(abs(L12))/hbar`, in order of the four fluxes | `0.033659, 0.037890, 0.028559, <3e-17` |
| Figure 3: 240×240 plateau, negative conductivity in e/(2π) | `0.27825232309197084` |
| Figure 3: independent 60×60 lattice plateau | `0.278252323085247` |
| Figure 3: maximum difference between the two meshes across the entire Fermi scan | `0.0092452 e/(2π)` |

The Figure 3 plateau agrees with the paper's plotted value of approximately `0.28`. The full curve reproduces the plateau, decay on either side, and particle-hole symmetry. The metallic regions still show finite-k sampling differences; the plateau's high-precision agreement should not be generalized to the whole scan.

Transitions occur at `phi/pi=±0.1715401788, ±0.8284598212`. Peaks near these points cannot yet be claimed converged on a uniform finite mesh. For example, at `phi=0.15*pi`, the negative conductivities on the 60, 120, 240, and 480 meshes are `2.497735, 1.192520, 1.110162, 1.109901 e/(2π)`. At the `phi=0.5*pi` plateau, the 240 and 480 meshes agree to the displayed precision. This is deterministic numerical convergence evidence, not a statistical confidence interval.

The paper does not publish original numerical tables or figure-generation code, so this deliverable does not claim pixel-for-pixel or point-for-point replication of all author sampling. Visual inspection used pages 3 and 4 of the repository PDF and the original figures accompanying the public arXiv source. Original images are cached only in `.local/` and are not delivered as new calculation results.

## Visualization and software provenance

Calculations use Matplotlib, NumPy, and SciPy; versions are recorded in the manifest. The scientific visualization workflow references:

Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026).
*Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*.
[arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065), version v2 at the time of retrieval.

This reference identifies the visualization workflow; it does not certify the calculation formulas or compliance with publication requirements.
