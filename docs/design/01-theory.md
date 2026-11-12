# 01 · Theory and computational specification

> **Historical design record:** This document preserves the original design. Claims about the Eq. (16b) sign, phases, Pezo truncation, the periodic position operator, units, and some benchmarks were corrected during implementation. [05-implementation-notes.md](05-implementation-notes.md) and the executable tests take precedence.

This document translates the theory in Gautam & Satpathy, *Phys. Rev. B* **114**, 214403 (2026), referred to below as "the paper," into a mathematical specification suitable for implementation. It contains no implementation and addresses three questions: what to compute, the computable form of each formula, and the pitfalls to avoid.

All formulas were independently checked in the three-way numerical comparisons recorded in `docs/design/04-validation-report.md`. Where the printed expressions appeared numerically inconsistent, this historical document gives and labels a proposed correction.

---

## 1. Quantities to compute

| Quantity | Symbol | Definition | Output shape |
|---|---|---|---|
| Orbital angular momentum matrix elements | $L_{nn'}(\mathbf k)$ | Paper Eq. (2) | One complex $N_b\times N_b$ matrix per k point |
| Orbital Hall conductivity | $\sigma^{\rm OH}_{xy}$ | Paper Eq. (22) | Scalar, or the tensor's $xy$ component |
| Orbital Hall conductivity, Fermi-energy scan | $\sigma^{\rm OH}_{xy}(E_F)$ | Eq. (22) at $\mu=E_F$ | One-dimensional array |
| Finite-lattice reference OAM | $L^{\rm lat}_{ij}$ | Eq. (1), evaluated directly in real space | Large matrix for validation only |

The unit convention throughout the library is $\hbar = 1$, $m_e = 1$, and $e = 1$, with lengths in lattice units. The output layer restores physical quantities in units such as $\hbar$ and $\mu_B$ using conversion factors. Setting the repeated $m_e/\hbar$ factors to one keeps only structural coefficients and reduces omissions.

---

## 2. Core operators and definitions

### 2.1 Symmetrized OAM operator (Eq. 1)

$$\mathbf L = \frac{m_e}{2}\,(\mathbf r\times\mathbf v - \mathbf v\times\mathbf r)$$

**Why symmetrization is required:** $\mathbf r$ and $\mathbf v$ do not commute, so $\mathbf r\times\mathbf v$ alone is not Hermitian. The antisymmetric combination satisfies $\mathbf L^\dagger = \mathbf L$. This is an implementation constraint: a non-Hermitian $L$ matrix is a bug and should fail an assertion.

For the $z$ component with $\hbar=m_e=1$:

$$L_z = \tfrac12\big[(x v_y - y v_x) - (v_x y - v_y x)\big]$$

### 2.2 Gauge-covariant form (Eq. 2)

$$L_{nn'}(\mathbf k) = \frac{m_e}{2}\big\langle u_{n\mathbf k}\big|(\mathbf r - \mathcal A_{n\mathbf k})\times \mathbf v_{\mathbf k} - \mathbf v_{\mathbf k}\times(\mathbf r - \mathcal A_{n'\mathbf k})\big|u_{n'\mathbf k}\big\rangle$$

Here $\mathcal A_{n\mathbf k} = i\langle u_{n\mathbf k}|\partial_{\mathbf k}u_{n\mathbf k}\rangle$ is the Berry connection.

**Key point:** $\mathcal A_n$ is not an optional modification; it is the unique choice that removes dependence on fractional shifts of the origin within a unit cell. The wave-packet argument of Niu and coworkers provides this offset. Omitting it, by setting $\mathcal A=0$, makes $L$ depend on the Wannier-center convention.

> **Pitfall 1, numerically checked:** $\mathcal A_n$ and the diagonal position element $r_{nn}$ must be handled **as a pair**. Eq. (2) contains $(\mathbf r - \mathcal A_{n})$, so at the matrix-element level:
> $$r_{nn'} \;\longrightarrow\; \begin{cases} \mathcal A_{n\mathbf k}, & n=n' \quad(\text{diagonal, i.e. the gauge part of } r)\\[2pt] \dfrac{-i\,\hbar\, v_{nn'}}{\varepsilon_n-\varepsilon_{n'}}, & n\neq n' \end{cases}$$
> In the derivation, the diagonal $r$ element $r_{nn}$ cancels $\mathcal A_n$ exactly; see part (iii) of the paper's Supplemental Material. Simply setting the diagonal to zero without replacing it by $\mathcal A_n$ is valid only in a gauge where that diagonal is already zero. **The recommended route constructs $r_{nn'}$ directly in the band basis as above, without explicitly computing $\mathcal A$**. In the $\mathbf r\times\mathbf v - \mathbf v\times\mathbf r$ structure of Eq. (2), the $n=n'$ terms are claimed here to cancel when $v_{nn}$ is real and can be omitted. See §4.3 for the computable expression.

---

## 3. Structure of the Wannier-basis derivation (Eqs. 3-15)

The code need not repeat this derivation, but understanding it explains the final expression.

Bloch and Wannier states are transformed into one another in Eq. (3). The Wannier position expectations in Eqs. (8)-(9) show a Berry-connection offset of the Wannier center relative to the lattice site: $\langle n\mathbf R|\mathbf r|n\mathbf R\rangle = \mathbf R + \frac{\Omega_c}{(2\pi)^3}\int_{\rm BZ} d\mathbf k\,\mathcal A_{n\mathbf k}$. This is the concrete meaning of a position operator that is well defined in the Wannier basis but ill defined in the Bloch basis.

Expanding the two terms of Eq. (2) in the Wannier basis gives Eq. (5):

$$\langle u_{n\mathbf k}|(\mathbf r\times\mathbf v_{\mathbf k} - \mathbf v_{\mathbf k}\times\mathbf r)|u_{n'\mathbf k}\rangle = \sum_{\mathbf R} e^{i\mathbf k\cdot\mathbf R}\langle n0|\mathbf r\times\mathbf v|n'\mathbf R\rangle - e^{-i\mathbf k\cdot\mathbf R}\langle n\mathbf R|\mathbf v\times\mathbf r|n'0\rangle$$

**First simplification from symmetrization:** The $\mathbf R\times\langle\mathbf v\rangle$ terms in $\mathbf r\times\mathbf v$ and $\mathbf v\times\mathbf r$ cancel exactly, leaving combinations involving only $\mathbf r-\mathbf R$. This is the practical advantage of the symmetric form emphasized in the paper: explicit dependence on $\mathbf R$ disappears from individual Wannier matrix elements.

Substituting the velocity operator in Eq. (6), $v = -\frac{i}{\hbar}[\mathbf r,H]$, makes the first term, Eq. (7), a triple product of two position matrix elements around a Hamiltonian. The periodic-gauge $\delta$ function in Eq. (11) removes one momentum integral; removing the lattice sum then gives Eqs. (12)-(14).

### 3.1 Main result: Eq. (16)

The paper gives two equivalent expressions. **Eq. (16a)** is the projected form:

$$L_{nn'}(\mathbf k)=\frac{i m_e}{\hbar}\Big\langle \partial_{\mathbf k}u_{n\mathbf k}\Big|\Big(\mathcal H_{\mathbf k}-\frac{\varepsilon_{n\mathbf k}+\varepsilon_{n'\mathbf k}}{2}\Big)\Big|\partial_{\mathbf k}u_{n'\mathbf k}\Big\rangle + \frac{m_e}{2}\Big[\mathcal A_{nn'\mathbf k}\times(\mathbf v_{nn\mathbf k}+\mathbf v_{n'n'\mathbf k})-(\mathcal A_{n\mathbf k}+\mathcal A_{n'\mathbf k})\times\mathbf v_{nn'\mathbf k}\Big]$$

The cross product acts in two-dimensional $\mathbf k$ space: $(\mathbf a\times\mathbf b)_z = a_x b_y - a_y b_x$. Also, $\mathcal A_{nn'\mathbf k}\equiv i\langle u_{n\mathbf k}|\partial_{\mathbf k}u_{n'\mathbf k}\rangle$, including when $n\ne n'$.

**Eq. (16b)** is the spectral-sum form:

$$L_{nn'}(\mathbf k)=\frac{i m_e\hbar}{2}\Bigg[\sum_{m\neq n,n'}\Big(\frac{1}{\varepsilon_{n\mathbf k}-\varepsilon_{m\mathbf k}}+\frac{1}{\varepsilon_{n'\mathbf k}-\varepsilon_{m\mathbf k}}\Big)\mathbf v_{nm\mathbf k}\times\mathbf v_{mn'\mathbf k}+\frac{(\mathbf v_{nn\mathbf k}+\mathbf v_{n'n'\mathbf k})\times\mathbf v_{nn'\mathbf k}}{\varepsilon_{n\mathbf k}-\varepsilon_{n'\mathbf k}}\Bigg]$$

> **Pitfall 2, historical claim of a printing error:** In Eq. (16b), this document proposes a **negative** sign before the sum:
> $$\sum_{m\neq n,n'}\Big(\tfrac{1}{\varepsilon_n-\varepsilon_m}+\tfrac{1}{\varepsilon_{n'}-\varepsilon_m}\Big)(-\mathbf v_{nm}\times\mathbf v_{mn'})$$
> Basis for the historical claim: expand $|\partial_{\mathbf k}u\rangle$ with Eq. (17) and compare the algebraic rearrangement of Eq. (16a) with the direct Eq. (2) definition on a random three-band model. The expression described here as printed differs by $O(1)$, while the proposed correction gives a $10^{-16}$ residual. In a two-band model the sum is empty, so **two-band tests alone cannot detect this difference**. This motivates mandatory regression tests with at least three bands.
>
> Eq. (16b) is used here only for $n\ne n'$; diagonal elements $L_{nn}$ are implemented separately in §4.2.

### 3.2 Symmetry constraints (Eq. 18)

Under time reversal $\mathcal T$ and spatial inversion $\mathcal I$:

$$L_{nn'}(-\mathbf k)=L_{nn'}(\mathbf k)\quad(\mathcal I\ \text{present}),\qquad L_{nn'}(-\mathbf k)=-L^{*}_{nn'}(\mathbf k)\quad(\mathcal T\ \text{present})$$

These are inexpensive checks that require no analytic solution. Numerically, the Haldane model at $\phi=0$, where $\mathcal T$ is present, satisfies $L_{11}(k)=-L_{11}(-k)$ to $10^{-15}$.

---

## 4. Computational specification

### 4.1 Input abstraction for a general tight-binding model

At a given $\mathbf k$, every model must provide two quantities:

$$\mathcal H(\mathbf k), \qquad \partial_{k_j}\mathcal H(\mathbf k),\quad j=x,y$$

For Wannier90 `_hr.dat`, both follow analytically:

$$\mathcal H(\mathbf k)=\sum_{\mathbf R} H(\mathbf R)\,e^{i\mathbf k\cdot\mathbf R},\qquad \partial_{k_j}\mathcal H(\mathbf k)=\sum_{\mathbf R} iR_j\,H(\mathbf R)\,e^{i\mathbf k\cdot\mathbf R}$$

**Do not use finite differences for $\partial_{k_j}\mathcal H$ when an analytic derivative is available.** With known $H(\mathbf R)$, the analytic derivative is exact and inexpensive. Finite differences amplify errors at small gaps, precisely where OAM is singular. They are a fallback for a black-box $\mathcal H(\mathbf k)$ interface only, using central differences, a sufficiently small step, and an explicit warning.

### 4.2 Gauge-invariant eigenbasis form (recommended main route)

In the eigenbasis $\{|u_{n\mathbf k}\rangle\}$ of $\mathcal H(\mathbf k)$, obtained with `numpy.linalg.eigh`, define

$$v^{j}_{nm}=\langle u_n|\partial_{k_j}\mathcal H|u_m\rangle,\qquad j\in\{x,y\}$$

This band velocity is **not corrected by the Berry connection**. Its diagonal $v^j_{nn}=\partial_{k_j}\varepsilon_n/\partial k_j$ is real, a key property.

**Diagonal OAM**, $n=n'$, has a closed form derived directly from the definition. Since $\mathbf v_{mn}=\mathbf v_{nm}^*$, $\mathrm{Im}[\mathbf v_{nm}\times\mathbf v_{mn}]=2\,\mathrm{Im}[v^x_{nm}(v^y_{nm})^*]$:

$$L_{nn}=-i\sum_{m\ne n}\frac{\mathbf v_{nm}\times\mathbf v_{mn}}{\varepsilon_n-\varepsilon_m}=2\sum_{m\ne n}\frac{\mathrm{Im}\big[v^x_{nm}\,(v^y_{nm})^*\big]}{\varepsilon_n-\varepsilon_m}$$

This historical discussion describes the second equality as having a **negative-sign form**. Writing $\mathrm{Im}[v^x v^y]$ instead can introduce an overall sign error because $\mathbf v_{mn}=\mathbf v_{nm}^*$ already supplies the conjugation. The expression above was checked pointwise against Eq. (2) to $10^{-16}$; writing $-\mathrm{Im}[v^x_{nm}v^y_{mn}]$ gives the opposite sign. With $\mathcal T$ present, $L_{nn}$ is real. Even when $\mathcal T$ is broken, its physical value is real because Hermitian $L_z$ has real diagonal entries.

**Off-diagonal OAM**, $n\ne n'$, uses the proposed corrected Eq. (16b) from §3.1:

$$L_{nn'}=\frac{i}{2}\Bigg[\frac{(\mathbf v_{nn}+\mathbf v_{n'n'})\times\mathbf v_{nn'}}{\varepsilon_n-\varepsilon_{n'}}-\sum_{m\ne n,n'}\Big(\frac{1}{\varepsilon_n-\varepsilon_m}+\frac{1}{\varepsilon_{n'}-\varepsilon_m}\Big)\mathbf v_{nm}\times\mathbf v_{mn'}\Bigg]$$

**Why the reality of $v_{nn}$ matters in the original argument:** It makes $\mathbf v_{nn}\times\mathbf v_{nn}=0$ zero. This historical argument therefore claims that $\mathcal A_n\times\mathbf v_{nn}$ in the $(\mathbf r-\mathcal A_n)\times\mathbf v$ term of Eq. (2) vanishes and that $r_{nn}$ disappears from the final expression. **The intended result is an algorithm requiring neither an explicit Berry connection nor gauge fixing.** This was the proposed basis for converting the paper's gauge-covariant formula into a gauge-invariant numerical route; see the implementation corrections for the revised explanation.

### 4.3 Matrix-product implementation without elementwise loops

Triple elementwise loops become a bottleneck for $N_b>20$. Both formulas above can be expressed as matrix products.

Define the derivative coefficient matrix $C_j$ through the off-diagonal part of $|\partial_{k_j}u_n\rangle = \sum_m (C_j)_{mn}|u_m\rangle$:

$$(C_j)_{mn}=\frac{v^j_{mn}}{\varepsilon_n-\varepsilon_m}\quad(m\ne n),\qquad (C_j)_{nn}=0$$

Then $\langle\partial_j u_n|\mathcal H|\partial_l u_{n'}\rangle = (C_j^\dagger \mathcal H C_l)_{nn'}$ and $\langle\partial_j u_n|\partial_l u_{n'}\rangle = (C_j^\dagger C_l)_{nn'}$. The first term, the projected part of Eq. (16a), is

$$T^{(1)} = i\Big[(C_x^\dagger \mathcal H C_y - \bar\varepsilon\, C_x^\dagger C_y) - (x\leftrightarrow y)\Big],\qquad \bar\varepsilon_{nn'}=\tfrac12(\varepsilon_n+\varepsilon_{n'})$$

The second term is $T^{(2)} = \frac12\big[\mathcal A_{nn'}\times(\mathbf v_{nn}+\mathbf v_{n'n'})\big]$, with $\mathcal A_{nn'}=-i\,v_{nn'}/(\varepsilon_{n'}-\varepsilon_n)$ for $n\ne n'$. These are $O(N_b^3)$ matrix products with no elementwise loops. This form was checked against the pointwise definition with a maximum difference of $2.5\times10^{-15}$.

### 4.4 Finite-lattice reference, for validation only

The paper's second method, Eqs. (20)-(21), constructs real-space $\mathcal H$ on a finite $N\times N$ lattice with periodic boundary conditions. Eigenvectors are coefficients in the $|i\rangle$ basis. The velocity elements are

$$\langle i|\mathbf v|j\rangle = \frac{1}{i\hbar}(\mathbf r_i-\mathbf r_j)\,t_{ij}$$

Only off-diagonal position elements are retained:

$$\langle\psi_n|\mathbf r|\psi_m\rangle = \frac{-i\hbar\langle\psi_n|\mathbf v|\psi_m\rangle}{\varepsilon_n-\varepsilon_m}\quad(n\ne m)$$

Then assemble $L$ directly using Eq. (1). **This step requires care:**

> **Pitfall 3, historical position-diagonal prescription:** The real-space $\mathbf r$ matrix is described here as requiring **diagonal subtraction before forming $\mathbf r\times\mathbf v$**. Retaining the diagonal gives spurious $L_{11}\sim O(10)$ and large $|L_{12}|$, unlike the Bloch result. After subtraction, the two were reported to agree pointwise. The original explanation is that $\mathbf r_{ii}$ contains absolute lattice coordinates that grow with system size, while $\mathcal A_i$ is the corresponding gauge quantity; these must not be mixed. The implementation notes clarify the required basis for this subtraction.
>
> A second pitfall: the Haldane next-nearest-neighbor hopping has the real contribution $d_0(\mathbf k)=2t_2\cos\phi\sum_i\cos(\mathbf k\cdot\mathbf b_i)$. **Omitting it gives $\mathrm{tr}\,\mathcal H=0$ and $\mathbf v_{11}+\mathbf v_{22}=0$, making the off-diagonal term in Eq. (16b) identically zero**, contrary to the paper's Figure 1 value of $|L_{12}|\approx0.04$. A two-band result with $L_{12}\equiv0$ should first prompt a check for missing Hamiltonian terms.

---

## 5. Orbital Hall conductivity (Eq. 22)

$$\sigma^{\rm OH}_{xy}=\frac{2e\hbar}{A}\,\mathrm{Im}\sum_{n}^{\rm occ}\sum_{m\ne n}\frac{\langle\psi_{n}|J_x|\psi_{m}\rangle\langle\psi_{m}|v_y|\psi_{n}\rangle}{(\varepsilon_n-\varepsilon_m)^2}$$

The orbital-current operator is

$$J_x = \tfrac12(v_x L_z + L_z v_x)$$

**Why this also requires symmetrization:** As in §2.1, $v_x$ and $L_z$ do not commute. The paper explicitly uses the symmetric form.

Implementation: insert $L$ and $\mathbf v$ from §4.2 directly into the band-basis expression. $A$ is the unit-cell area in 2D or volume in 3D. Occupation is $\varepsilon_n<\mu$. For an insulator, $\mu$ may lie anywhere in the gap; for a metal, use the specified $E_F$.

> The value of $\sigma^{\rm OH}$ is sensitive to the sign of $L$. This historical note claims that leaving pitfall 2 uncorrected flips $\sigma^{\rm OH}$ from $-0.17$ to $+0.17$ in units of $e/2\pi$. It therefore proposes **OHC polarity as a sensitive $L$ sign-error test** for CI. The implementation notes supersede the two-band polarity claim.

---

## 6. Validation benchmarks from the paper

The paper uses the Haldane model, Eq. (19), with $t_1=1$, $t_2=0.15$, and $M=0.4$, scanning $\phi$ over $0\to0.5\pi$. The original quantitative targets were:

| Quantity | Paper location | Original target |
|---|---|---|
| Chern-number transition | Figure 1 inset | $0\to-1$ at $\phi_c=\arcsin[M/(3\sqrt3 t_2)]/\pi\approx0.1715\pi$ |
| $L_{12}(\mathbf k)$ along $\Gamma$-K | Figure 1, top | Peak near $0.04\hbar$; zero at $\phi=0$ |
| $L_{11}(\mathbf k)$ | Figure 1, bottom | At $\phi=0$, $L_{11}(k)=-L_{11}(-k)$; magnitude $\le1\hbar$ |
| $\sigma^{\rm OH}_{xy}(\phi)$ | Figure 2 | Diverges at the transition; varies with $\phi$ in the trivial region |
| $\sigma^{\rm OH}_{xy}(E_F)$ | Figure 3 | Plateau near $E_F=0$ at $\phi=\pi/2$, with magnitude about $0.28\,e/2\pi$ |

> **Interpretation:** The $|L_{12}|$ peak in Figure 1 and the absolute amplitudes in Figures 2-3 depend on $k$-mesh density, $E_F$, and the selected $\phi$. Acceptance is based on **three independent paths agreeing to $10^{-12}$ together with symmetry and topology constraints**, rather than pixel-level figure matching. Paper values serve as scale and shape checks.

---

## 7. Scope and limitations

The paper states that the formalism applies to **metals and insulators**. For an insulator, $\mu$ lies in the gap and $\sigma^{\rm OH}$ is a nonsingular Fermi-sea integral. Metals have partially occupied bands, and $(\varepsilon_n-\varepsilon_m)^{-2}$ is singular at band crossings, explaining the transition divergence in Figure 2. The library should detect near-degeneracy and expose a warning threshold rather than silently returning large values.

The spectral Eq. (16b) costs $O(N_b^3)$ per k point at large $N_b$, and the projected Eq. (16a) likewise costs $O(N_b^3)$. Both are practical for $N_b\lesssim 100$; larger systems require $k$-point parallelism and possibly a window near the occupied bands.
