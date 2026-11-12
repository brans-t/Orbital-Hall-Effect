# 04 · Formula verification report

> **Historical design record:** This document preserves the original design. Claims about the Eq. (16b) sign, phases, Pezo truncation, the periodic position operator, units, and some benchmarks were corrected during implementation. [05-implementation-notes.md](05-implementation-notes.md) and the executable tests take precedence.

This document records independent design-phase numerical checks. Their purpose was to identify a reference formula for implementation, not reproduce figures. All checks used `numpy` independently of the library, which was not yet implemented.

Original conclusions:

1. **Eq. (16a) agrees with the Eq. (2) definition** pointwise to about $10^{-15}$.
2. **Historical claim of an Eq. (16b) printing error:** a negative sum sign was proposed; the other sign differed from Eq. (16a) by $O(1)$.
3. **The diagonal needs its own closed form:** Eq. (16b) was used only for $n\ne n'$.
4. **Subtract diagonal position elements in the real-space oracle:** retaining them was reported to give $O(10)$ artifacts.
5. **Retain the Haldane next-nearest-neighbor real term $d_0$:** omitting it gives $L_{12}\equiv0$, contrary to Figure 1.
6. **Two-band models cannot test the Eq. (16b) sum:** use at least three bands.

---

## 1. Verification environment and models

**Haldane model**, paper Eq. (19):

$$\mathcal H = M\sum_i\xi_i c^\dagger_i c_i + t_1\sum_{\langle i,j\rangle}c^\dagger_ic_j + t_2\sum_{\langle\langle i,j\rangle\rangle}e^{i\nu_{ij}\phi}c^\dagger_ic_j$$

Parameters: $t_1=1$, $t_2=0.15$, $M=0.4$. In the two-band $\{A,B\}$ basis:

$$\mathcal H(\mathbf k)=\begin{pmatrix} d_0+d_z & t_1 f(\mathbf k)\\ t_1 f^*(\mathbf k) & d_0-d_z\end{pmatrix},\quad f=\sum_{i=1}^{3}e^{i\mathbf k\cdot\boldsymbol\delta_i} $$

where

$$d_0(\mathbf k)=2t_2\cos\phi\sum_{i=1}^{3}\cos(\mathbf k\cdot\mathbf b_i),\qquad d_z(\mathbf k)=M-2t_2\sin\phi\sum_{i=1}^{3}\sin(\mathbf k\cdot\mathbf b_i)$$

Here $\boldsymbol\delta_i$ are nearest-neighbor vectors and $\mathbf b_i$ next-nearest-neighbor vectors. **The $d_0$ term is essential for completeness**; see §5.

**Random multiband model**, to exercise the sum: four bands on a square lattice, random complex $t(\mathbf R)$ hoppings along six directions $+x,+y,+xy,-x+y,+x-y,2x$, and random Hermitian onsite terms. The sum of $\mathbf v_{nn}$ is not identically zero, and the intermediate-band sum is nonempty.

---

## 2. Eq. (16a) versus Eq. (2): passed

Compare projected Eq. (16a) with the direct Eq. (2) definition using $r_{nn'}=-iv_{nn'}/(\varepsilon_n-\varepsilon_{n'})$, on Haldane and random three-/four-band models along $\Gamma$-K and at random $\mathbf k$ points:

| Model | $\max|L_{16a}-L_{\rm eq2}|$ | Hermiticity residual |
|---|---|---|
| Haldane, four k points | $<10^{-8}$, with finite-difference $A_n$ | $<10^{-8}$ |
| Random three-band, two k points | $9\times10^{-16}$, $3\times10^{-16}$ | $<10^{-10}$ |
| Random four-band, two k points | $2.5\times10^{-15}$, $2.2\times10^{-15}$ | $<10^{-15}$ |

Conclusion: **Eq. (16a) is a reliable reference.** Its matrix-product implementation, theory §4.3, matches the pointwise definition and can be the main route.

The original checks shifted $\mathcal A_n$ by an arbitrary constant vector $\mathbf a$ and reported unchanged $L$, attributing this to cancellation of $\mathcal A_n$ terms when $v_{nn}$ is real. The implementation notes clarify the gauge interpretation.

---

## 3. Eq. (16b): historical sign-error claim

### 3.1 Observation

The original comparison implemented the expression described as printed, with a **+** sum sign, against Eq. (2):

| Model | $\max|L_{16b(+)}-L_{\rm eq2}|$ | $\max|L_{16b(-)}-L_{\rm eq2}|$ |
|---|---|---|
| Random three-band, k=(0.21,0.37) | $1.24$ | $1.8\times10^{-16}$ |
| Random three-band, k=(1.1,-0.8) | $0.83$ | $5.7\times10^{-17}$ |
| Random three-band, k=(2.7,0.15) | $0.93$ | $1.3\times10^{-16}$ |
| Random four-band, k=(0.23,0.41) | $1.51$ | $1.7\times10^{-8}$ |

Changing the sum sign to **negative** reduced residuals to machine precision. The report identified this as the correct Eq. (16b); the implementation notes explain the denominator-transcription error behind the finding.

### 3.2 Why two bands cannot distinguish the signs

In a two-band model, $m$ in $\sum_{m\ne n,n'}$ has only two values, already occupied by $n,n'$. The sum is empty, and Eq. (16b) reduces to:

$$L_{nn'}\big|_{\rm 2band}=\frac{i}{2}\cdot\frac{(\mathbf v_{nn}+\mathbf v_{n'n'})\times\mathbf v_{nn'}}{\varepsilon_n-\varepsilon_{n'}}$$

The signs are indistinguishable here. The original report used this to explain the paper's Haldane-only validation and motivate three-band tests. Its attribution to a printing error is superseded by the implementation notes.

### 3.3 Claimed effect on OHC

The historical report claimed that the sign error **reverses $\sigma^{\rm OH}$**, rather than merely changing its magnitude:

| Variant | $\sigma^{\rm OH}_{xy}$ at $\phi=\pi/4$ | $\phi=\pi/2$ |
|---|---|---|
| Eq. (2) definition | $-0.095$ | $-0.030$ |
| Eq. (16b), proposed corrected sign | $-0.171$ | $-0.338$ |
| Eq. (16b), described as printed | $+0.166$ | $+0.338$ |

Units: $e/2\pi$, on a $60\times60$ mesh. The report compared these values with the negative region of Figure 2 and proposed polarity as a sensitive CI check. The implementation notes supersede this two-band sign-flip claim.

> **Historical interpretation:** the table gives different Eq. (2) and proposed Eq. (16b) results at $\phi=\pi/2$: $-0.030$ versus $-0.338$. It attributed this to the empty two-band sum in Eq. (2) and strong contributions near a claimed touching at $\phi=\pi/2$, concluding that **two-band formulas underestimate OHC near degeneracy** and recommending full Eq. (16a)/(16b). These interpretations are retained as history, not current validation results.

---

## 4. Closed-form diagonal

Eq. (16b) was used for off-diagonal elements. The diagonal follows from the definition:

$$L_{nn}=2\sum_{m\ne n}\frac{\mathrm{Im}\big[v^x_{nm}(v^y_{nm})^*\big]}{\varepsilon_n-\varepsilon_m}$$

Compare each band of a random three-band model to the Eq. (2) diagonal. Of $+2\mathrm{Im}$, $-2\mathrm{Im}$, and $-\mathrm{Im}[\text{cross}]$, only $+2\mathrm{Im}[v^x(v^y)^*]$ agreed to $10^{-16}$; the others had the opposite sign. **Preserve this sign in the implementation.**

---

## 5. The Haldane $d_0$ pitfall

Omitting the real next-nearest-neighbor contribution $d_0(\mathbf k)$ while retaining only $d_z$ makes $\mathcal H$ traceless:

$$\mathrm{tr}\,\mathcal H=0\;\Longrightarrow\;\mathbf v_{11}+\mathbf v_{22}=0\;\Longrightarrow\;L_{12}\equiv0\ (\text{from the first term of 16b})$$

Numerical comparison:

| $k$ along $\Gamma$-K | With $d_0$: $|L_{12}|$ | Without $d_0$: $|L_{12}|$ | $\|\mathbf v_{11}+\mathbf v_{22}\|$, with/without |
|---|---|---|
| $0.2\,\Gamma K$ | $0.0064$ | $0$ | $1.15$ / $0$ |
| $0.4\,\Gamma K$ | $0.0448$ | $0$ | $1.74$ / $0$ |
| $0.6\,\Gamma K$ | $0.1204$ | $0$ | $1.54$ / $0$ |
| $0.8\,\Gamma K$ | $0.2011$ | $0$ | $0.79$ / $0$ |

The upper panel of Figure 1 has a $|L_{12}|$ peak near $0.04\hbar$, comparable in scale to the result with $d_0$. **Retain the $d_0$ term.**

The original report proposed a $\Gamma$-point check: $\varepsilon(\Gamma)=\pm3.0265$ with $d_0$ versus $\pm3.0$ without it, describing $d_z$ at $\Gamma$ as near $M+4t_2\cos\phi=\sqrt2-0.4$. These historical energy values are corrected in the implementation notes.

---

## 6. Two finite-lattice oracle pitfalls

Construct $L$ directly from Eq. (1) on an $N\times N$ periodic real-space lattice and compare to Bloch at allowed $\mathbf k=(m_1,m_2)/N$ points.

**Pitfall A: retaining diagonal position elements gives artifacts.** The real-space $\mathbf r$ diagonal contains absolute coordinates of order $O(N)$, not $\mathcal A_n$. Without subtraction, the checks reported:

| $N=24$ | $L_{11}$ | $\max|L_{12}|$ |
|---|---|---|
| Retain diagonal $\mathbf r$ | $-39.1$ artifact, varying with $N$ | Fluctuates from $10.5$ to $269$ |
| Subtract diagonal $\mathbf r$ | $0.0000$, correct at $\Gamma$ | Agrees with Bloch |

**Pitfall B: lattice Hermiticity.** Adding only half of $(i,j,R)$, such as `H[J,I]+=t` without `H[I,J]+=conj(t)`, can give Hermiticity error $1.0$. Include conjugate-transpose partners of directed hoppings.

After subtraction, the report found real-space/Bloch energy agreement, with $\varepsilon$ differences $<10^{-5}$ attributed to sorting. It also found $L_{11}$ equal to $0$ at $\Gamma$, as required by symmetry.

> **Original recommendation:** construct the oracle's $\mathbf r$ as `r_offdiag = r - diag(diag(r))`, document the reason, and provide no option to retain the diagonal. The implementation notes specify the correct basis.

---

## 7. Symmetry and topology constraints

| Check | Result |
|---|---|
| $L_{11}(k)=-L_{11}(-k)$ at $\phi=0$, where $\mathcal T$ is present | Agreement for $\pm0.01109$ and $\pm0.16389$ to $10^{-16}$ |
| $|L_{12}|\to0$ at $\phi=\pi/2$ | Order $10^{-17}$ |
| $L_{11}=L_{22}$ in a two-band model | Exactly equal, difference $0$ |
| Hermiticity $L^\dagger=L$ | $<10^{-15}$ |
| Chern transition at $\phi_c=\arcsin[M/(3\sqrt3t_2)]/\pi\approx0.1715\pi$ | $0\to-1$, described here as agreeing with the Figure 1 inset |

**Original FHS note:** sum logarithms of the Berry phase on each plaquette. The report claimed that $-\arg$ rather than $\log$ introduces integer multiples of $2\pi$. It obtained $-1$ over $0.20\pi$-$0.50\pi$ and $0$ over $0$-$0.17\pi$, consistent with $\phi_c=0.1715\pi$. The implementation notes correct the arg/log statement.

---

## 8. Quantitative comparison with the paper

At $\phi=\pi/4$, scan $k$ along $\Gamma\to K\to K'$ with $401$ points:

| Quantity | Original verification | Paper Figure 1 |
|---|---|---|
| $\max|L_{12}|$ along the path | $0.22\hbar$ | Peak near $0.04\hbar$ in the displayed range |
| $L_{11}$ range | $[-14.4, 0]$ | Approximately $\approx[-1,1]$ displayed |

**Original explanation:** $L$ becomes singular near $K$ as $\varepsilon_n-\varepsilon_m\to$ becomes small, so near-degenerate regions exceed the plotted range. The report suggested Figure 1 excludes part of the region near $K$ or clips its $y$ axis. It prioritized three-path agreement to $10^{-12}$ over amplitudes that depend on sample selection. Later reproduction notes address physical scale assumptions.

**Original OHC comparison:** on a $60\times60$ mesh at $E_F=0$ and $\phi=\pi/4$, the report obtained $\sigma^{\rm OH}\approx-0.17\,e/2\pi$, with sign and scale described as matching Figure 2 at that $\phi$. Checking the $E_F$ scan and Figure 3 plateau was deferred to the complete implementation and denser meshes.

---

## 9. Constraints on the architecture

Original design constraints from these checks:

1. Make `OAMVariant` pluggable to compare the proposed sign correction and truncated approximation.
2. Use matrix Eq. (16a) as the main route and Eq. (16b) as a cross-check with the **proposed corrected sign**.
3. Implement the closed-form diagonal separately from Eq. (16b).
4. Subtract the diagonal from $\mathbf r$ in `FiniteLattice`, with no opt-out.
5. Retain full $d_0(\mathbf k)$ in `HaldaneModel` and lock $\varepsilon(\Gamma)$ in tests.
6. **Include random models with at least three bands** to detect the Eq. (16b) sum-sign difference.
7. Assert the sign of $\sigma^{\rm OH}$ independently, negative for the proposed correction.
8. Flag near-degeneracies explicitly rather than silently truncating.

---

## Appendix: verification script location

The checks used one-off scripts not added to the repository. They lived in a temporary `/tmp` directory and expired with the session. This document records their methods, parameters, and reported results. The original plan was to turn those values into executable assertions in `tests/test_oam_consistency.py` during P1.
