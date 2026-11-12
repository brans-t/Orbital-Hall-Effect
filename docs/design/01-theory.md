# 01 · 理论梳理与可计算规范

> **历史设计记录：** 本文保留原始设计。关于 Eq. (16b) 符号、相位、Pezo 截断、周期位置算符、单位及部分基准值的结论已在实现阶段校订，以 [05-implementation-notes.md](05-implementation-notes.md) 和可执行测试为准。

本文件把 Gautam & Satpathy, *Phys. Rev. B* **114**, 214403 (2026)（下称"论文"）中的理论，翻译成一份可以直接照着写代码的数学规范。不含实现，只回答三个问题：要算什么、每一条公式的可计算形式是什么、实现时必须避开哪些陷阱。

文中所有公式都在 `docs/design/04-validation-report.md` 记录的三方数值对照中被独立验证过。凡论文印刷形式与数值不自洽之处，本文给出修正后的形式并标注。

---

## 1. 要计算的物理量

| 量 | 符号 | 定义 | 输出形状 |
|---|---|---|---|
| 轨道角动量矩阵元 | $L_{nn'}(\mathbf k)$ | 论文 Eq. (2) | $N_b\times N_b$ 复矩阵，每 k 点一个 |
| 轨道霍尔电导 | $\sigma^{\rm OH}_{xy}$ | 论文 Eq. (22) | 标量（或张量的 $xy$ 分量） |
| 轨道霍尔电导（费米能扫描） | $\sigma^{\rm OH}_{xy}(E_F)$ | Eq. (22) 在 $\mu=E_F$ 处 | 一维数组 |
| 有限格点基准 OAM | $L^{\rm lat}_{ij}$ | Eq. (1) 在实空间直接求 | 大矩阵，仅用于验证 |

单位约定贯穿全库：$\hbar = 1$、$m_e = 1$、$e = 1$，长度以晶格常数为单位。所有物理量在输出层再乘以单位换算因子还原为 $\hbar$ 和 $\mu_B$ 等。这样做是因为论文中 $m_e/\hbar$ 因子反复出现且极易漏掉，统一取 1 后公式只保留结构性系数。

---

## 2. 核心算符与定义式

### 2.1 对称化 OAM 算符（Eq. 1）

$$\mathbf L = \frac{m_e}{2}\,(\mathbf r\times\mathbf v - \mathbf v\times\mathbf r)$$

**为什么必须对称化**：$\mathbf r$ 与 $\mathbf v$ 不对易，$\mathbf r\times\mathbf v$ 单独使用不是厄米的。取反对称组合后 $\mathbf L^\dagger = \mathbf L$。这一条在实现上是硬约束——任何非厄米的 $L$ 矩阵都是 bug，应作为断言。

在 $\hbar=m_e=1$ 下，只关心 $z$ 分量时：

$$L_z = \tfrac12\big[(x v_y - y v_x) - (v_x y - v_y x)\big]$$

### 2.2 规范协变形式（Eq. 2）

$$L_{nn'}(\mathbf k) = \frac{m_e}{2}\big\langle u_{n\mathbf k}\big|(\mathbf r - \mathcal A_{n\mathbf k})\times \mathbf v_{\mathbf k} - \mathbf v_{\mathbf k}\times(\mathbf r - \mathcal A_{n'\mathbf k})\big|u_{n'\mathbf k}\big\rangle$$

其中 $\mathcal A_{n\mathbf k} = i\langle u_{n\mathbf k}|\partial_{\mathbf k}u_{n\mathbf k}\rangle$ 是 Berry 联络。

**关键点**：$\mathcal A_n$ 在这里不是可选修饰，而是使结果与"单位胞内原点平移"无关的唯一选择。Niu 等人的波包论证给出这个偏移。实现时如果忽略它（即令 $\mathcal A=0$），得到的 $L$ 会依赖 Wannier 中心的选取。

> **陷阱 1（已在数值中确认）**：$\mathcal A_n$ 与对角位置矩阵元 $r_{nn}$ 必须**成对**处理。Eq. (2) 中出现的组合是 $(\mathbf r - \mathcal A_{n})$，因此矩阵元层面成立的是
> $$r_{nn'} \;\longrightarrow\; \begin{cases} \mathcal A_{n\mathbf k}, & n=n' \quad(\text{对角，即 } r \text{ 的规范部分})\\[2pt] \dfrac{-i\,\hbar\, v_{nn'}}{\varepsilon_n-\varepsilon_{n'}}, & n\neq n' \end{cases}$$
> 推导中 $r$ 的对角元 $r_{nn}$ 与 $\mathcal A_n$ 精确相消（论文 SM 第 (iii) 部分）。实现时若"把对角元置零"而不替换成 $\mathcal A_n$，只有在对角元本就为零的规范下才正确。**推荐做法是在能带基下直接令 $r_{nn'}$ 为上式，不需要显式算 $\mathcal A$**——$n=n'$ 项在 Eq. (2) 的 $\mathbf r\times\mathbf v - \mathbf v\times\mathbf r$ 结构中最终 $v_{nn}$ 为实数时自动消去，可安全省略。详见 §4.3 的可计算形式。

---

## 3. Wannier 基下的推导结构（Eqs. 3–15）

这段推导本身不必在代码中复现，但理解它才能明白最终公式为何长这样。

Bloch 态与 Wannier 态互换（Eq. 3）；Wannier 态的位置期望 Eq. (8)(9) 表明 Wannier 中心相对格点偏移了 Berry 联络：$\langle n\mathbf R|\mathbf r|n\mathbf R\rangle = \mathbf R + \frac{\Omega_c}{(2\pi)^3}\int_{\rm BZ} d\mathbf k\,\mathcal A_{n\mathbf k}$。这正是"位置算符在 Wannier 基下良定义、在 Bloch 基下不良定义"的具体含义。

把 Eq. (2) 的两项分别在 Wannier 基展开，得 Eq. (5)：

$$\langle u_{n\mathbf k}|(\mathbf r\times\mathbf v_{\mathbf k} - \mathbf v_{\mathbf k}\times\mathbf r)|u_{n'\mathbf k}\rangle = \sum_{\mathbf R} e^{i\mathbf k\cdot\mathbf R}\langle n0|\mathbf r\times\mathbf v|n'\mathbf R\rangle - e^{-i\mathbf k\cdot\mathbf R}\langle n\mathbf R|\mathbf v\times\mathbf r|n'0\rangle$$

**对称化带来的第一处化简**：$\mathbf r\times\mathbf v$ 与 $\mathbf v\times\mathbf r$ 中各自出现的"$\mathbf R\times\langle\mathbf v\rangle$"项精确对消，剩下只含 $\mathbf r-\mathbf R$ 的组合。这就是论文强调"用对称形式是一个优势"的实际含义——它消掉了单个 Wannier 矩阵元中对 $\mathbf R$ 的显式依赖。

速度算符 Eq. (6) $v = -\frac{i}{\hbar}[\mathbf r,H]$，代入后第一项 Eq. (7) 是两个位置矩阵元夹一个 Hamiltonian 的三重积。利用周期规范下的 $\delta$ 函数 Eq. (11) 消掉一重动量积分，再消掉格点和，得 Eq. (12)(13)(14)。

### 3.1 主结果：Eq. (16)

论文最终给出两条等价形式。**Eq. (16a)** 是"投影形式"：

$$L_{nn'}(\mathbf k)=\frac{i m_e}{\hbar}\Big\langle \partial_{\mathbf k}u_{n\mathbf k}\Big|\Big(\mathcal H_{\mathbf k}-\frac{\varepsilon_{n\mathbf k}+\varepsilon_{n'\mathbf k}}{2}\Big)\Big|\partial_{\mathbf k}u_{n'\mathbf k}\Big\rangle + \frac{m_e}{2}\Big[\mathcal A_{nn'\mathbf k}\times(\mathbf v_{nn\mathbf k}+\mathbf v_{n'n'\mathbf k})-(\mathcal A_{n\mathbf k}+\mathcal A_{n'\mathbf k})\times\mathbf v_{nn'\mathbf k}\Big]$$

其中叉积作用在二维 $\mathbf k$ 空间里：$(\mathbf a\times\mathbf b)_z = a_x b_y - a_y b_x$，且 $\mathcal A_{nn'\mathbf k}\equiv i\langle u_{n\mathbf k}|\partial_{\mathbf k}u_{n'\mathbf k}\rangle$（$n\ne n'$ 时也定义）。

**Eq. (16b)** 是"求和形式"：

$$L_{nn'}(\mathbf k)=\frac{i m_e\hbar}{2}\Bigg[\sum_{m\neq n,n'}\Big(\frac{1}{\varepsilon_{n\mathbf k}-\varepsilon_{m\mathbf k}}+\frac{1}{\varepsilon_{n'\mathbf k}-\varepsilon_{m\mathbf k}}\Big)\mathbf v_{nm\mathbf k}\times\mathbf v_{mn'\mathbf k}+\frac{(\mathbf v_{nn\mathbf k}+\mathbf v_{n'n'\mathbf k})\times\mathbf v_{nn'\mathbf k}}{\varepsilon_{n\mathbf k}-\varepsilon_{n'\mathbf k}}\Bigg]$$

> **陷阱 2（已在数值中确认，论文印刷有误）**：在 Eq. (16b) 中，求和项前的符号应为**负**。即
> $$\sum_{m\neq n,n'}\Big(\tfrac{1}{\varepsilon_n-\varepsilon_m}+\tfrac{1}{\varepsilon_{n'}-\varepsilon_m}\Big)(-\mathbf v_{nm}\times\mathbf v_{mn'})$$
> 依据：把它作为 Eq. (16a) 的代数变形（用 Eq. 17 展开 $|\partial_{\mathbf k}u\rangle$），在 3 带随机模型上与 Eq. (2) 的直接定义对照，印刷符号偏差达 $O(1)$，改正后残差 $10^{-16}$；在 2 带模型上两项无法区分（求和无项），所以**只用 2 带模型做测试无法发现这个错误**。这是本库必须内置 3 带以上回归测试的直接原因。
>
> 同时 Eq. (16b) 只覆盖了 $n\ne n'$ 的情形；对角元 $L_{nn}$ 需另写（见 §4.2）。

### 3.2 对称性约束（Eq. 18）

时间反演 $\mathcal T$ 与空间反演 $\mathcal I$ 下：

$$L_{nn'}(-\mathbf k)=L_{nn'}(\mathbf k)\quad(\mathcal I\ \text{存在}),\qquad L_{nn'}(-\mathbf k)=-L^{*}_{nn'}(\mathbf k)\quad(\mathcal T\ \text{存在})$$

这两条是最便宜的测试：任何实现都能在没有解析解的情况下用它们自查。数值中已确认：Haldane 模型 $\phi=0$（$\mathcal T$ 存在）时 $L_{11}(k)=-L_{11}(-k)$ 精确成立到 $10^{-15}$。

---

## 4. 可计算规范（实现者照此写）

### 4.1 通用紧束缚模型的输入抽象

任何模型归结为在给定 $\mathbf k$ 下能给出两样东西：

$$\mathcal H(\mathbf k), \qquad \partial_{k_j}\mathcal H(\mathbf k),\quad j=x,y$$

对 Wannier90 的 `_hr.dat`，两者都可解析得到：

$$\mathcal H(\mathbf k)=\sum_{\mathbf R} H(\mathbf R)\,e^{i\mathbf k\cdot\mathbf R},\qquad \partial_{k_j}\mathcal H(\mathbf k)=\sum_{\mathbf R} iR_j\,H(\mathbf R)\,e^{i\mathbf k\cdot\mathbf R}$$

**不要用有限差分求 $\partial_{k_j}\mathcal H$**。$H(\mathbf R)$ 已知时解析导数精确且零成本；有限差分会在小能隙处放大误差，而 OAM 恰恰在小能隙处是奇异的。有限差分只在"只有黑盒 $\mathcal H(\mathbf k)$"的接口下作为兜底，并且必须用中心差分 + 足够小的步长，同时给出警告。

### 4.2 规范不变的本征基形式（推荐主路径）

在 $\mathcal H(\mathbf k)$ 的本征基 $\{|u_{n\mathbf k}\rangle\}$ 下（`numpy.linalg.eigh`），令

$$v^{j}_{nm}=\langle u_n|\partial_{k_j}\mathcal H|u_m\rangle,\qquad j\in\{x,y\}$$

注意这是**未经 Berry 联络修正**的"能带速度"；对角元 $v^j_{nn}=\partial_{k_j}\varepsilon_n/\partial k_j$ 是实数，这是关键性质。

**对角 OAM**（$n=n'$）有一个闭式，直接可由定义推得（$\mathbf v_{mn}=\mathbf v_{nm}^*$，故 $\mathrm{Im}[\mathbf v_{nm}\times\mathbf v_{mn}]=2\,\mathrm{Im}[v^x_{nm}(v^y_{nm})^*]$）：

$$L_{nn}=-i\sum_{m\ne n}\frac{\mathbf v_{nm}\times\mathbf v_{mn}}{\varepsilon_n-\varepsilon_m}=2\sum_{m\ne n}\frac{\mathrm{Im}\big[v^x_{nm}\,(v^y_{nm})^*\big]}{\varepsilon_n-\varepsilon_m}$$

注意第二个等号的**符号为负号形式**：写成 $\mathrm{Im}[v^x v^y]$ 时会差一个整体负号，因为 $\mathbf v_{mn}=\mathbf v_{nm}^*$ 已把复共轭引入。实测两种写法与 Eq. (2) 逐点吻合到 $10^{-16}$（上式），而写成 $-\mathrm{Im}[v^x_{nm}v^y_{mn}]$ 会得到整体反号的结果。当 $\mathcal T$ 存在时 $L_{nn}$ 为纯实数；$\mathcal T$ 破缺时一般仍取实部为物理值（$L_z$ 厄米 ⇒ 对角元必实）。

**非对角 OAM**（$n\ne n'$）：用 §3.1 的修正版 Eq. (16b)，

$$L_{nn'}=\frac{i}{2}\Bigg[\frac{(\mathbf v_{nn}+\mathbf v_{n'n'})\times\mathbf v_{nn'}}{\varepsilon_n-\varepsilon_{n'}}-\sum_{m\ne n,n'}\Big(\frac{1}{\varepsilon_n-\varepsilon_m}+\frac{1}{\varepsilon_{n'}-\varepsilon_m}\Big)\mathbf v_{nm}\times\mathbf v_{mn'}\Bigg]$$

**为什么 $v_{nn}$ 是实数很重要**：它使 $\mathbf v_{nn}\times\mathbf v_{nn}=0$，因此 Eq. (2) 中 $(\mathbf r-\mathcal A_n)\times\mathbf v$ 里的 $\mathcal A_n\times\mathbf v_{nn}$ 项为零；对角位置算符 $r_{nn}$ 因此在最终表达式中完全不出现。**这意味着算法不需要计算 Berry 联络，也不需要固定规范**。这是把论文"规范协变"变成"规范不变数值算法"的关键一步，也是本库主路径的立足点。

### 4.3 矩阵乘法实现（避免逐元素循环）

逐元素三重循环在 $N_b>20$ 时会成为瓶颈。上面两条公式都能写成矩阵乘：

定义"导数系数矩阵" $C_j$，满足 $|\partial_{k_j}u_n\rangle = \sum_m (C_j)_{mn}|u_m\rangle$ 的非对角部分：

$$(C_j)_{mn}=\frac{v^j_{mn}}{\varepsilon_n-\varepsilon_m}\quad(m\ne n),\qquad (C_j)_{nn}=0$$

于是 $\langle\partial_j u_n|\mathcal H|\partial_l u_{n'}\rangle = (C_j^\dagger \mathcal H C_l)_{nn'}$，$\langle\partial_j u_n|\partial_l u_{n'}\rangle = (C_j^\dagger C_l)_{nn'}$。第一项（Eq. 16a 的投影部分）为

$$T^{(1)} = i\Big[(C_x^\dagger \mathcal H C_y - \bar\varepsilon\, C_x^\dagger C_y) - (x\leftrightarrow y)\Big],\qquad \bar\varepsilon_{nn'}=\tfrac12(\varepsilon_n+\varepsilon_{n'})$$

第二项为 $T^{(2)} = \frac12\big[\mathcal A_{nn'}\times(\mathbf v_{nn}+\mathbf v_{n'n'})\big]$，其中 $\mathcal A_{nn'}=-i\,v_{nn'}/(\varepsilon_{n'}-\varepsilon_n)$（$n\ne n'$）。全部是 $O(N_b^3)$ 的矩阵乘，无循环。此形式已与逐点定义对照验证，最大偏差 $2.5\times10^{-15}$。

### 4.4 有限格点基准（仅验证用）

论文的第二种方法（Eqs. 20–21）在 $N\times N$ 有限格点 + PBC 上构造实空间 $\mathcal H$，本征矢为 $|i\rangle$ 基下的系数。速度矩阵元

$$\langle i|\mathbf v|j\rangle = \frac{1}{i\hbar}(\mathbf r_i-\mathbf r_j)\,t_{ij}$$

位置矩阵元只取非对角：

$$\langle\psi_n|\mathbf r|\psi_m\rangle = \frac{-i\hbar\langle\psi_n|\mathbf v|\psi_m\rangle}{\varepsilon_n-\varepsilon_m}\quad(n\ne m)$$

然后直接按 Eq. (1) 组装 $L$。**这一步必须小心**：

> **陷阱 3（已在数值中确认）**：实空间版的 $\mathbf r$ 矩阵必须**先扣除对角元再参与 $\mathbf r\times\mathbf v$ 的乘积**。若保留对角元，有限格点给出 $L_{11}\sim O(10)$ 的伪值且 $|L_{12}|$ 巨大，与 Bloch 结果彻底不符；扣除后两者逐点一致。原因是实空间对角元 $\mathbf r_{ii}$ 是格点绝对坐标（随体系尺寸增长），而 $\mathcal A_i$ 才是它的规范对应物，二者不能混用。
>
> 第二条陷阱：Haldane 模型的次近邻跃迁有实部 $d_0(\mathbf k)=2t_2\cos\phi\sum_i\cos(\mathbf k\cdot\mathbf b_i)$。**漏掉它会使 $\mathrm{tr}\,\mathcal H=0$，从而 $\mathbf v_{11}+\mathbf v_{22}=0$，使 Eq. (16b) 的非对角项恒为零**——正好与论文 Fig. 1 的 $|L_{12}|\approx0.04$ 矛盾。任何"二带模型 $L_{12}\equiv0$"的结果都应先怀疑模型漏项，而不是公式。

---

## 5. 轨道霍尔电导（Eq. 22）

$$\sigma^{\rm OH}_{xy}=\frac{2e\hbar}{A}\,\mathrm{Im}\sum_{n}^{\rm occ}\sum_{m\ne n}\frac{\langle\psi_{n}|J_x|\psi_{m}\rangle\langle\psi_{m}|v_y|\psi_{n}\rangle}{(\varepsilon_n-\varepsilon_m)^2}$$

轨道流算符

$$J_x = \tfrac12(v_x L_z + L_z v_x)$$

**为什么也要对称化**：与 §2.1 同理，$v_x$ 与 $L_z$ 不对易。论文明确写出对称形式。

实现要点：此式在能带基下直接代入 §4.2 得到的 $L$、$\mathbf v$ 即可；$A$ 是单位胞面积（二维）或体积（三维）。占据判据用 $\varepsilon_n<\mu$；对绝缘体 $\mu$ 可放在带隙内任意处，对金属则取给定的 $E_F$。

> 数值实现中，$\sigma^{\rm OH}$ 的值强烈依赖 $L$ 的正确符号。§3.1 陷阱 2 若未修正，$\sigma^{\rm OH}$ 的符号会整个翻转（实测：修正后 $-0.17$，印刷符号 $+0.17$，单位 $e/2\pi$）。因此 **OHC 的符号是检测 $L$ 符号错误的最灵敏测试**，应作为 CI 的必过项。

---

## 6. 论文给出的验证基准

论文用 Haldane 模型（Eq. 19）作检验场，参数 $t_1=1$、$t_2=0.15$、$M=0.4$，$\phi$ 扫描 $0\to0.5\pi$。可复现的定量目标：

| 量 | 论文位置 | 本库应达到 |
|---|---|---|
| Chern 数拓扑相变 | Fig. 1 inset | $\phi_c=\arcsin[M/(3\sqrt3 t_2)]/\pi\approx0.1715\pi$ 处 $0\to-1$ |
| $L_{12}(\mathbf k)$ 沿 $\Gamma$–K | Fig. 1 上 | 峰值约 $0.04\hbar$，$\phi=0$ 时为零 |
| $L_{11}(\mathbf k)$ | Fig. 1 下 | $\phi=0$ 时 $L_{11}(k)=-L_{11}(-k)$；量级 $\le1\hbar$ |
| $\sigma^{\rm OH}_{xy}(\phi)$ | Fig. 2 | 在相变点发散，非拓扑区随 $\phi$ 变化 |
| $\sigma^{\rm OH}_{xy}(E_F)$ | Fig. 3 | $\phi=\pi/2$ 时在 $E_F=0$ 附近形成平台，幅值约 $0.28\,e/2\pi$ |

> **说明**：Fig. 1 的 $|L_{12}|$ 峰值与 Fig. 2/3 的绝对幅值与 $k$ 网格密度、$E_F$ 取值、以及 $\phi$ 的具体点有关。本库的验收标准不是"像素级复现论文图"，而是"三条独立计算路径互相一致到 $10^{-12}$，且对称性与拓扑约束成立"。论文数值只作为量级与形状的 sanity check。

---

## 7. 适用边界

论文强调其形式对**金属与绝缘体都成立**。实现时需要区分：绝缘体时 $\mu$ 在带隙内，$\sigma^{\rm OH}$ 是费米海积分，无发散；金属时存在部分占据带，$(\varepsilon_n-\varepsilon_m)^{-2}$ 在带交叉处奇异（这是 Fig. 2 在相变点发散的原因）。库应显式检测近简并并给出告警阈值，而不是静默返回大数。

此外 Eq. (16b) 的求和形式在 $N_b$ 大时是 $O(N_b^3)$ 每 k 点；Eq. (16a) 的投影形式同样是 $O(N_b^3)$。对 $N_b\lesssim 100$ 两者都可接受；更大体系需要考虑 $k$ 点并行与只算占据带附近的窗口。
