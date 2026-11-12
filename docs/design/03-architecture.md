# 03 · 算法库架构设计

> **历史设计记录：** 本文保留原始设计。关于 Eq. (16b) 符号、相位、Pezo 截断、周期位置算符、单位及部分基准值的结论已在实现阶段校订，以 [05-implementation-notes.md](05-implementation-notes.md) 和可执行测试为准。

目标：一个可复用的 Python 库，输入任意紧束缚 / Wannier90 模型，输出规范不变的 OAM 矩阵元与轨道霍尔电导，并以有限格点实空间计算作为自验证基准。

设计的第一原则来自理论分析的直接结论：**$L_{nn'}$ 不是一个唯一确定的公式，而是一族定义。** 论文 Eq. (16) 只是其中之一（且印刷形式有误）。因此架构必须把"$L$ 的定义"当作可插拔对象，而不是把公式写死在某个函数里。

---

## 1. 分层

```
┌──────────────────────────────────────────────────────────┐
│  L5  CLI / 脚本层        ohc-run, oam-scan, validate      │  ← 用户入口
├──────────────────────────────────────────────────────────┤
│  L4  报告层              Report: tensors, figures, checks │
├──────────────────────────────────────────────────────────┤
│  L3  输运层              Kuboχ: conductivity, Fermi scan  │
├──────────────────────────────────────────────────────────┤
│  L2  OAM 内核层          OAMVariant 族 + 组装器            │  ← 本库核心
├──────────────────────────────────────────────────────────┤
│  L1  能带量层            EigenBasis: ε, v, ∂H, C_j         │
├──────────────────────────────────────────────────────────┤
│  L0  模型层              Model → H(k), dH(k)              │
└──────────────────────────────────────────────────────────┘
```

依赖严格向下：L2 只依赖 L0/L1 的公共接口，L3 只依赖 L2 产出的 $L$ 与 L1 的 $v$。这样任何一层都可以单独替换（例如把 L0 换成 WannierBerri 后端）。

**关键的数据契约是 L1**：几乎所有计算都需要"给定 $\mathbf k$，给出 $\varepsilon$、$v^j_{nm}$、$\partial_{k_j}\mathcal H$ 或等价的 $C_j$"。把这个对象定义清楚，上面所有层都只跟它打交道。

---

## 2. 模块与职责

### L0 `ohk.model`

| 类/函数 | 职责 |
|---|---|
| `TightBindingModel`（抽象基类） | 唯一必需方法：`H(k) -> (N_b,N_b) ndarray`；`dH(k) -> (2,N_b,N_b)`；`lattice()` → 实空间基矢；`num_bands` |
| `AnalyticModel` | 从跃迁列表 `(i, j, R, amplitude)` 构造，解析给出 `H` 与 `dH`。Haldane 等手搭模型用它 |
| `Wannier90Model` | 读 `_hr.dat`（或委托 PythTB `W90`），同样给出解析 `dH` |
| `FiniteLattice` | 构造 $N\times N$ PBC 实空间 $\mathcal H$、位置矩阵、速度矩阵。**只供 L5 验证使用** |
| `HaldaneModel` | 论文基准模型的内置实现，含 $d_0(\mathbf k)$ 全项（见理论文档陷阱 3） |

设计要求：`dH` 必须是解析的。基类提供 `dH_finite_difference` 作为兜底并在文档中标注"仅当模型无法解析求导时使用"。

### L1 `ohk.eigen`

| 类/函数 | 职责 |
|---|---|
| `EigenBasis`（数据类） | 一次 `at(k)` 调用的产物：`eps`、`U`、`v`（$v^j_{nm}$）、`dH`；惰性缓存 `C_j` |
| `derivative_coefficients(eps, v)` | 构造 $C_j$（见理论 §4.3） |
| `velocity_matrix(H, dH, U)` | $v^j = U^\dagger \partial_{k_j}H\,U$ |
| `berry_connection(eps, v)` | 可选，仅当某变体需要 $\mathcal A_n$ 时计算；主路径不需要 |

**规范处理的硬约定**：`EigenBasis` 保证 `U` 只用于构造规范不变量（$v$、$L$、$\sigma$），任何下游都不得直接使用 `U` 的相位。`berry_connection` 显式标注"相位敏感"，并在测试中检查它随规范跳变而 $L$ 不变。

### L2 `ohk.oam` — 核心

这是本库区别于一切现有工具的一层。

```
oam/
├── base.py           OAMVariant 抽象类
├── covariant.py      Eq.(16a) 投影形式（主路径）
├── spectral.py       Eq.(16b) 求和形式（修正符号）
├── atomic.py         ACA：原子中心近似（作为反例/对照）
├── pezo.py           Pezo-Manchon 截断式（论文 ref 26）
└── assemble.py       组装完整的 L 矩阵（对角 + 非对角）
```

`OAMVariant` 抽象接口：

```python
class OAMVariant(ABC):
    name: str
    @abstractmethod
    def matrix(self, basis: EigenBasis) -> np.ndarray:
        """返回 (N_b, N_b) 复矩阵 L_{nn'}，单位 hbar=1, m_e=1。"""
```

四个具体变体的差异仅在于 `matrix` 的实现：

| 变体 | 公式 | 用途 |
|---|---|---|
| `CovariantProjector` | Eq. (16a) 矩阵形式（§4.3） | **主路径**，规范不变，$O(N_b^3)$ |
| `SpectralSum` | Eq. (16b) 修正符号 + 对角闭式 | 交叉验证；对小 $N_b$ 更快 |
| `AtomCentered` | $\sum_i \delta_{\text{cell}}\,\mathbf L_i$ | 作为物理反例：对 s 轨道系统恒为零（论文 §IV 的关键对照） |
| `PezoTruncated` | 只保留 Eq. (16b) 第一项 | 复现论文指出的"二带 $L_{12}\equiv0" 缺陷 |

**为什么要把错误的变体也实现**：论文的核心论点是"不同文献给出不同结果"。把这些变体都实现，才能在同一模型上定量展示差异（例如 PezoTruncated 对 Haldane 给出 $L_{12}=0$ 而主路径给 0.04）。这既是科研价值，也是最强的回归测试。

**组装规则**：`assemble.full_matrix()` 负责"对角用闭式、非对角用变体"，并强制厄米化检查（`L - L†` 的范数超过阈值即报错）。

### L3 `ohk.transport`

| 类/函数 | 职责 |
|---|---|
| `orbital_current(L, v, axis)` | $J_x = \frac12(v_xL_z+L_zv_x)$，对称化 |
| `KuboConductivity` | 输入 `L_variant`、`velocity`、占据判据，输出 $\sigma^{\rm OH}_{xy}$ |
| `fermi_scan(L_variant, energies)` | 对 $E_F$ 扫描，返回数组 |
| `sum_rule_check` | 检查 $\mathrm{tr}\,L=0$、厄米性、对称性 |

Kubo 积分器的网格策略与"网格生成"分离：`BZGrid` 对象负责产生 $\mathbf k$ 点与权重（支持均匀 Monkhorst-Pack、Haldane 的非标准格子、以及任意自定义网格）。这是论文 Fig. 2/3 能复现的前提。

### L4 `ohk.report`

把张量、检查结果、可选的图打包成 `Report`。设计上让"计算"与"呈现"解耦：`Report` 只吃纯数据，不 import matplotlib（matplotlib 作为可选依赖）。这样 CI 可以只跑计算不画图。

### L5 `ohk.cli`

三个命令：`ohc-run`（算单个模型的 $\sigma$）、`oam-scan`（沿 k 路径输出 $L_{nn'}$）、`validate`（跑全套一致性检查）。CLI 是薄封装，所有逻辑在 L0–L4。

---

## 3. 数据流（以一次 OHC 计算为例）

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

一次 `EigenBasis` 构造被同一 k 点上的所有变体复用——这是"多变体并列对比"能便宜实现的原因：对角化只做一次。

---

## 4. 张量与单位约定

全库内部**不做单位换算**，统一 $\hbar=m_e=e=1$，长度以晶格常数为单位。理由：论文公式里 $m_e/\hbar$ 因子反复出现，混用极易出错；而在 $e=\hbar=1$ 下库仑单位制与 SI 的换算只在输出层做一次。

| 量 | 内部约定 | 输出换算 |
|---|---|---|
| $L$ | $\hbar$ 为单位，`hbar=1` | 乘 $\hbar$ |
| $v$ | $[\text{能量}]\cdot[\text{长度}]$，`hbar=1` | 乘 $1/\hbar$ |
| $\varepsilon$ | 模型自然单位 | 乘能量单位 |
| $\sigma^{\rm OH}$ | $e/2\pi$ 为单位（二维），$A$ 为胞面积 | 乘 $e/2\pi$ 或 $e/(2\pi a)$ |

张量下标约定：`L[n, n']`、`v[j, n, m]`（$j\in\{x,y\}$）、`dH[j, n, m]`。**统一 $v$ 的下标顺序为 `[方向, 带, 带]`**——这与 `EigenBasis` 一致，避免"有时 $v_{nm}$ 有时 $v_{mn}$"的混乱。

---

## 5. 数值策略

**对角化**：`numpy.linalg.eigh`，每 k 点一次。对 $N_b\le 200$ 无压力。若未来需要更大体系，`scipy.linalg.eigh` 的 `subset_by_index` 可只算占据带附近窗口。

**简并处理**：$L_{nn'}$ 分母 $\varepsilon_n-\varepsilon_m$ 在带简并处发散。策略：

1. 每个 k 点计算最小能隙 $\delta_{\min}$；
2. 若 $\delta_{\min}<\tau_{\text{deg}}$（默认 $10^{-8}$），将该点标记为 `degenerate`；
3. `Report` 显式列出简并点占比；超过阈值（默认 0.1%）时给出警告；
4. 可选提供"简并子空间内正则化"——在简并子空间内对角化 $\mathbf v$ 的某个组合，取本征值差的倒数截断。**这一步不设为默认**，因为它改变了物理量的定义；默认行为是"如实报告"。

**近相变点**：论文 Fig. 2 在 Chern 数翻转处 $\sigma$ 发散，这不是数值问题而是物理。库应在 \(\phi\) 接近 $\phi_c$ 时给出"near topological transition"标记，而不是试图压平。

**k 网格收敛**：提供 `convergence_study(grid_sizes)`，返回 $\sigma$ 随网格的变化与外推值。论文 Fig. 2 的实线是 $k$ 连续极限，圆圈是 $60\times60$ 格点；我们的测试应记录两者的差，作为"有限格点 vs 连续"的定量对照。

**BZ 积分**：默认等权 Monkhorst-Pack。对 Haldane 的蜂窝格子，`BZGrid` 必须用正确的倒格矢（由实空间基矢 $2\pi A^{-T}$ 得到）——这里极易写错，测试需包含"$a_i\cdot b_j=2\pi\delta_{ij}$"的自检。

---

## 6. 可插拔点（扩展位）

| 扩展位 | 接口 | 未来用途 |
|---|---|---|
| 新模型格式 | 子类化 `TightBindingModel` | TBmodels HDF5、pyqula 模型、非共线自旋 |
| 新 $L$ 定义 | 子类化 `OAMVariant` | Sastges-Baek 式、含反常位置的 Go-Lee 式 |
| 新响应函数 | L3 增加 `KuboX` 子类 | 轨道 Edelstein、轨道磁化、自旋 Hall 对照 |
| 新积分器 | 替换 `BZGrid` | 自适应网格、四面体法、WannierBerri 后端 |
| 新后端 | 实现 `TightBindingModel` 适配 WannierBerri `System` | 大规模第一性原理体系 |

这些接口的存在是"通用引擎"与"一次性脚本"的分界线。

---

## 7. 包结构

```
orbital-hall/
├── pyproject.toml
├── README.md
├── docs/design/                    01-theory · 02-tooling-survey · 03-architecture · 04-validation
├── src/ohk/
│   ├── __init__.py                 __version__，顶层导出
│   ├── units.py                    hbar, m_e, 换算因子
│   ├── model/
│   │   ├── base.py                 TightBindingModel, ModelSpec
│   │   ├── analytic.py             AnalyticModel（跃迁列表 → H, dH）
│   │   ├── haldane.py              HaldaneModel（基准）
│   │   ├── wannier90.py            Wannier90Model（_hr.dat）
│   │   └── finite.py               FiniteLattice（oracle 用）
│   ├── eigen.py                    EigenBasis, derivative_coefficients
│   ├── oam/
│   │   ├── base.py                 OAMVariant
│   │   ├── covariant.py            Eq.(16a)
│   │   ├── spectral.py             Eq.(16b) 修正
│   │   ├── atomic.py               ACA 反例
│   │   ├── pezo.py                 Pezo 截断式
│   │   └── assemble.py             对角+非对角组装
│   ├── transport/
│   │   ├── current.py              轨道流算符
│   │   ├── kubo.py                 KuboConductivity, fermi_scan
│   │   └── grids.py                BZGrid
│   ├── validation/
│   │   ├── checks.py               厄米、trace、对称性
│   │   ├── oracle.py               有限格点对照
│   │   └── benchmarks.py           Haldane 基准数值
│   ├── report.py                   Report
│   └── cli.py
└── tests/
    ├── conftest.py                 Haldane / 随机模型 fixtures
    ├── test_model.py
    ├── test_eigen.py
    ├── test_oam_consistency.py     三路径互相对照（核心）
    ├── test_symmetry.py            T / I 约束
    ├── test_oracle.py              有限格点 vs Bloch
    ├── test_transport.py           OHC 符号与量级
    └── test_units.py
```

命名理由：包名 `ohk`（orbital Hall kit）短且不冲突；模块按"物理对象"而非"数学操作"划分，使新用户能凭物理直觉找到代码。

---

## 8. 测试策略（架构的一部分）

本库的正确性不能靠"论文图看起来对"来保证。测试必须建立在**多条独立计算路径互相印证**上：

| 测试 | 内容 | 预期 |
|---|---|---|
| `test_oam_consistency` | 同一 k 点：`CovariantProjector` vs `SpectralSum` vs 逐点定义 | 差 $<10^{-12}$ |
| `test_oracle` | 有限格点实空间 $L$ vs Bloch 公式 $L$ | 在可用 k 点差 $<10^{-10}$ |
| `test_symmetry` | $\phi=0$ 时 $L_{11}(k)=-L_{11}(-k)$；$|L_{12}|$ 在 $\phi=\pi/2$ 时为零 | $<10^{-12}$ |
| `test_degenerate_bands` | 二带模型（$L_{12}$ 无法检验求和项）→ 必须用 3 带随机模型 | 3 带残差 $<10^{-12}$ |
| `test_transport` | $\sigma^{\rm OH}$ 符号为正/负的定标 | 修正版为负，印刷版为正 |
| `test_chern` | Haldane Chern 数在 $\phi_c$ 翻转 | $0\to-1$ |
| `test_units` | $\hbar$ 恢复后 $L$ 的量级 | 与 $\hbar$ 同量级 |

**最重要的一条**：`test_degenerate_bands` 中那个"3 带随机模型"测试，正是因为二带模型对 Eq. (16b) 的求和项是盲的。任何只用 Haldane（二带）做验证的实现，都会漏掉论文印刷版的符号错误。这条测试必须写进 CI 且不可删除。

---

## 9. 分阶段实施（本次不实现）

| 阶段 | 内容 | 验收 |
|---|---|---|
| P0 | L0 模型层 + L1 `EigenBasis` + 测试 | Haldane 能带、Chern 数与论文一致 |
| P1 | L2 `CovariantProjector` + 有限格点 oracle | 三路径一致性 $<10^{-12}$ |
| P2 | L2 其余变体 + 差异展示 | PezoTruncated 在 Haldane 上 $L_{12}=0$ |
| P3 | L3 Kubo + BZGrid | $\sigma^{\rm OH}(\phi)$ 复现论文 Fig. 2 形状与符号 |
| P4 | L4 Report + L5 CLI | `validate` 一键跑全套检查 |
| P5 | Wannier90 接入 | 用真实材料 `_hr.dat` 跑通 |

本次交付只到"骨架"：接口、数据类、抽象基类、模块 docstring 就位，函数体为 `NotImplementedError`。P0 的实现可以在下一轮直接填充，不需要改架构。
