# 02 · 现有计算工具调研

结论先行：**没有一个现成工具可以直接算出论文 Eq. (16) 定义的 $L_{nn'}(\mathbf k)$ 全套矩阵元并配套做验证。** 但在"造轮子"和"从零写"之间存在一个明确的最优路径：用成熟的紧束缚/Wannier 生态提供 $\mathcal H(\mathbf k)$ 与 $\partial_{k_j}\mathcal H(\mathbf k)$，自己实现 $L$ 与 OHC 这一层。下面按这个判断组织。

## 1. 生态全景

Wannier 社区维护了一份自描述的代码注册表（`wannier-developers/wannier-ecosystem-registry`，配综述 arXiv:2312.10769），按功能分类列出所有与 Wannier 相关的软件。其中"线性响应"一栏的代码可以提供 Kubo 框架下的 AHE/SHE 等，**但没有一项列出轨道霍尔效应或轨道角动量矩阵元**。这是本调研最重要的一条事实。

## 2. 逐个评估

### 2.1 Wannier90

用途：从第一性原理构造最大局域化 Wannier 函数，输出 `_hr.dat`（实空间 $H(\mathbf R)$）、`_mmn`、`_eig`、`_chk`、`_uHu`、`_uIu` 等。

对本项目的价值：**很大，但不是计算引擎**。`_hr.dat` 里的 $H(\mathbf R)$ 就是我们需要的全部输入——$\mathcal H(\mathbf k)$ 与其解析导数都能从它得到。后处理模块 `wannier90.x -berry`（或独立 `postw90`）能做 AHE、SHE、Berry curvature。

关键局限：官方 Berry 模块**没有 OHE 的 current vertex**。它内置的流算符是自旋流或电流，轨道流 $J_x=\frac12(v_xL_z+L_zv_x)$ 需要自己补。另外 `_uHu` 文件（$\langle u_{mk}|\hat H|u_{n,k+b}\rangle$）是给轨道磁矩用的，与论文的 $L$ 不是同一个量。

实用提示：`uHu_formatted = .true.` 可让 Wannier90 输出文本版 `_uHu`，便于核对文件布局；Wannier90 也不存在完整的逐列格式规范页，实际布局需从生成器源码确认。

### 2.2 WannierBerri

用途：高性能 Wannier 插值与 BZ 积分，内置 Berry curvature、AHC、SHC、轨道磁化等计算器。

对本项目的价值：**A 级参考实现与可借用的后端**。

- `wannierberri.utils.mmn2uHu` 从 `_mmn`、`_eig` 生成 `_uHu`，正文明确说是"计算轨道磁矩所需的矩阵元"。这说明它的轨道磁矩计算走的是 $\langle u|\mathbf r\times\mathbf v|u\rangle$ 的现代理论路线，与论文 §2.2 的规范协变观点同源。
- `spin_Hall` 计算器接受 `spin_current_type` 参数（`"ryoo"` 或 `"qiao"`），并支持对称化——这正是我们要在 OHE 上复刻的架构模式：**一个通用的"流-顶点 + 对称化"计算器，把流算符当参数注入**。
- 它把"公式层"（`formula/covariant.py`）与"计算器层"（`calculators/`）分离，还把 `spin_Hall` 的流算符 $J$ 单独抽象。这个分层是我们架构设计的直接蓝本。

**是否已有 OHC？** 检索未发现 WannierBerri 文档列出轨道霍尔计算器。它有能力算 $L$（轨道磁矩路径）也有能力算 Kubo 响应（SHC 路径），但把两者拼成 OHC 的 current vertex 那一层，公开文档里没有。因此本项目要做的是**补上这个连接层**，而不是重写插值器。

### 2.3 WannierTools

用途：基于 Wannier90 紧束缚模型做拓扑性质计算（表面态、Wilson loop、AHC/SHC 等）。

关键事实：**GitHub releases（约 2025-06 的 v2.7.2）明确写了 "Added orbital Hall effects calculations"**。这是目前唯一有公开 OHE 计算能力的成品软件。它读 `wannier90_hr.dat` 加 `wt.in`（NAMELIST + INPUT_CARDS 格式）。

对本项目的价值与定位：它是**竞品也是参照**。它的存在说明 OHE 在 Wannier 体系里是成熟需求；但它基于 Fortran + 固定输入文件，不适合做"论文公式的独立验证"——用它验证等于用一个未开源的实现去检验另一个。我们的定位应是：**以论文 Eq. (16) 为准的 Python 参考实现 + 有限格点 oracle**，WannierTools 的 OHE 结果可作为第三方交叉校验（若其 vertex 定义与论文一致）。

### 2.4 PythTB

用途：Python 紧束缚模型构造与拓扑/量子几何分析。

- `tb_model`/`TBModel`/`W90` 三层：`W90("path").model` 直接从 Wannier90 输出读模型，`W90._read_hr_dat` 内部处理 `_hr.dat`；`min_hopping_norm` 可裁剪小跃迁。
- 有 `_model.get_hoppings()` 之类的接口可拿到跃迁列表，位置由 `lattice()`/`sites()` 提供。

局限：**没有任何 OHE/OAM 功能**，连 Kubo 输运都没内置。它的价值是**模型定义的便利性**——尤其适合写 Haldane 这类手搭模型，以及作为"用户自定义模型的输入格式"。它的 `W90` 读取器可作为我们读 `_hr.dat` 的参考（或直接用）。

### 2.5 TBmodels

用途：从 Wannier90 输出构造紧束缚模型，存 HDF5。

价值：`Model` 类直接读 `*_hr.dat`（需 Wannier90 侧 `write_hr = .true.`），提供 `eval`/`bands`/`wannier_centers`。适合作为"从文件到模型对象"的轻量适配层。与 PythTB 二选一即可——我们倾向于 **PythTB 的 `W90` 路径**，因为它同时提供了手搭模型的 API（Haldane 基准需要）。

### 2.6 其他相关代码

| 代码 | 与 OHE 的关系 | 结论 |
|---|---|---|
| **Rustb** (docs.rs/Rustb) | 文档明确写"calculates the orbital Hall conductivity and the orbital angular moment" | 有 OHC，但 Rust + 文档极简，输入参数体系不明；作为第三方交叉验证的候选 |
| **pyqula** (joselado) | Python TB 库，有 Kubo/Hall 类线性响应 | 可借鉴其 Kubo 实现；无 OAM 顶点 |
| **TBPLaS 2.0** | Kubo-Bastin 大体系 Hall 输运 | 方法（实空间 Kubo-Bastin）与我们的 k 空间路线互补，可用于强无序体系 |
| **TightBinding++** | Kubo-Greenwood 电导 | 同上，参考 |
| **Go & Lee, PRB 109, 174435 (2024)** | 论文核心文献：Wannier 插值算 OHC，强调"反常位置的轨道依赖" | **未见公开代码仓库**。但其物理结论（位置算符的轨道依赖修正不可忽略）直接支持论文的规范协变处理 |

### 2.7 文献中的平行形式（重要）

调研发现两篇与论文几乎同期的工作，它们给出的 $L$ 公式**互不相同**——这正是论文引言强调的"conflicting results"：

- **Pezo, Ovalle, Manchon, PRB 106, 104414 (2022)** [论文 ref 26]：非对角项只保留第一项（Eq. 16b 的第一项）。论文明确指出：**该式对任意二带系统给出恒为零的 $L_{12}$**，而 Haldane 直接计算显示 $L_{12}\ne0$。
- **Cysne, Souza, Rappoport, PRR 8, 033086 (2026)** [论文 ref 27]：给出的式子在最终形式上与论文 Eq. (16) 一致（论文自述"agrees with [27]"）。
- **Sastges & Baek, arXiv:2604.08280**：从轨道磁化现代理论出发的另一种 Wannier 基 OAM 构造。

这对架构的含义：**"用哪条公式"是必须显式建模的一等公民，不能硬编码单一公式。** 库必须支持多个 $L$ 定义（至少 paper-16a、paper-16b-corrected、Pezo-style-truncated、Rappoport-style），并能在同一组 $\mathbf v$、$\varepsilon$ 上并列计算、对比。这是本库区别于"一次性脚本"的核心设计决策。详见架构文档的 `oam_variants` 模块。

## 3. 能力缺口汇总

| 能力 | 现成方案 | 缺口 |
|---|---|---|
| $\mathcal H(\mathbf k)$ 与解析 $\partial_k\mathcal H$ | Wannier90 `_hr.dat`、PythTB、TBmodels | 无 |
| Wannier 插值、BZ 积分、k 并行 | WannierBerri | 无（可直接借用） |
| 轨道磁矩 $\langle u|\mathbf r\times\mathbf v|u\rangle$ | WannierBerri `mmn2uHu`/轨道磁化路径 | 只有对角（带内） |
| **OAM 全套矩阵元（含非对角）** | **无通用实现** | **本库核心** |
| **轨道流顶点 $J_x=\frac12(v_xL_z+L_zv_x)$** | **无** | **本库核心** |
| **OHC Kubo 积分** | WannierTools（闭源实现）、Rustb | 需自建参考实现以做独立验证 |
| **有限格点实空间 oracle** | 无 | 本库验证层 |
| 多套 $L$ 公式并列对比 | 无 | 本库设计特色 |

## 4. 技术选型建议

`Python 3.10+`，`numpy` + `scipy` 起家；`PythTB` 作为可选依赖（读 `_hr.dat`、手搭模型），不设为硬依赖以免绑定其数据结构。BZ 积分先自写（Haldane 需要非标准格子与不规则网格），大规模体系再考虑接 WannierBerri 的积分器。不做 C++/Fortran 内核——$N_b\lesssim100$ 时 NumPy 的 $O(N_b^3)$ 矩阵乘足够，先把正确性做对。

`pytest` 作为测试框架；`h5py` 可选（缓存大模型的 $\mathcal H(\mathbf k)$ 采样）。
