# Orbital Hall kit (ohk)

按照原设计的 L0–L5 分层实现的 Python 科学计算库：从紧束缚或 Wannier90 Hamiltonian 计算轨道角动量矩阵与轨道霍尔电导，并提供独立的有限周期格点参考计算。

基于 Gautam & Satpathy, [Phys. Rev. B **114**, 214403 (2026)](https://doi.org/10.1103/pd3s-q2z7)。原设计中的分母转写、相位、周期边界与单位问题已校订，详见 [实现与验收记录](docs/design/05-implementation-notes.md)。

## 实现状态

| 阶段 | 已实现 |
|---|---|
| P0 | 通用跃迁模型、完整 Haldane 模型、解析导数、本征基、FHS Chern 数 |
| P1 | 投影形式 OAM、独立有限周期格点参考计算 |
| P2 | 谱求和、Pezo 截断、ACA，以及三／四带交叉验证 |
| P3 | Kubo/BZ 积分、费米能扫描、网格收敛序列、显式展宽 |
| P4 | JSON 报告、CSV 路径扫描、三个命令行入口 |
| P5 | Wannier90 原生读写、简并权重、三维模型、真实 Fe 数据端到端冒烟验证 |

真实材料测试验证文件与计算链路。材料预测仍需检查网格、展宽收敛及位置算符近似。

## 安装与验证

Python 3.10 或更新版本：

~~~bash
python -m pip install -e ".[dev]"
python -m pytest
ohk-validate --n-k 30 --output validation.json
~~~

当前 Windows 工作区已准备好项目环境，也可以直接运行：

~~~powershell
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\ohk-validate.exe --n-k 30
~~~

## Python API

~~~python
import numpy as np
from ohk import (
    HaldaneModel, EigenBasis, CovariantProjector,
    BZGrid, full_matrix, orbital_hall_conductivity, fermi_scan,
)

model = HaldaneModel(t1=1.0, t2=0.15, M=0.4, phi=np.pi / 4)
variant = CovariantProjector()
basis = EigenBasis.at(model, [0.7, 1.1])
L, diagnostics = full_matrix(variant, basis)

grid = BZGrid.uniform(model.lattice(), n=60)
result = orbital_hall_conductivity(model, variant, grid, mu=0.0)
print(result.sigma, result.units)
# -1.0616491612557013, e/(2*pi)

scan = fermi_scan(model, variant, grid, [-0.5, 0.0, 0.5])
~~~

默认 `variant.matrix(basis)` 的投影／截断变体只提供非对角部分；用 `full_matrix` 填充对角元并检查厄米性。也可使用 `CovariantProjector(include_diagonal=True)`。

`k` 使用笛卡尔坐标，晶格基矢按行排列，`v[j,n,m]` 的首个下标是方向。Haldane 模型包含胞内轨道位置；`AnalyticModel` 的位置默认为零，可用 `orbital_positions` 指定。

## 命令行

~~~bash
ohc-run --phi 0.7853981633974483 --n-k 60 --mu 0 --output ohc.json
ohc-run --n-k 30 --fermi -0.5 0 0.5 --convergence 20 40 60
oam-scan --points-per-segment 80 --output oam.csv
ohk-validate --n-k 30 --output validation.json
python examples/haldane.py --n-k 30 --output output/haldane.json
~~~

`--phi` 使用弧度。`--variant` 可选 `covariant`、`spectral`、`pezo`、`aca`。路径扫描默认 Γ–K–K′，也可用 `--vertices "0,0;0.333333,0.666667"` 指定倒空间分数坐标。CSV 包含能带与全部 OAM 矩阵元的实、虚部。

## Wannier90 与材料示例

`_hr.dat` 不包含晶格，必须另提供按行排列基矢的文本矩阵：

~~~bash
ohc-run --hr material_hr.dat --lattice lattice.txt --n-k 12 --mu 9.1839 --output material.json
~~~

可下载固定版本的公开 WannierTools Fe 示例。脚本验证压缩包和数据文件的 SHA-256，仅保存在不纳入版本控制的 `.local/fe/`：

~~~bash
python scripts/fetch_fe_example.py
ohc-run --hr .local/fe/wannier90_hr.dat --lattice .local/fe/lattice.txt --n-k 3 --mu 9.1839 --eta 0.01 --output .local/fe/smoke.json
~~~

此命令是粗网格冒烟测试。`--eta` 显式启用正则化。Hamiltonian-only 入口采用胞周期规范与对角位置近似，未导入 `_r.dat`、完整异常位置矩阵或 WS-distance 修正，不能直接等同于完整第一性原理 OHC。Python API 另提供受限的 legacy PythTB 自旋无关、全周期模型适配器。

## 数值与单位

- 内部设 `hbar=m_e=e=1`。非对角 OAM 随本征态相位**协变**；对角期望及 Kubo 收缩不变。
- 电导数值按二维 `e/(2*pi)`、三维 `e/(2*pi*length)` 表示。SI 单位分别为 C、C/m，不是电荷电导的 S。
- eV/Å 模型恢复 `L/hbar` 还需因子 `m_e E0 a0²/hbar²`；用 `ohk.units.UnitSystem` 指定尺度。
- 精确简并默认抛出 `DegeneracyError`；近简并保留原值并标记。展宽是显式近似，不默认截断。
- `convergence_study` 返回实测序列及相邻差，不从少量网格自动声称收敛或给出外推。
- 完整 Haldane 模型的 `d0` 移动带隙中心，`mu=0` 不保证半填充绝缘态。相变发散研究须选择位于带隙内的化学势。

## 论文关键图复现

```powershell
.\.venv\Scripts\python.exe scripts/reproduce_paper.py --output output/paper-reproduction
```

生成图 1 的 OAM 路径、图 2 的通量扫描和图 3 的费米能扫描，同时输出 PNG/SVG、原始 CSV 和收敛表。
默认对照 240×240 Bloch 网格与独立 60×60 有限周期格点。通量方向、eV/Å 尺度假设以及相变附近的
收敛局限，见 [复现说明](docs/reproduction.md)。新环境需先安装 `python -m pip install -e ".[plot]"`。

## 结构与测试

[03 架构设计](docs/design/03-architecture.md) 保留原始分层；历史记录中与实现冲突的内容以 [05 校订记录](docs/design/05-implementation-notes.md) 为准。

测试覆盖三路径公式对照、三／四带负对照、相位与对称性、有限格点、Chern 相变、独立 Kubo 双重求和、扫描复用对角化、Wannier90 往返与错误输入、单位、报告和 CLI。CI 配置覆盖 Linux/Windows 与 Python 3.10/3.13。外部材料下载不作为离线测试依赖。

MIT，见 [LICENSE](LICENSE)。
