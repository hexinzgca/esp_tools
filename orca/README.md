# ORCA：wB97M-V + 弥散 def2 基组 + 10×原子数 ESP 格点

## 1. 基组

默认使用 `def2-TZVPD`（三-zeta + 弥散，OMol25 同级别）；需要双-zeta + 弥散
时使用 `def2-SVPD`。

| 关键字 | 级别 | 用途 |
|---|---|---|
| `def2-SVPD` | 双-zeta + 弥散 | 备选 |
| `def2-TZVPD` | 三-zeta + 弥散 | **默认** |
| `def2-TZVPPD` | 三-zeta 双极化 + 弥散 | 更高精度 |
| `def2-QZVPD` | 四-zeta + 弥散 | 参考值 |

本目录提供两个可直接运行的输入：

- `wb97m-v_def2-tzvpd.inp`（默认）
- `wb97m-v_def2-svpd.inp`（双-zeta 弥散）

抽取脚本与基组无关。

## 2. SCF 输入

`wb97m-v_def2-tzvpd.inp`：

```text
! wB97M-V def2-TZVPD AutoAux RIJCOSX VeryTightSCF

%pal nprocs 4 end
%maxcore 4000

%method
  Grid 5
end

* xyzfile 0 1 water.xyz
```

- `AutoAux` 自动生成辅助基；`RIJCOSX` 加速交换。
- `VeryTightSCF`：ESP 对密度收敛敏感，不建议用默认收敛限。
- `Grid 5`：与现有项目一致；追求更高精度可改 `Grid 6`（更慢）。
- `xyzfile` 的电荷/多重度按体系修改；阴离子/阳离子同时改 `--gbw` 对应的 SCF 输入。
- ORCA 6.1 在 SCF 后自动写 `water.gbw` 和 `water.densities`，`orca_vpot` 需要两者。
- 默认 `%pal nprocs 4`；按机器核数调整。`run_orca_esp.sh` 设置了
  `OMPI_MCA_rmaps_base_oversubscribe=1`，避免请求核数超过物理核时直接失败。

运行（容器）：

```bash
./run_orca_esp.sh wb97m-v_def2-tzvpd.inp water.xyz
```

该脚本会：

1. 在 `/mnt/e/apptainer_sifs/build/orca-6.1.sif` 中运行 ORCA；
2. 生成 `wb97m-v_def2-tzvpd.out`、`.gbw`、`.densities`；
3. 调用 `extract_esp.py` 采样 10×3=30 个点并写出 `*_esp.json`。

> ORCA 并行时必须用 `orca` 的**完整路径**启动；脚本内部用
> `ORCA=$(command -v orca); exec "$ORCA"` 满足该要求。

## 3. 采样

`extract_esp.py` 使用共享的 `common/sampling.py`：

- 总点数 = `10 × natoms`，严格保证（水：30 点）；
- 四箱按最近原子距离 `[0.25,0.5) [0.5,1.5) [1.5,5) [5,∞)` Å，比例 `1:4:4:1`；
- 第一箱只占 10%，且最小距离 0.25 Å，**不会把大量点堆在核附近**；
- 远场点位于质心球壳 `r_max+5 … r_max+10` Å；
- 种子 = 固定基值 + 几何哈希：同一结构可复现，不同结构不共享采样模式。

这与旧 COWORK 版本按 DFT 积分壳层取点的做法不同；后者会把点推进核区。本实现不生成
DFT 积分格点，只生成用于 `orca_vpot` 评估的探针集。

## 4. 抽取

```bash
# 采样 + 用已有 GBW 提取（容器内 ORCA）
python orca/extract_esp.py \
  --xyz orca/water.xyz \
  --gbw /path/to/water.gbw \
  --orca-container /mnt/e/apptainer_sifs/build/orca-6.1.sif \
  --output water_esp.json

# 原生 ORCA 安装（目录内需有 orca_vpot）
python orca/extract_esp.py --xyz orca/water.xyz --gbw water.gbw \
  --orca-home /opt/orca/orca-6.1.0-f.0_linux_x86-64/bin --output water_esp.json

# 只采样，不做 SCF
python orca/extract_esp.py --xyz orca/water.xyz --output water_grid.json

# 同一组探针复用于不同电荷态
python orca/extract_esp.py --points water_grid.json --gbw anion.gbw \
  --orca-container /path/orca-6.1.sif --output anion_esp.json
```

`orca_vpot` 的密度名不靠猜：先运行 `orca_vpot water.densities` 列出容器内容，
再取唯一的 `.scfp`（本机 ORCA 6.1.0 实测可用）。若目录中是独立 `.scfp`（旧版本），自动回退。
网格和原始结果写在 GBW 同目录，需有写权限。

`--xyz` 必须与 GBW 的核坐标一致；`--points` 复用网格时脚本不做几何校验，必须保证该网格
来自同一结构（例如同一几何的 neutral/anion/cation 三次 SCF）。

## 5. 输出与单位

```text
esp.positions_Angstrom
esp.phi_molecular_eV_per_e        核 + SCF 电子
esp.phi_background_eV_per_e       背景点电荷（若 --background）
esp.phi_uniform_field_eV_per_e    均匀外场（若 --field）
esp.phi_eV_per_e                  三项之和
provenance                        GBW、密度名、orca_vpot 命令、原始文件
sampler                           四箱 requested/actual 点数、种子
```

`orca_vpot` 原始输出是 Bohr 坐标 + Hartree/e；JSON 统一转为 Å 和 eV/e（数值上等于 V）。
分子项已包含核与电子，**不要**再额外加一次 `Q/r`。

`--background` 的 JSON 是 `[[x_A, y_A, z_A, q_e], ...]`，与 ORCA `%pointcharges` 文件的
`q x y z` 列序不同。`--field` 单位是原子单位，原点 `(0,0,0)`。

这些选项只**添加外源势**；要包含电子对外场/背景电荷的响应，必须在 SCF 输入里同时写
`%pointcharges` 或 `%scf EField`，然后用该次 SCF 的 GBW 提取。

## 6. 验证

本仓库用 ORCA 6.1.0 容器、`wB97M-V/def2-TZVPD/AutoAux RIJCOSX`、水分子实际跑通：
SCF 正常结束，`orca_vpot water.densities` 列出 `water.scfp`，`orca_vpot` 对 30 个探针
（四箱 `[3,12,12,3]`）返回有限的分子 ESP。抽取脚本已用该 GBW/densities 端到端验证。

## 7. 常见问题

- **`orca_vpot` 找不到密度/segfault**：确认 `.gbw` 与 `.densities` 同目录，且不要手工指定
  不存在的 `.scfp`；让脚本通过 `orca_vpot *.densities` 自动发现。
- **容器内路径**：`--orca-container` 会自动 bind GBW 所在目录；GBW 与 `.densities` 必须同目录。
- **并行报 “ORCA has to be called with full pathname”**：用脚本，或显式调用完整路径。
- **带电体系**：`phi` 的零点不是无穷远真空；与 VASP 一样需明确参考约定。

## 8. 统一记录

```bash
cd ..
python3 unify.py orca --esp water_esp.json --out wb97m-v_def2-tzvpd.out \
  --engrad wb97m-v_def2-tzvpd.engrad --inp wb97m-v_def2-tzvpd.inp \
  --xyz water.xyz --output record.json
```

ORCA 分子没有应力张量，统一记录会写 `stress.available=false`；字段定义见
`../skill/esp-potential/references/unified-record.md`。
