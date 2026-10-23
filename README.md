# esp_tools

VASP / ORCA 静电势（ESP）格点采样与抽取工具。两个分支共用同一套采样逻辑：

- **VASP**：`LVHAR=.TRUE.` 生成 `LOCPOT`，在 `10 × natoms` 个探针点上做周期性三线性插值；
- **ORCA**：`wB97M-V` + 弥散 def2 基组做 SCF，用 `orca_vpot` 在同一套探针点上评估分子 ESP。

只抽取探针点，**不导出完整 LOCPOT 网格，也不做全空间 cube**。

## 目录

```text
esp_tools/
├── common/sampling.py       共享采样：10×N、1:4:4:1、几何哈希种子
├── vasp/
│   ├── README.md            INCAR/LOCPOT/参考零点/PAW 注意事项
│   ├── INCAR.example        LVHAR 静电势设置
│   ├── KPOINTS.example      孤立分子 Γ-only
│   ├── POSCAR.example       16 Å 水盒
│   ├── locpot.py            LOCPOT 解析 + 周期性三线性插值
│   └── extract_esp.py       CLI：LOCPOT → esp JSON
├── orca/
│   ├── README.md            基组选择、SCF 输入、orca_vpot 用法
│   ├── wb97m-v_def2-tzvpd.inp   默认输入（推荐）
│   ├── wb97m-v_def2-svpd.inp    双-zeta 弥散备选
│   ├── water.xyz
│   ├── run_orca_esp.sh      容器内跑 SCF + 抽取
│   └── extract_esp.py       CLI：XYZ/points + GBW → esp JSON
├── unify.py                 统一记录：结构/方法/能量/力/应力/电势
├── schema/
│   └── record-v1.schema.json 统一记录 JSON Schema
├── skill/
│   ├── README.md            skill 注册说明
│   └── esp-potential/
│       ├── SKILL.md         主 skill
│       ├── references/      vasp / orca / unified-record 细节
│       └── scripts/check_update.sh  远程更新预检
├── tests/                   采样、LOCPOT、插值、unify 测试
└── requirements.txt
```

## 安装

```bash
python -m pip install -r requirements.txt   # 仅需 numpy
```

要求 Python ≥ 3.10。VASP/ORCA 本身通过 Apptainer 镜像调用，不纳入 Python 依赖。

## 采样规则

默认 **每原子 10 点**，总数严格等于 `10 × natoms`。按到最近原子的距离分四箱：

| bin | 距离 / Å | 目标比例 | 水分子（3 原子） |
|---|---|---|---|
| 0 | `[0.25, 0.5)` | 1 | 3 |
| 1 | `[0.5, 1.5)` | 4 | 12 |
| 2 | `[1.5, 5.0)` | 4 | 12 |
| 3 | `[5.0, ∞)` | 1 | 3 |

- 第一箱从 **0.25 Å** 起，点不会落在核上；相比旧版按 DFT 积分壳层取点、大量点贴近核的分布，
  本方案把 80% 点放在 `0.5–5 Å` 的化学相关区域，仅 10% 在最内层。
- 种子 = `20260916` 与几何哈希混合；同结构完全可复现，不同结构不共享采样模式。
- 非周期（ORCA）：远场点在质心球壳 `r_max+5 … r_max+10 Å`。
- 周期（VASP）：所有点折回晶胞，距离用最小镜像；远场点在晶胞内均匀抽样并要求 `>=5 Å`。
  因此 VASP 盒子必须留足够真空（水建议 16 Å 立方盒）。
- VASP 周期盒中若 `[5,∞)` 无法填充，自动把缺口并入 `[1.5,5)`，总点数仍为 `10×N`，
  并打印 WARNING；JSON 记录 `requested_bin_counts` / `actual_bin_counts` /
  `shortfall_bins` / `reallocated_bin_counts`。只有近场箱也无法填充时才报错。
- ORCA 的远场球壳总能填充，默认严格按 1:4:4:1 执行，不启用重分配。

采样分布取自参考 `esp_grid_share/esp_grid.py` 的四箱方案；本仓库只把每原子点数从
“20 向上取整到 14 的倍数”改为**精确 10**，并增加 VASP 周期最小镜像模式。
VASP 侧沿用 `dft_additional` 的 `LVHAR=.TRUE.` 输出约定。

## VASP 快速开始

```bash
cd vasp
cp INCAR.example INCAR
cp KPOINTS.example KPOINTS
cp POSCAR.example POSCAR
# 按 POSCAR 元素顺序拼接 POTCAR（示例：O 然后 H）
apptainer exec --bind "$PWD:/work" \
  /mnt/e/apptainer_sifs/dft1/vasp-cpu-intel.sif \
  bash -c 'cd /work && mpirun -np 4 /opt/vasp/bin/vasp_gam > vasp.out 2>&1'

python extract_esp.py --locpot LOCPOT --poscar POSCAR --output esp_vasp.json
```

关键点：`LVHAR=.TRUE.` 写的是 `V_ionic + V_Hartree`（单位 eV），不是含 XC 的 `LVTOT`。
LOCPOT 是**电子势能**，抽取程序取负号得到与 ORCA 一致的静电势。默认参考零点取 `>=5 Å`
远场点的平均势（中性分子 + 大盒子）；可用 `--reference none` / `--reference <eV>` 改。
详见 `vasp/README.md`。

## ORCA 快速开始

```bash
cd orca
./run_orca_esp.sh wb97m-v_def2-tzvpd.inp water.xyz
# 或分步：
apptainer exec --bind "$PWD:/work" \
  /mnt/e/apptainer_sifs/build/orca-6.1.sif \
  bash -c 'cd /work; ORCA=$(command -v orca); exec "$ORCA" wb97m-v_def2-tzvpd.inp' \
  > wb97m-v_def2-tzvpd.out 2>&1
python extract_esp.py --xyz water.xyz --gbw wb97m-v_def2-tzvpd.gbw \
  --orca-container /mnt/e/apptainer_sifs/build/orca-6.1.sif \
  --output water_esp.json
```

ORCA 默认 `wB97M-V/def2-TZVPD/AutoAux RIJCOSX`（三-zeta 弥散，OMol25 同级别）；
若需要双-zeta 弥散，使用 `wb97m-v_def2-svpd.inp`。详见 `orca/README.md`。

## 输出

两个 CLI 都输出 JSON，单位统一为 Å 和 eV/e（数值等于 V）：

- `esp.positions_Angstrom`：探针坐标；
- VASP：`esp.phi_locpot_eV`（原始电子势能）、`esp.phi_eV_per_e`
  （静电势 `-LOCPOT - reference`）；
- ORCA：`phi_molecular_eV_per_e`、`phi_background_eV_per_e`、
  `phi_uniform_field_eV_per_e`、`phi_eV_per_e`；
- `sampler`：分箱、点数、种子、是否周期；
- `provenance` / `locpot`：输入文件哈希、网格尺寸、密度名、`orca_vpot` 命令、日志路径。

输出文件默认拒绝覆盖；需要覆盖时加 `--force`。

## 统一记录：结构 + 方法 + 能量/力/应力 + 电势

建议每个结构最终产出一条自包含记录（`esp_tools/record-v1`），避免后续再解析
OUTCAR/ORCA 输出：

```bash
# VASP：ESP JSON + OUTCAR/INCAR/KPOINTS/POTCAR/POSCAR
python3 unify.py vasp --esp esp_vasp.json --outcar OUTCAR \
  --incar INCAR --kpoints KPOINTS --potcar POTCAR --poscar POSCAR \
  --output record.json

# ORCA：ESP JSON + out/engrad/inp/xyz
python3 unify.py orca --esp water_esp.json --out water.out \
  --engrad water.engrad --inp water.inp --xyz water.xyz \
  --output record.json
```

记录包含：

- `structure`：元素、坐标、晶格、pbc、电荷、多重度/自旋；
- `method`：泛函、基组、辅助基组、软件版本、VASP POTCAR 赝势（TITEL/ZVAL）、
  INCAR / ORCA 输入设置；
- `energy`：eV、Hartree、每原子 eV；
- `forces`：eV/Å（ORCA 由 `-gradient(Eh/bohr) × 51.422067` 转换）；
- `stress`：VASP 同时给 kBar 与 eV/Å³；ORCA 非周期分子明确标记不可用；
- `electrostatic_potential.points`：`[x, y, z, q]`，q 为总电势 eV/e；
- `method.software`：VASP/ORCA 版本、build/tag、libXC、ORCA 容器路径；
- `method.pseudopotentials`：逐元素 POTCAR `TITEL`（含日期/版本）和 `ZVAL`，
  以及 `source.files.potcar` 的 SHA-256；
- `provenance.esp_tools`：`unify.py` 所在仓库的 git commit 与 dirty 状态。

详见 `skill/esp-potential/references/unified-record.md`。

## Skill

`skill/esp-potential/SKILL.md` 是可直接注册的 omp skill，包含 VASP/ORCA
输入设置、抽取命令、统一 JSON 字段和 QA 清单。注册方式见 `skill/README.md`：
把 `esp_tools/skill` 加入 `skills.customDirectories`，使用
`skill://esp-potential` 或 `/skill:esp-potential`。

skill 预检会先检查远程更新：

```bash
bash skill/esp-potential/scripts/check_update.sh /path/to/esp_tools
```

脚本只执行 `git fetch` 并比较 `HEAD` 与 `origin/<branch>`；落后时提示
`git pull --ff-only`。

## 验证状态

- ORCA 6.1.0（Apptainer）：水分子 `wB97M-V/def2-TZVPD/AutoAux RIJCOSX` SCF 正常，
  `orca_vpot` 通过 `.densities` 自动发现 `water.scfp`，30 点（`[3,12,12,3]`）端到端抽取成功。
- VASP 6.5.1（Apptainer）：水分子 16 Å 盒、Γ-only、PBE PAW 生成 `240×240×240` 的 `LOCPOT`；
  抽取 30 点，四箱 `[3,12,12,3]`，远场参考约 `0.017 eV`。合成 LOCPOT 测试验证了 VASP 的
  x-fastest 存储顺序、双 block 平均与三线性插值。
- VASP 致密盒（4 Å 单 O）：无 `>=5 Å` 远场，`requested [1,4,4,1] → actual [1,4,5,0]`，
  远场缺口自动并入 `[1.5,5)`，参考自动回退到 `none`，正常退出并打印 WARNING。
- `unify.py` 已用真实 ORCA 水分子 `.out/.engrad/.inp` 和 VASP 16 原子 diamond
  `OUTCAR` 验证能量、力、应力、统一记录写出。
- `tests/` 可在无 VASP/ORCA 的情况下运行（16 个测试，含统一记录解析）。

## 已知限制

- **VASP 近核**：PAW 的 `0.25–0.5 Å` bin 是赝势势，不是全电子 ESP；与 ORCA 全电子值逐点比较时
  应排除该 bin。工具保留 1:4:4:1 比例，不静默删点。
- **VASP 参考零点**：周期 LOCPOT 只定义到常数；`--reference vacuum` 依赖足够真空与中性体系。
  带电体系/背景电荷没有有限真空能级，需显式参考。
- **ORCA 外场/背景电荷**：`--background` / `--field` 只加外源势；电子响应必须来自包含
  `%pointcharges` / `%scf EField` 的真实 SCF。
- **`orca_vpot` 密度容器**：`.gbw` 与 `.densities` 必须同目录；旧版本独立 `.scfp` 自动回退。
- **致密固体/小盒子**：没有距所有原子 `>=5 Å` 的远场区域时，VASP 会把远场配额自动并入
  `[1.5,5)` 并打印 WARNING，总点数不变；`--reference vacuum` 会自动回退到 `none`。
  只有近场箱也无法填充时才报错。
