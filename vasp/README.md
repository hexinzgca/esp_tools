# VASP：LOCPOT 电势 + 10×原子数格点抽取

## 1. `LVHAR` 到底输出什么

| INCAR | LOCPOT 内容 | 是否用于 ESP |
|---|---|---|
| `LVHAR = .TRUE.` | $V_{\rm ionic}(\mathbf r) + V_{\rm Hartree}(\mathbf r)$ | **是**，静电势 |
| `LVTOT = .TRUE.` | 再加 $V_{\rm XC}(\mathbf r)$ | 否，不是纯静电势 |

`LVHAR=.TRUE.` 会强制 `LVTOT=.FALSE.`。LOCPOT 的数值单位是 **eV**，但它是**电子势能**
（核附近为负，对正电荷是排斥势）。ORCA `orca_vpot` 输出的是静电势（核附近为正），因此
抽取程序统一做

$$
\phi_{\rm ESP} = -\frac{V_{\rm LOCPOT}}{e},\qquad e=1\ \text{(atomic units)}
$$

JSON 同时保留原始 `phi_locpot_eV`，便于核对。`ISPIN=2` 时 LOCPOT 含两个 block（spin up/down）；
`LVHAR` 下两块的 Hartree/ionic 部分应相同，抽取程序取平均并把最大差值写入
`spin_block_max_abs_diff_eV` 供检查。

## 2. INCAR 设置

复制 `INCAR.example`，核心只有：

```text
PREC   = Accurate
ENCUT  = 520
EDIFF  = 1E-6
ISMEAR = 0
SIGMA  = 0.05
IBRION = -1
NSW    = 0
ISPIN  = 1
ISYM   = 0
LREAL  = .FALSE.
LWAVE  = .FALSE.
LCHARG = .FALSE.
LVHAR  = .TRUE.
LVTOT  = .FALSE.
NCORE  = 4
KPAR   = 1
```

- `PREC=Accurate` 决定细 FFT 网格；LOCPOT 写在 `NGXF×NGYF×NGZF` 上。
- `LCHARG` / `LWAVE` 对 ESP 不是必需；关闭可显著减少磁盘占用。
- 孤立分子用 **Γ-only**（`KPOINTS.example`，`vasp_gam`）。周期性体系按常规 K 点。
- 盒子要留足够真空：水分子建议 ≥ 16 Å 立方盒；`>=5 Å` 远场 bin 才有足够体积和近似真空参考。
- 只对**中性**体系使用 `--reference vacuum`。带电体系/背景电荷（`LCHRED`）没有有限真空能级，
  用 `--reference none` 或显式数值，并明确参考零点。
- 若周期盒没有距所有原子 `>=5 Å` 的区域（小盒或致密固体），`[5,∞)` 的配额会自动并入
  `[1.5,5)`，总点数仍严格为 `10×N`，并打印 **WARNING**。JSON 同时记录
  `requested_bin_counts` / `actual_bin_counts` / `shortfall_bins` / `reallocated_bin_counts`。
  只有 `[1.5,5)` 等近场箱也无法填充时才报错。
- 这种情况下 `--reference vacuum` 没有远场点，会自动回退到 `--reference none`
  （LOCPOT 零点），并打印 WARNING；JSON 的 `locpot.reference.fallback=true` 记录该行为。

## 3. 运行

```bash
# 目录内准备 POSCAR / POTCAR / KPOINTS / INCAR
# POTCAR 顺序必须与 POSCAR 的元素顺序一致（示例：O 然后 H）
apptainer exec --bind "$PWD:/work" \
  /mnt/e/apptainer_sifs/dft1/vasp-cpu-intel.sif \
  bash -c 'cd /work && mpirun -np 4 /opt/vasp/bin/vasp_gam > vasp.out 2>&1'
```

本机镜像说明：

- `vasp-cpu-gcc.sif` 可直接运行 VASP，但其内置 POTCAR 库不完整；本仓库验证时把 Intel 镜像里的
  `O`/`H` POTCAR 复制出来，再用 GCC 镜像运行（Intel 镜像在 WSL 下 `vasp_gam` 多进程曾 SIGSEGV）。
- `vasp-cpu-intel.sif` 有完整 POTCAR 库（`/opt/vasp_pot/<Element>/POTCAR`），但本机 WSL 下运行
  `vasp_gam` 不稳定；若你的平台正常可直接用 Intel 镜像。
- 计算本身与抽取程序无关，LOCPOT 格式相同。

确认得到 `LOCPOT`，且其头部 `NGXF NGYF NGZF` 与预期一致。

## 4. 抽取：只算 10×N 个点

```bash
python vasp/extract_esp.py --locpot LOCPOT --output esp_vasp.json
# 可选：核对 POSCAR 与 LOCPOT 的晶格/坐标一致
python vasp/extract_esp.py --locpot LOCPOT --poscar POSCAR --output esp_vasp.json
```

默认参数：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--sample-mode` | `dft-grid` | 直接使用 LOCPOT 格点；`interpolate` 为任意点 + 三线性 |
| `--points-per-atom` | 10 | 总点数 = 10 × natoms，严格保证 |
| 分箱 | `1:4:4:1` | 最近原子距离 `[0.25,0.5) [0.5,1.5) [1.5,5) [5,∞)` Å |
| `--reference` | `vacuum` | 取 `>=5 Å` 远场点平均作为真空参考 |
| `--seed` | 20260916 | 与几何哈希混合，同结构可复现 |
| `--min-distance` | 0.25 Å | 不允许点落到核上 |

`dft-grid` 模式下每个探针都是 LOCPOT 的 `(i,j,k)` 格点，`phi_locpot_eV` 直接读取
`potential.ravel(order="C")[index]`，**没有插值误差**。JSON 增加
`esp.grid_indices_ijk`，`sampler.sampling_mode="dft-grid"`。网格分辨率由
`NGXF/NGYF/NGZF` 决定：`PREC=Accurate` 通常约 0.05–0.1 Å；若某距离箱格点不足，
仍按规则把远场缺口并入 `[1.5,5)` 并 warning。

`interpolate` 模式保留原来的任意点 + 周期性三线性，用于与旧结果对比。

脚本只处理这 10×N 个探针（默认直接读取 LOCPOT 格点，可选三线性插值），
不会导出完整 LOCPOT 网格。`--reference vacuum` 取远场 LOCPOT 平均后换算成静电势零点；
`--reference none` 以 LOCPOT 自身零点为参考；`--reference 5.2` 表示静电势零点为 5.2 eV
（即 $\phi=-V_{\rm LOCPOT}-5.2$）。

### 4.1 插值 vs 直接格点的数值测试

同一 8 Å 水盒、同一结构，分别用 `NGXF=96³` 与 `NGXF=192³` 计算，取 192³ 的 30 个
直接格点作为参考，比较 96³ 三线性插值在这些点上的误差：

| 距离箱 | 点数 | 平均绝对误差 / eV | 最大绝对误差 / eV |
|---|---:|---:|---:|
| `[0.25,0.5)` | 3 | 0.563 | 0.846 |
| `[0.5,1.5)` | 12 | 0.0276 | 0.0948 |
| `[1.5,5)` | 12 | 0.000160 | 0.00131 |
| `[5,∞)` | 3 | 0.000008 | 0.000022 |
| 全部 | 30 | 0.0675 | 0.846 |

同时验证：在自身格点上三线性插值与直接值完全一致（最大差 0.0）。因此：

- 插值误差主要集中在 `[0.25,0.5)` 近核箱，那里势的曲率最大；`[1.5,5)` 与远场
  误差已很小。
- 默认 `dft-grid` **完全消除插值误差**；剩余的是 FFT 格点离散误差（96³ 对 192³
  的直接格点比较：全部平均 0.0098 eV、最大 0.223 eV，同样集中在近核箱）。
- 需要严格可复现的 DFT 格点值时用 `dft-grid`；需要任意空间点或与旧数据对齐时
  才用 `interpolate`。

复现该对比：

```bash
python vasp/accuracy_check.py \
  --coarse-locpot coarse/LOCPOT \
  --fine-locpot fine/LOCPOT \
  --fine-esp fine/grid.json \
  --output accuracy.json
```

`fine/grid.json` 必须由 `--sample-mode dft-grid` 生成；脚本会同时检查“细网格
自身格点上插值 = 直接值”。

输出 JSON 关键字段：

```text
structure           晶格、坐标、元素
sampler             sampling_mode、requested/actual 四箱点数、种子、分箱参数
                    shortfall_bins、reallocated_bin_counts、bin_reallocation
locpot              grid_shape、block 数、spin block 差、参考值/回退、
                    interpolation、符号约定
esp.positions_Angstrom
esp.grid_indices_ijk  dft-grid 模式下的 LOCPOT 格点索引
esp.phi_locpot_eV   LOCPOT 原始电子势能（格点直接值或插值）
esp.phi_eV_per_e    静电势 = -LOCPOT - reference
```

## 5. 采样分布与近核问题

参考实现把 10% 点放在 `0.25–0.5 Å`、40% 在 `0.5–1.5 Å`、40% 在 `1.5–5 Å`、10% 在 `>=5 Å`。
这已经避免旧方案中大量点贴近核的分布；所有点距最近原子 **≥0.25 Å**。

VASP 的 PAW 势在核附近是赝势势，不是全电子 ESP；对多数元素，`0.25–0.5 Å` 可能落在
augmentation sphere 内。因此：

- 该 bin 适合与同一 PAW 设置下的数据做机器学习一致性训练；
- 若目标是与全电子 ESP（如 ORCA `orca_vpot`）逐点对比，应把第一 bin 排除或仅使用 `>=0.5 Å` 的点。
  正常分子盒中保留 1:4:4:1；只有远场 `[5,∞)` 不可用时才按 WARNING 把缺口并入
  `[1.5,5)`，不会静默改比例或静默删点。

## 6. 验证

本仓库用 VASP 6.5.1 + PBE PAW（O/H）、16 Å 水盒、Γ-only 实际生成 `LOCPOT`
（`240×240×240`，240 MB），抽取 30 点（10×3）：

- 四箱实际点数 `[3, 12, 12, 3]`；
- 远场（`>=5 Å`）LOCPOT 平均 `0.0174 eV`，换算成静电势零点 `-0.0174 eV`；
- 静电势范围 `-0.497 … 53.52 eV`（近核 bin 最大，符合电子势能取负后的预期）；
- `ISPIN=1`，LOCPOT block 数 1。

合成 LOCPOT 测试（线性势 + 单点 delta）验证了 VASP 的 x-fastest 数据顺序、双 block 平均、
三线性插值精度，以及 `dft-grid` 输出点确实是对应 `(i,j,k)` 格点且值直接读取。
8 Å 水盒 96³/192³ 的插值误差实验见 §4.1：96³ 插值全部平均 `0.0675 eV`、最大 `0.846 eV`；
`[1.5,5)` 与远场分别 ≤`0.0013`、`0.000022 eV`。

另外用真实 VASP 跑了一个 **4 Å 单 O 盒**（不存在 `>=5 Å` 远场）验证回退逻辑：输出
`requested [1,4,4,1] → actual [1,4,5,0]`，`reallocated_bin_counts [0,0,1,0]`，
`--reference vacuum` 自动回退为 `none`，程序正常退出并打印两条 WARNING。

## 7. 常见问题

- **`LOCPOT` 不存在**：`LVHAR` 未开，或计算未正常结束。
- **参考值随盒子变化**：真空不足或体系带净电荷；加大盒子/改用显式参考。
- **出现 `bin_reallocation` / `reallocated_bin_counts`**：盒子没有 `>=5 Å` 远场；
  远场配额已并入 `[1.5,5)`，总点数不变。若 `reference.fallback=true`，真空参考也回退到了
  LOCPOT 零点。
- **近核值异常**：PAW 赝势势的正常表现；见第 5 节。
- **`spin_block_max_abs_diff_eV` 很大**：可能误用了 `LVTOT`，或 LOCPOT 不是 `LVHAR` 输出。

## 8. 统一记录

```bash
cd ..
python3 unify.py vasp --esp esp_vasp.json --outcar OUTCAR \
  --incar INCAR --kpoints KPOINTS --potcar POTCAR --poscar POSCAR \
  --output record.json
```

字段定义见 `skill/esp-potential/references/unified-record.md`。
