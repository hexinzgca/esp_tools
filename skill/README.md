# esp-potential skill

本目录按 omp / Agent Skills 的约定部署：

```text
skill/
├── README.md                 # 本文件：注册说明
└── esp-potential/
    ├── SKILL.md              # 主说明（frontmatter + 工作流）
    ├── references/
    │   ├── vasp.md           # VASP 输入/输出/抽取细节
    │   ├── orca.md           # ORCA 基组/输入/orca_vpot 细节
    │   └── unified-record.md # 统一 JSON 字段定义
    └── scripts/
        └── check_update.sh   # 远程更新预检
```

## 注册

把 `skill/` 目录加入 omp 的 custom skills 目录，例如：

```yaml
skills:
  customDirectories:
    - /absolute/path/to/esp_tools/skill
```

重启/重新加载后，skill 名为 `esp-potential`：

- 模型按需读取：`skill://esp-potential`
- 交互式调用：`/skill:esp-potential`
- 子文件：`skill://esp-potential/references/vasp.md`

目录扫描是**非递归**的，因此必须是 `skill/<skill-name>/SKILL.md`；不要把 `SKILL.md` 直接放在 `skill/` 下。

## 脚本位置

skill 只包含说明，不复制脚本。脚本在本仓库中：

```text
esp_tools/vasp/extract_esp.py
esp_tools/orca/extract_esp.py
esp_tools/orca/run_orca_esp.sh
esp_tools/unify.py
```

从 `esp_tools/` 目录运行。若把 skill 单独复制到别处，请同时复制脚本，或修改说明中的路径。

## 更新预检

运行 skill 前先检查远程仓库：

```bash
bash skill/esp-potential/scripts/check_update.sh /path/to/esp_tools
```

脚本只执行 `git fetch` 并比较本地 `HEAD` 与 `origin/<branch>`；落后时提示
`git pull --ff-only`，不会自动修改工作区。

## 验证

- 16 个单元测试：`python3 -m unittest discover -s tests`
- 真实 ORCA 6.1.0 与 VASP 6.5.1 流程见主 README 的“验证状态”。
