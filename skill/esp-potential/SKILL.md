---
name: esp-potential
description: Extract DFT electrostatic potentials on an exact 10 x natoms probe grid with VASP (LVHAR/LOCPOT) or ORCA (wB97M-V/def2-TZVPD + orca_vpot), then emit one JSON record with structure, method, energy, forces, stress, and [x,y,z,phi] points.
globs:
  - "**/INCAR"
  - "**/LOCPOT"
  - "**/*.inp"
  - "**/*.gbw"
  - "**/*.engrad"
---

# ESP potential extraction with VASP / ORCA

## When to use

Use this skill when you need one or more of:

- an electrostatic potential (ESP) on a **small probe set**, not the full FFT grid;
- a dataset record for ML/analysis containing structure + method + energy + forces +
  stress + ESP;
- to check that a VASP `LOCPOT` or ORCA `orca_vpot` workflow is physically consistent
  (units, sign, reference, spin/charge).

Scripts live in the `esp_tools` repository. Run commands from the repository root:

```text
vasp/extract_esp.py       # LOCPOT -> ESP JSON
orca/extract_esp.py       # XYZ + GBW -> ESP JSON
orca/run_orca_esp.sh      # ORCA SCF in Apptainer + extraction
unify.py                  # ESP JSON + OUTCAR/ORCA outputs -> unified record
```

## Preflight: remote update check

Before extracting, check the checkout:

```bash
bash skill/esp-potential/scripts/check_update.sh /path/to/esp_tools
```

The script only runs `git fetch` and compares `HEAD` with `origin/<branch>`.
If it reports `behind`, update first:

```bash
git -C /path/to/esp_tools pull --ff-only
```

If the working tree is dirty, do not pull over local changes; reconcile them
first. Do not run extraction with a stale checkout.

## Non-negotiable rules

| Item | Rule |
|---|---|
| Probe count | exactly `10 × natoms` |
| Distance bins | `[0.25,0.5) [0.5,1.5) [1.5,5) [5,∞)` Å |
| Target ratio | `1 : 4 : 4 : 1` |
| Minimum distance | 0.25 Å; never sample a nucleus |
| VASP distances | minimum-image in the periodic cell |
| ORCA distances | Euclidean, non-periodic |
| VASP far shortfall | merge `[5,∞)` quota into `[1.5,5)` first; warn, do not abort |
| Full grid | never export the complete LOCPOT/FFT grid |
| VASP potential sign | `LOCPOT` is electron potential energy; ESP `= -LOCPOT` |
| ORCA density | `.gbw` **and** matching `.densities` are required |

The sampler uses a geometry-hashed seed: the same structure reproduces the same
probe set, different structures do not share the same pattern.

## VASP: input settings and extraction

### Required physics

`LVHAR = .TRUE.` writes $V_{\rm ionic}+V_{\rm Hartree}$ to `LOCPOT` in eV. This is
the electrostatic potential **energy of an electron**, negative near nuclei. The
tool converts it to the electrostatic potential:

```text
phi_eV_per_e = -LOCPOT_eV - reference_eV
```

Do **not** use `LVTOT = .TRUE.` for ESP: it adds the XC potential. Do not compare
raw `LOCPOT` values directly with ORCA `orca_vpot` without the sign flip.

### Minimal INCAR

```text
PREC   = Accurate
ENCUT  = 520
EDIFF  = 1E-6
NELM   = 200
NELMIN = 4
ISMEAR = 0
SIGMA  = 0.05
ISPIN  = 1
ISYM   = 0
LREAL  = .FALSE.
IBRION = -1
NSW    = 0
LWAVE  = .FALSE.
LCHARG = .FALSE.
LVHAR  = .TRUE.
LVTOT  = .FALSE.
NCORE  = 4
KPAR   = 1
```

- `PREC=Accurate` controls the fine FFT grid used by `LOCPOT`.
- `LCHARG`/`LWAVE` are not needed for ESP; disable them to save disk.
- Isolated molecule: Γ-only (`KPOINTS` Gamma, `vasp_gam`). Use a large box
  (16 Å cubic for water) so the `>=5 Å` far bin exists.
- Neutral systems only for `--reference vacuum`. Charged/background systems have
  no finite vacuum level; use `--reference none` or an explicit value.
- If the cell has no `>=5 Å` region, the far quota is **reallocated to `[1.5,5)`**
  with a WARNING; the total remains `10×N`. `--reference vacuum` then falls back
  to `LOCPOT` zero and records `reference.fallback=true`.
- POTCAR order must match POSCAR species order.

### Run and extract

```bash
cd esp_tools
python3 vasp/extract_esp.py --locpot LOCPOT --poscar POSCAR --output esp_vasp.json
```

Useful options:

| Option | Default | Meaning |
|---|---|---|
| `--points-per-atom` | 10 | total points = value × natoms |
| `--reference` | `vacuum` | mean of `>=5 Å` points, `none`, or explicit eV |
| `--seed` | 20260916 | geometry-hashed base seed |
| `--force` | off | overwrite output |

The JSON contains `structure`, `sampler`, `locpot` (grid shape, blocks, reference,
sign convention), and `esp`:

```text
esp.positions_Angstrom
esp.phi_locpot_eV      raw electron potential energy from LOCPOT
esp.phi_eV_per_e       ESP = -LOCPOT - reference
```

`sampler.actual_bin_counts` may differ from `requested_bin_counts` only when
`sampler.bin_reallocation=true`; check `shortfall_bins` and
`reallocated_bin_counts`.

## ORCA: input settings and extraction

### Basis

Default to `def2-TZVPD` (triple-zeta + diffuse, OMol25 level); use `def2-SVPD`
for a double-zeta + diffuse alternative.

| Keyword | Level | Note |
|---|---|---|
| `def2-SVPD` | double-zeta + diffuse | alternative |
| `def2-TZVPD` | triple-zeta + diffuse | **default** |
| `def2-TZVPPD` | triple-zeta, double polarization + diffuse | higher accuracy |
| `def2-QZVPD` | quadruple-zeta + diffuse | reference quality |

### Minimal input

```text
! wB97M-V def2-TZVPD AutoAux RIJCOSX VeryTightSCF

%pal nprocs 4 end
%maxcore 4000

%method
  Grid 5
end

* xyzfile 0 1 water.xyz
```

- `AutoAux` generates the auxiliary basis automatically; `RIJCOSX` accelerates
  exchange. `VeryTightSCF` is recommended because ESP is sensitive to density
  convergence.
- ORCA 6.1 writes `water.gbw` and `water.densities`; `orca_vpot` needs both.
- When `%pal nprocs > 1`, ORCA must be launched with its full pathname. The
  provided `run_orca_esp.sh` handles this.
- For anions/cations, change the charge/multiplicity in the `* xyzfile` line and
  rerun SCF. The probe grid can be reused (`--points`) only for the same geometry.

### Run and extract

```bash
cd esp_tools/orca
./run_orca_esp.sh wb97m-v_def2-tzvpd.inp water.xyz
```

Or step by step:

```bash
cd esp_tools
python3 orca/extract_esp.py --xyz orca/water.xyz --gbw water.gbw \
  --orca-container /mnt/e/apptainer_sifs/build/orca-6.1.sif \
  --output water_esp.json
```

Native ORCA:

```bash
python3 orca/extract_esp.py --xyz orca/water.xyz --gbw water.gbw \
  --orca-home /opt/orca/orca-6.1.0-f.0_linux_x86-64/bin \
  --output water_esp.json
```

Sample only (no GBW), or reuse points:

```bash
python3 orca/extract_esp.py --xyz orca/water.xyz --output water_grid.json
python3 orca/extract_esp.py --points water_grid.json --gbw anion.gbw \
  --orca-container /path/orca-6.1.sif --output anion_esp.json
```

`--background` is a JSON array `[[x_A,y_A,z_A,q_e], ...]`; `--field` is a uniform
field in atomic units with origin `(0,0,0)`. These options only **add external
source potentials**; the electronic response must come from a SCF that actually
included `%pointcharges` / `%scf EField`.

The ORCA JSON contains `phi_molecular_eV_per_e`, `phi_background_eV_per_e`,
`phi_uniform_field_eV_per_e`, and `phi_eV_per_e` (sum). ORCA references the
molecular potential to infinity.

## Unified JSON record

Yes: for a dataset, emit **one record per structure** with structure, method,
energy, forces, stress, and the ESP points. Use `unify.py`:

```bash
cd esp_tools

# VASP
python3 unify.py vasp --esp esp_vasp.json --outcar OUTCAR \
  --incar INCAR --kpoints KPOINTS --potcar POTCAR --poscar POSCAR \
  --output record.json

# ORCA
python3 unify.py orca --esp water_esp.json --out water.out \
  --engrad water.engrad --inp water.inp --xyz water.xyz \
  --output record.json
```

The format is `esp_tools/record-v1`; the machine-readable schema is
`schema/record-v1.schema.json` in the repository.

```text
{
  "format": "esp_tools/record-v1",
  "source":     {code, files{path,sha256,bytes}},
  "structure":  {natoms, symbols, counts, formula, positions_Angstrom,
                 lattice_Angstrom, pbc, charge, multiplicity, spin},
  "method":     {code, functional, basis, aux_basis, software,
                 pseudopotentials, approximations, solvent, settings, ...},
  "energy":     {total_eV, total_hartree, per_atom_eV, unit, source},
  "forces":     {values_eV_per_Angstrom, max_norm..., unit, source},
  "stress":     {available, values_kBar, values_eV_per_Angstrom3,
                 pressure_kBar, sign_convention} | {available:false, reason},
  "electrostatic_potential": {
      points: [[x,y,z,q], ...],
      point_columns: ["x_Angstrom","y_Angstrom","z_Angstrom","phi_eV_per_e"],
      units, reference, sampler, components_eV_per_e, source_report},
  "units":      {...},
  "provenance": {tool, esp_tools{commit,dirty}, ...}
}
```

- `method.software`: VASP/ORCA version, build/tag/variant, libXC version, and
  ORCA container path when available.
- `method.pseudopotentials`: per-element POTCAR `TITEL` (includes the
  pseudopotential date/version) and `ZVAL`. The POTCAR file SHA-256 is in
  `source.files.potcar`.
- `provenance.esp_tools.commit/dirty`: the extraction tool revision.
- `[x,y,z,q]`: `x/y/z` in Å; `q = phi_eV_per_e` is the **total** ESP in eV/e
  (numerically V). For ORCA this includes molecular + background + field.
- VASP stress is stored both as printed (`kBar`) and converted
  (`eV/Å³`). Sign convention is recorded; do not silently flip it.
- ORCA non-periodic molecules have no stress tensor; `stress.available=false`
  with a reason. Do not fabricate zeros.
- ORCA forces are `-gradient(Eh/bohr) × 51.422067` eV/Å.
- Energy is stored in eV and Hartree, plus per-atom eV.

For a JSONL dataset, write one `record.json` per structure and concatenate the
JSON objects one per line; keep the original ESP report for provenance.

## QA checklist before trusting a record

0. `check_update.sh` reports up to date (otherwise update and rerun).
1. `npoints == 10 × natoms`.
2. `sampler.requested_bin_counts` is `[N,4N,4N,N]`; if not, inspect
   `bin_reallocation`, `shortfall_bins`, `reallocated_bin_counts`.
3. VASP: `LVHAR=.TRUE.`, `LOCPOT` exists, `locpot.n_blocks` and grid shape sane.
4. VASP: `phi` is positive near nuclei after the sign flip; far-field `phi ≈ 0`
   when using `--reference vacuum`.
5. ORCA: SCF terminated normally; `.gbw` and `.densities` are in the same
   directory; `orca_vpot` found exactly one `.scfp`.
6. All energy/force/stress/phi values finite.
7. Forces rows == natoms; stress units explicit; charge/multiplicity consistent
   with the input.
8. No stale output: extraction refuses overwrite unless `--force`.

## Troubleshooting

| Symptom | Action |
|---|---|
| ORCA `has to be called with full pathname` | use `run_orca_esp.sh` or full path |
| `orca_vpot` cannot find density / segfault | ensure `.gbw` + `.densities` same dir; let the script discover `.scfp` via `orca_vpot *.densities` |
| VASP `vasp_gam` SIGSEGV in WSL | try the GCC image or fewer ranks; computation is independent of extraction |
| VASP far bin short | expected in dense/small cells; check WARNING, `actual_bin_counts`, and `reference.fallback` |
| `LOCPOT` missing | `LVHAR` not enabled or SCF did not finish |
| VASP vs ORCA values differ strongly near nuclei | PAW pseudo-potential vs all-electron; exclude the `0.25–0.5 Å` bin for cross-code comparison |
| stress suspicious | read `sign_convention`; VASP prints external-pressure sign, not a continuum-mechanics tensile convention |

## References

- `skill://esp-potential/references/vasp.md` — INCAR/POSCAR/LOCPOT details.
- `skill://esp-potential/references/orca.md` — basis, input blocks, `orca_vpot`.
- `skill://esp-potential/references/unified-record.md` — full record schema.
