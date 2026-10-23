# Unified record `esp_tools/record-v1`

Machine-readable schema: `../../schema/record-v1.schema.json` in the repository.

## Why one record

A dataset entry should be self-contained: structure, method, energy, forces,
stress, and ESP points in one JSON. This removes sidecar parsing and makes the
units/sign/reference explicit. Build it with `unify.py` from an existing ESP
report plus the code's native output.

```bash
cd esp_tools

python3 unify.py vasp --esp esp_vasp.json --outcar OUTCAR \
  --incar INCAR --kpoints KPOINTS --potcar POTCAR --poscar POSCAR \
  --output record.json

python3 unify.py orca --esp water_esp.json --out water.out \
  --engrad water.engrad --inp water.inp --xyz water.xyz \
  --output record.json
```

`unify.py` refuses to overwrite unless `--force`.

## Top-level fields

| Field | Meaning |
|---|---|
| `format` | always `esp_tools/record-v1` |
| `source.code` | `vasp` or `orca` |
| `source.files` | input/output paths, sizes, SHA-256 |
| `structure` | geometry, cell, charge, spin |
| `method` | functional, basis, auxiliary basis, settings |
| `energy` | total/per-atom energy |
| `forces` | Cartesian forces |
| `stress` | stress tensor or explicit unavailability |
| `electrostatic_potential` | probe points and potential |
| `units` | unit declaration |
| `provenance` | tool, timestamp, arguments |

## `structure`

| Field | Unit/type |
|---|---|
| `natoms` | int |
| `symbols` | list[str] |
| `counts` | list[int], same order as symbols |
| `formula` | string or null |
| `positions_Angstrom` | `natoms × 3` |
| `lattice_Angstrom` | `3 × 3` or null |
| `pbc` | 3 booleans |
| `charge` | float or null |
| `multiplicity` | int or null |
| `spin` | VASP: `{ispin, magmom}`; ORCA: `{multiplicity}` |

VASP charge is computed only when both POTCAR `ZVAL` and OUTCAR `NELECT` are
available:

```text
charge = sum(count_i * ZVAL_i) - NELECT
```

Otherwise it is `null`. ORCA charge/multiplicity come from the `.inp`
(`* xyz/xyzfile charge mult`) or, if absent, from the `.out`.

## `method`

VASP:

```text
code: "vasp"
functional: from POTCAR TITEL (PAW_PBE -> PBE) or INCAR GGA
basis: "PAW"
software: {name: "vasp", version, tag, build, variant}   # from OUTCAR header
pseudopotentials: [{symbol, title, zval}, ...]          # TITEL includes the date
potcar_titles: [...]
settings: parsed INCAR key/value strings
kpoints: raw KPOINTS text
```

ORCA:

```text
code: "orca"
functional: first non-keyword token of the ! line (e.g. wB97M-V)
basis: def2-TZVPD
aux_basis: AutoAux
software: {name: "orca", version, libxc_version, container}
approximations: [RIJCOSX, VeryTightSCF, ...]
solvent: CPCM or null
settings: simple input tokens, nprocs, maxcore, grid, efield,
          point_charges_file, cpcm, scf_converged, terminated_normally
```

Parsing is intentionally conservative. If a custom input uses an unusual
functional/basis token, inspect `method.settings.simple_input` rather than
trusting the heuristic.

## `energy`

```text
total_eV            VASP TOTEN (last ionic step) or ORCA final energy × 27.211386245988
total_hartree       total_eV / 27.211386245988
per_atom_eV         total_eV / natoms
unit                "eV"
source              human-readable source
```

Missing energy is `null`, never zero.

## `forces`

```text
values_eV_per_Angstrom   natoms × 3
max_norm_eV_per_Angstrom max row norm
unit                     "eV/Angstrom"
source                   VASP OUTCAR or ORCA engrad conversion
```

ORCA `.engrad` contains the gradient in Eh/bohr:

```text
force_eV_per_A = -gradient_Eh_per_bohr * 51.422067476
```

The negative sign is part of the force definition; do not drop it.

## `stress`

VASP:

```text
available: true
values_kBar: 3×3
values_eV_per_Angstrom3: 3×3
pressure_kBar: external pressure (OUTCAR; fallback trace/3)
unit: "kBar and eV/Angstrom^3"
sign_convention: "VASP OUTCAR values as printed; positive trace corresponds
                  to external pressure in the OUTCAR convention"
```

Conversion:

```text
1 kBar = 0.1 GPa
1 eV/A^3 = 160.21766208 GPa
eV/A^3 = kBar * 0.1 / 160.21766208
```

The six printed components are `xx yy zz xy yz zx`; the matrix is

```text
[[xx, xy, zx],
 [xy, yy, yz],
 [zx, yz, zz]]
```

ORCA:

```text
available: false
reason: "ORCA non-periodic molecular calculations do not define a stress tensor"
```

Do not fabricate zeros for missing ORCA stress.

## `electrostatic_potential`

| Field | Meaning |
|---|---|
| `points` | list of `[x, y, z, q]` |
| `point_columns` | `["x_Angstrom","y_Angstrom","z_Angstrom","phi_eV_per_e"]` |
| `point_format` | human-readable description |
| `units` | `{positions: Angstrom, phi: eV/e (= V)}` |
| `reference` | VASP vacuum/none/value or ORCA infinity |
| `sampler` | requested/actual bins, seed, reallocation fields |
| `components_eV_per_e` | VASP raw LOCPOT + total; ORCA molecular/background/field/total |
| `source_report` | original ESP JSON path, format, SHA-256 |

`q` is the **total** electrostatic potential in eV/e, numerically volt. For
ORCA it is `molecular + background + uniform field`. For VASP it is
`-LOCPOT - reference`.

The point count is `10 × natoms` unless the sampler explicitly reports a
shortfall reallocation. Check:

```text
sampler.requested_bin_counts
sampler.actual_bin_counts
sampler.shortfall_bins
sampler.reallocated_bin_counts
sampler.bin_reallocation
```

## Example (truncated)

```json
{
  "format": "esp_tools/record-v1",
  "source": {
    "code": "vasp",
    "files": {
      "esp": {"path": "/data/water/esp_vasp.json", "sha256": "...", "bytes": 1234}
    }
  },
  "structure": {
    "natoms": 3,
    "symbols": ["O", "H", "H"],
    "counts": [1, 2],
    "formula": "OH2",
    "positions_Angstrom": [[8.0, 8.0, 8.119262], [8.0, 8.763239, 7.522953], [8.0, 7.236761, 7.522953]],
    "lattice_Angstrom": [[16, 0, 0], [0, 16, 0], [0, 0, 16]],
    "pbc": [true, true, true],
    "charge": 0.0,
    "multiplicity": null,
    "spin": {"ispin": 1, "magmom": null}
  },
  "method": {
    "code": "vasp",
    "functional": "PBE",
    "basis": "PAW",
    "software": {"name": "vasp", "version": "6.5.1", "tag": "10Mar25", "build": "Jun 15 2026 06:12:10", "variant": "gamma-only"},
    "pseudopotentials": [
      {"symbol": "O", "title": "PAW_PBE O 08Apr2002", "zval": 6.0},
      {"symbol": "H", "title": "PAW_PBE H 15Jun2001", "zval": 1.0}
    ],
    "potcar_titles": ["PAW_PBE O 08Apr2002", "PAW_PBE H 15Jun2001"],
    "settings": {"PREC": "Accurate", "ENCUT": "520", "LVHAR": ".TRUE."},
    "kpoints": "Gamma-only\n0\nGamma\n1 1 1\n0 0 0\n"
  },
  "energy": {
    "total_eV": -14.218630,
    "total_hartree": -0.5226,
    "per_atom_eV": -4.739543,
    "unit": "eV",
    "source": "OUTCAR free energy TOTEN (last ionic step)"
  },
  "forces": {
    "values_eV_per_Angstrom": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    "max_norm_eV_per_Angstrom": 0.0,
    "unit": "eV/Angstrom",
    "source": "OUTCAR TOTAL-FORCE (last ionic step)"
  },
  "stress": {
    "available": true,
    "values_kBar": [[100.0, 10.0, 30.0], [10.0, 200.0, 20.0], [30.0, 20.0, 300.0]],
    "values_eV_per_Angstrom3": [[0.0624, 0.0062, 0.0187], [0.0062, 0.1248, 0.0125], [0.0187, 0.0125, 0.1872]],
    "pressure_kBar": 200.0,
    "unit": "kBar and eV/Angstrom^3",
    "sign_convention": "VASP OUTCAR values as printed; positive trace corresponds to external pressure in the OUTCAR convention"
  },
  "electrostatic_potential": {
    "points": [[7.8, 8.1, 8.3, 1.2], [9.0, 9.1, 9.2, 0.3]],
    "point_columns": ["x_Angstrom", "y_Angstrom", "z_Angstrom", "phi_eV_per_e"],
    "point_format": "[x, y, z, q]; q = total electrostatic potential in eV/e (= V)",
    "units": {"positions": "Angstrom", "phi": "eV/e (= V)"},
    "reference": {"mode": "vacuum", "value_eV": -0.0174, "fallback": false},
    "sampler": {"npoints": 30, "requested_bin_counts": [3, 12, 12, 3], "actual_bin_counts": [3, 12, 12, 3]},
    "components_eV_per_e": {"total_eV_per_e": [1.2, 0.3], "locpot_electron_potential_energy_eV": [-1.2, -0.3]},
    "source_report": {"format": "esp_tools/vasp-esp-v1", "path": "/data/water/esp_vasp.json", "sha256": "..."}
  },
  "units": {
    "positions": "Angstrom",
    "energy": "eV",
    "forces": "eV/Angstrom",
    "stress": "kBar (raw) and eV/Angstrom^3",
    "phi": "eV/e (= V)"
  },
  "provenance": {
    "tool": "esp_tools/unify.py",
    "format": "esp_tools/record-v1",
    "created_utc": "2026-10-23T12:00:00+00:00",
    "source_code": "vasp",
    "esp_tools": {"commit": "<git sha>", "dirty": false}
  }
}
```

## `provenance`

```text
tool            "esp_tools/unify.py"
format          "esp_tools/record-v1"
created_utc     ISO-8601 UTC
source_code     "vasp" or "orca"
esp_report      original ESP JSON path
esp_tools       {commit, dirty}          # git revision of this repository
arguments       CLI arguments used to build the record
```

`method.software` records the electronic-structure code version; `source.files`
records the SHA-256 of every input/output used. VASP pseudopotential identity is
in `method.pseudopotentials` (TITEL + ZVAL) and the POTCAR hash. ORCA uses
all-electron def2 basis sets for H–Kr and ORCA's automatic def2-ECP for heavier
elements; `method.basis` records the requested basis name.

## Dataset assembly

- One `record.json` per structure.
- For JSONL: write each record as one line
  (`json.dumps(record, separators=(",", ":"))`).
- Keep the original ESP report and raw code outputs; `source.files` and
  `source_report` carry hashes for reproducibility.
- Do not split the ESP points into a separate file if the consumer expects a
  self-contained record.
- If a field is unavailable, use `null` or `{available:false, reason}`; never
  invent a value.
