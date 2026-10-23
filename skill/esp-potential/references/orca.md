# ORCA reference: wB97M-V + diffuse def2 + orca_vpot

## 1. Basis

Default to `def2-TZVPD` (triple-zeta + diffuse, OMol25 level); use `def2-SVPD`
for a double-zeta + diffuse alternative.

| Keyword | Level | Recommendation |
|---|---|---|
| `def2-SVPD` | double-zeta + diffuse | alternative |
| `def2-TZVPD` | triple-zeta + diffuse | **default** |
| `def2-TZVPPD` | triple-zeta, double polarization + diffuse | higher accuracy |
| `def2-QZVPD` | quadruple-zeta + diffuse | reference |

## 2. SCF input

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
  exchange.
- `VeryTightSCF`: recommended; ESP is sensitive to the density.
- `Grid 5`: matches the existing project; `Grid 6` is more accurate/slower.
- `%pal`: 4 is a safe default; increase to the machine's physical cores.
- ORCA requires a **full pathname** when `%pal nprocs > 1`; the runner uses
  `ORCA=$(command -v orca); exec "$ORCA" ...`.
- ORCA 6.1 writes `water.gbw` and `water.densities` after SCF. `orca_vpot` needs
  both files in the same directory.

For an anion/cation, edit the charge/multiplicity in `* xyzfile`. Reuse the same
probe points (`--points`) only for the same geometry.

## 3. Run

```bash
cd esp_tools/orca
./run_orca_esp.sh wb97m-v_def2-tzvpd.inp water.xyz
```

The runner:

1. runs ORCA in `/mnt/e/apptainer_sifs/build/orca-6.1.sif`;
2. writes `wb97m-v_def2-tzvpd.out`, `.gbw`, `.densities`;
3. calls `extract_esp.py` with `--orca-container`.

Step-by-step:

```bash
apptainer exec --bind "$PWD:/work" \
  /mnt/e/apptainer_sifs/build/orca-6.1.sif \
  bash -c 'cd /work; export OMPI_MCA_rmaps_base_oversubscribe=1; ORCA=$(command -v orca); exec "$ORCA" wb97m-v_def2-tzvpd.inp' \
  > wb97m-v_def2-tzvpd.out 2>&1

cd ..
python3 orca/extract_esp.py --xyz orca/water.xyz \
  --gbw orca/wb97m-v_def2-tzvpd.gbw \
  --orca-container /mnt/e/apptainer_sifs/build/orca-6.1.sif \
  --output orca/water_esp.json
```

Native ORCA:

```bash
python3 orca/extract_esp.py --xyz orca/water.xyz --gbw water.gbw \
  --orca-home /opt/orca/orca-6.1.0-f.0_linux_x86-64/bin \
  --output water_esp.json
```

## 4. How the density is found

ORCA 6.1 stores SCF densities in a `.densities` container. The script runs:

```bash
orca_vpot water.densities
```

to list the contents, selects the unique `.scfp`, then evaluates:

```bash
orca_vpot water.gbw water.scfp grid_bohr.txt esp_bohr.txt
```

This avoids guessing a density name and avoids the common segfault from passing
an invalid `.scfp`. For older ORCA versions with a standalone `water.scfp`, the
script automatically falls back to that file.

The grid file has a count header followed by `x y z` in Bohr; the output has a
count header followed by `x y z V` with `V` in Hartree/e. The script converts to
Angstrom and eV/e:

```text
phi_eV_per_e = V_Hartree_per_e * 27.211386245988
```

The molecular potential is referenced to infinity. ORCA's `phi_molecular` already
contains nuclei + SCF electrons; do not add another `Q/r` for the nuclei.

## 5. Probe grid

- Exactly `10 × natoms`; bins `[0.25,0.5) [0.5,1.5) [1.5,5) [5,∞)` Å;
  target `1:4:4:1`.
- Far points lie on a centroid shell `r_max+5 … r_max+10 Å`.
- Non-periodic Euclidean distances; the first bin starts at 0.25 Å, so no point
  is on a nucleus.
- Same geometry reproduces the same grid (geometry-hashed seed).

Sample only / reuse:

```bash
python3 orca/extract_esp.py --xyz orca/water.xyz --output water_grid.json
python3 orca/extract_esp.py --points water_grid.json --gbw water.gbw \
  --orca-container /path/orca-6.1.sif --output water_esp.json
```

`--points` does not re-sample and does not validate the geometry; the grid must
come from the same structure as the GBW.

## 6. Background charges and electric field

- `--background background.json`: JSON `[[x_A,y_A,z_A,q_e], ...]`. This is
  **not** ORCA's `q x y z` point-charge file order.
- `--field EX EY EZ`: uniform field in atomic units, gauge origin `(0,0,0)`.
- The extractor adds the external **source** potential. The electronic response
  must come from a SCF that included `%pointcharges` or `%scf EField`; otherwise
  the result is only a superposition, not the self-consistent response.

Output keys:

```text
esp.positions_Angstrom
esp.phi_molecular_eV_per_e
esp.phi_background_eV_per_e
esp.phi_uniform_field_eV_per_e
esp.phi_eV_per_e              # sum
provenance                    # GBW, density name, commands, raw files
sampler                       # bins, seed, reallocation (not used by default)
```

## 7. Validation in this repository

- ORCA 6.1.0 in Apptainer, water, `wB97M-V/def2-TZVPD/AutoAux RIJCOSX
  VeryTightSCF`, Grid 5: SCF energy `-76.433367884141 Ha`, normal termination.
- `.densities` listed `water.scfp`; `orca_vpot` evaluated 30 probes; bins
  `[3,12,12,3]`; finite `phi` from `-0.087077` to `24.677244 eV/e`.
- `def2-SVPD` input also runs normally.

## 8. Common failures

| Symptom | Cause/action |
|---|---|
| `ORCA has to be called with full pathname` | use the runner/full path |
| no `.densities` | SCF did not finish, or ORCA version too old; fallback `.scfp` |
| `orca_vpot` segfault / density dimension error | do not invent `.scfp`; let the script list `.densities` |
| `no point in >=5 Å far bin` (ORCA) | far shell should always exist; check sampler params |
| values differ from VASP near nuclei | PAW vs all-electron; exclude `0.25–0.5 Å` |
| charged system | potential reference is infinity only for the isolated system; state the convention |

## 9. Version metadata in the unified record

`unify.py` reads `Program Version` and the libXC version from the ORCA `.out`:

```text
Program Version 6.1.0-f.0  -   RELEASE   -
Your ORCA version has been built with support for libXC version: 7.0.0
```

and writes `method.software = {name, version, libxc_version, container}`.
The container path comes from the ESP report provenance when
`--orca-container` was used. `provenance.esp_tools.commit/dirty` records the
extraction tool revision.
