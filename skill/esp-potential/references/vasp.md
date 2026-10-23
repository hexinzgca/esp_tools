# VASP reference: LOCPOT and ESP extraction

## 1. What `LVHAR` writes

| INCAR | LOCPOT content | ESP? |
|---|---|---|
| `LVHAR = .TRUE.` | $V_{\rm ionic} + V_{\rm Hartree}$ in eV | **yes** |
| `LVTOT = .TRUE.` | adds $V_{\rm XC}$ | no |

`LVHAR=.TRUE.` forces `LVTOT=.FALSE.`. The values are the **potential energy of
an electron**, negative near nuclei. The extraction script converts:

$$
\phi_{\rm ESP}[{\rm eV}/e] = -V_{\rm LOCPOT}[{\rm eV}] - V_{\rm ref}[{\rm eV}]
$$

The raw values are kept as `esp.phi_locpot_eV`; the physical ESP is
`esp.phi_eV_per_e`. Compare with ORCA `orca_vpot` only after this sign flip.

`ISPIN=2` writes two blocks; under `LVHAR` their ionic/Hartree parts should be
identical. The parser averages blocks and records
`locpot.spin_block_max_abs_diff_eV`. A large difference suggests `LVTOT` or a
different LOCPOT source.

## 2. Input files

### POSCAR

- Isolated molecule: large cubic box, e.g. 16 Å for water. The `>=5 Å` far bin
  needs volume; 16 Å is comfortable. Center the molecule.
- Element order in POSCAR determines POTCAR order.
- Cartesian or Direct coordinates are both supported by the parser.

### POTCAR

Concatenate in **POSCAR species order**. Example for `O H`:

```text
/opt/vasp_pot/O/POTCAR + /opt/vasp_pot/H/POTCAR
```

The parser reads `TITEL` and `ZVAL` from POTCAR. If POTCAR is supplied to
`unify.py`, it computes:

```text
charge = sum(count_i * ZVAL_i) - NELECT
```

`NELECT` is read from the last `number of electron` line in OUTCAR. If POTCAR or
NELECT is missing, charge is `null` rather than guessed.

### INCAR

See `vasp/INCAR.example`. Core tags:

| Tag | Value | Why |
|---|---|---|
| `PREC` | `Accurate` | fine FFT grid for LOCPOT |
| `ENCUT` | 520 (or ≥1.3×ENMAX) | converged basis |
| `EDIFF` | `1E-6` | SCF convergence |
| `NELM` | 200 | enough SCF cycles |
| `ISMEAR` / `SIGMA` | 0 / 0.05 | molecules |
| `ISPIN` | 1 or 2 | spin treatment |
| `ISYM` | 0 | avoid symmetry reductions for molecules |
| `LREAL` | `.FALSE.` | accurate real-space projection |
| `IBRION` / `NSW` | -1 / 0 | static |
| `LWAVE` / `LCHARG` | `.FALSE.` | not needed for ESP |
| `LVHAR` | `.TRUE.` | write electrostatic potential |
| `LVTOT` | `.FALSE.` | do not add XC |
| `NCORE` / `KPAR` | 4 / 1 | adjust to machine |

For a molecule, use Γ-only:

```text
Gamma-only
0
Gamma
1 1 1
0 0 0
```

and `vasp_gam`. For a periodic system with k-points, use `vasp_std`.

### Box, charge, spin, dipole

- Use a **neutral** cell for `--reference vacuum`. A charged cell has no finite
  vacuum level.
- A large box is required so the far-field potential is flat. If the cell has a
  net dipole, the potential may have a slope; a dipole correction changes the
  potential and must be applied consistently if used.
- `ISPIN=2` is collinear spin. VASP does not provide a molecular multiplicity;
  the unified record stores `spin.ispin` and `spin.magmom`, with
  `multiplicity=null`.
- `NELECT`/background-charge calculations (`LCHRED`) do not have a vacuum
  reference; use `--reference none` and state the zero.

## 3. Run

```bash
apptainer exec --bind "$PWD:/work" \
  /mnt/e/apptainer_sifs/dft1/vasp-cpu-intel.sif \
  bash -c 'cd /work && mpirun -np 4 /opt/vasp/bin/vasp_gam > vasp.out 2>&1'
```

The Intel image has the full POTCAR library under `/opt/vasp_pot/<Element>/POTCAR`.
On this WSL host, `vasp_gam` in the Intel image has shown SIGSEGV; the GCC image
works if POTCARs are supplied separately. The extraction is independent of the
image.

Check:

```bash
grep -E "free  energy|reached required accuracy" vasp.out
grep -n "NGXF" OUTCAR | head
head -10 LOCPOT
```

## 4. Extract

```bash
cd esp_tools
python3 vasp/extract_esp.py --locpot LOCPOT --poscar POSCAR --output esp_vasp.json
```

Options:

| Option | Default | Notes |
|---|---|---|
| `--sample-mode` | `dft-grid` | actual LOCPOT grid points; no interpolation |
| `--points-per-atom` | 10 | total = value × natoms |
| `--reference` | `vacuum` | mean of `>=5 Å` points; `none`; or explicit eV |
| `--seed` | 20260916 | geometry-hashed |
| `--min-distance` | 0.25 Å | minimum distance to any atom |
| `--far-span` | 5 Å | `interpolate` mode only |
| `--force` | off | overwrite |

`dft-grid` scans the LOCPOT grid, classifies points by minimum-image distance,
and samples the 1:4:4:1 quotas. Each point is an exact `(i,j,k)` grid point and
`phi_locpot_eV` is read directly. `interpolate` generates arbitrary points and
uses periodic trilinear interpolation.

Accuracy test (8 Å water, 96³ vs 192³, 30 fine-grid references):

| Bin | Points | Mean abs error / eV | Max abs error / eV |
|---|---:|---:|---:|
| `[0.25,0.5)` | 3 | 0.563 | 0.846 |
| `[0.5,1.5)` | 12 | 0.0276 | 0.0948 |
| `[1.5,5)` | 12 | 0.000160 | 0.00131 |
| `[5,∞)` | 3 | 0.000008 | 0.000022 |
| all | 30 | 0.0675 | 0.846 |

Interpolation at its own grid nodes is exact (max difference 0.0). Use
`dft-grid` when exact DFT grid values are required; use `interpolate` only for
arbitrary-point comparisons.

Reproduce the comparison:

```bash
python vasp/accuracy_check.py \
  --coarse-locpot coarse/LOCPOT \
  --fine-locpot fine/LOCPOT \
  --fine-esp fine/grid.json \
  --output accuracy.json
```

The fine ESP report must come from `--sample-mode dft-grid`.

The sampler:

- uses `10 × natoms`, bins `[0.25,0.5) [0.5,1.5) [1.5,5) [5,∞)` with target
  `1:4:4:1`;
- wraps all points into the cell and uses minimum-image distances;
- samples far-field points uniformly in the cell and accepts only `>=5 Å`;
- if `[5,∞)` cannot be filled, moves the deficit to `[1.5,5)` first, then lower
  bins only if necessary; total remains `10×N`; prints WARNING;
- records `requested_bin_counts`, `actual_bin_counts`, `shortfall_bins`,
  `reallocated_bin_counts`, `bin_reallocation`.

If `--reference vacuum` finds **zero** far-field points, it falls back to
`--reference none` (LOCPOT zero), prints WARNING, and records
`reference.requested_mode="vacuum"`, `reference.mode="none"`,
`reference.fallback=true`. It does not abort.

## 5. Output

```text
structure            lattice, positions, symbols, counts
sampler              bin counts, seed, edges, reallocation
locpot               grid_shape, n_blocks, spin diff, reference, sign_convention
esp.positions_Angstrom
esp.grid_indices_ijk dft-grid mode: LOCPOT (i,j,k) indices
esp.phi_locpot_eV    raw electron potential energy (direct or interpolated)
esp.phi_eV_per_e     ESP = -LOCPOT - reference
```

`locpot.reference.value_eV` is the **electrostatic** reference. For
`--reference vacuum`, it equals `-mean(LOCPOT far)`. `locpot_mean_eV` stores the
raw LOCPOT mean for transparency.

## 6. PAW near-core caveat

VASP `LOCPOT` is built from PAW/pseudopotential potentials. The first bin
`[0.25,0.5)` can lie inside augmentation spheres, so it is not an all-electron
ESP. Consequences:

- use it for ML consistency within the same PAW setup;
- for point-by-point comparison with ORCA `orca_vpot`, exclude the first bin or
  only use `>=0.5 Å` points;
- do not claim the VASP near-core values are all-electron.

## 7. Validation in this repository

- VASP 6.5.1, water in a 16 Å box, Γ-only, PBE PAW: `LOCPOT` grid
  `240×240×240`; 30 probes; bins `[3,12,12,3]`; far reference `0.0174 eV`.
- VASP 6.5.1, 4 Å single-O box: no `>=5 Å` region; `requested [1,4,4,1] ->
  actual [1,4,5,0]`; `reallocated_bin_counts [0,0,1,0]`; `--reference vacuum`
  fell back to `none`; process exited 0.

## 8. Common failures

| Symptom | Cause/action |
|---|---|
| no `LOCPOT` | `LVHAR` not on, or SCF failed |
| `spin_block_max_abs_diff_eV` large | likely `LVTOT` or non-LVHAR LOCPOT |
| `reference` changes with box | insufficient vacuum or charged cell |
| `bin_reallocation=true` | dense/small cell; inspect counts and fallback |
| near-core phi huge/negative | PAW pseudo-potential; see §6 |
| POSCAR/LOCPOT mismatch | use the same geometry; script refuses mixed inputs |

## 9. Version metadata in the unified record

`unify.py` reads the VASP version from the OUTCAR header:

```text
vasp.6.5.1 10Mar25 (build Jun 15 2026 06:12:10) gamma-only
```

and writes `method.software = {name, version, tag, build, variant}`.
POTCAR `TITEL` and `ZVAL` are paired with the POSCAR symbols in
`method.pseudopotentials`; the POTCAR file SHA-256 is in
`source.files.potcar`. `provenance.esp_tools.commit/dirty` records the
extraction tool revision.
