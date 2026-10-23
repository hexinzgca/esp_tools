#!/usr/bin/env python3
"""Merge an ESP extraction report with VASP/ORCA outputs into one JSON record.

The extraction scripts answer "what is the potential on the 10 x natoms grid".
This tool answers the dataset question: what structure, method settings,
energy, forces, stress, and potential points belong to that grid.  It keeps
units explicit and never invents a quantity that the code did not produce
(for example, ORCA stress is reported as unavailable).

Examples
--------
VASP::

    python unify.py vasp --esp esp_vasp.json --outcar OUTCAR \
        --incar INCAR --kpoints KPOINTS --potcar POTCAR \
        --output record.json

ORCA::

    python unify.py orca --esp water_esp.json --out water.out \
        --engrad water.engrad --inp water.inp --xyz water.xyz \
        --output record.json

The output format is ``esp_tools/record-v1``.  Potential points are stored as
``[x, y, z, q]`` with ``x/y/z`` in Angstrom and ``q = phi_eV_per_e`` (the
total electrostatic potential in eV/e, numerically V).
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common.sampling import BOHR_ANGSTROM, HARTREE_EV, read_xyz  # noqa: E402
from vasp.locpot import read_poscar, sha256_file  # noqa: E402

# 1 Eh / bohr in eV / Angstrom
ORCA_GRAD_TO_EV_PER_ANGSTROM = HARTREE_EV / BOHR_ANGSTROM
# 1 kBar = 0.1 GPa; 1 eV/A^3 = 160.21766208 GPa
KBAR_TO_EV_PER_ANGSTROM3 = 0.1 / 160.21766208
RECORD_FORMAT = "esp_tools/record-v1"


def _float_or_none(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _last_float(pattern: str, text: str, flags: int = 0) -> float | None:
    matches = re.findall(pattern, text, flags)
    return float(matches[-1]) if matches else None


def _sha(path: str | Path) -> str:
    return sha256_file(path)


def _file_entry(path: str | Path) -> dict:
    path = Path(path).expanduser().resolve()
    return {"path": str(path), "sha256": _sha(path), "bytes": path.stat().st_size}


# ---------------------------------------------------------------------------
# VASP
# ---------------------------------------------------------------------------


def parse_incar(path: str | Path) -> dict:
    """Parse INCAR into an upper-case key -> string value mapping."""
    settings: dict[str, object] = {}
    for raw in Path(path).read_text(errors="replace").splitlines():
        line = raw.split("!", 1)[0].split("#", 1)[0].strip()
        if not line:
            continue
        if "=" in line:
            key, value = line.split("=", 1)
            settings[key.strip().upper()] = value.strip()
        else:
            parts = line.split()
            settings[parts[0].upper()] = " ".join(parts[1:]) or True
    return settings


def parse_potcar(path: str | Path) -> dict:
    text = Path(path).read_text(errors="replace")
    titles = re.findall(r"^\s*TITEL\s*=\s*(.+?)\s*$", text, re.MULTILINE)
    zvals = [
        float(x)
        for x in re.findall(r"ZVAL\s*=\s*([-+0-9.Ee]+)", text)
    ]
    functional = None
    if titles and all(title.startswith("PAW_PBE") for title in titles):
        functional = "PBE"
    elif titles and all(title.startswith("PAW_LDA") for title in titles):
        functional = "LDA"
    return {"titles": titles, "zvals": zvals, "functional": functional}


def _parse_magmom(value: object) -> list[float] | None:
    if value is None or value is True:
        return None
    tokens = str(value).split()
    moments: list[float] = []
    for token in tokens:
        if "*" in token:
            count_text, moment_text = token.split("*", 1)
            try:
                moments.extend([float(moment_text)] * int(count_text))
            except ValueError:
                return None
        else:
            try:
                moments.append(float(token))
            except ValueError:
                return None
    return moments


def _last_force_block(lines: list[str]) -> list[list[float]] | None:
    starts = [i for i, line in enumerate(lines) if "TOTAL-FORCE (eV/Angst)" in line]
    if not starts:
        return None
    index = starts[-1] + 1
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped and any(ch not in "- " for ch in stripped):
            break
        index += 1
    forces: list[list[float]] = []
    while index < len(lines):
        stripped = lines[index].strip()
        if not stripped or all(ch in "- " for ch in stripped):
            break
        parts = stripped.split()
        if len(parts) < 6:
            break
        try:
            forces.append([float(x) for x in parts[3:6]])
        except ValueError:
            break
        index += 1
    return forces or None


def parse_outcar(path: str | Path) -> dict:
    text = Path(path).read_text(errors="replace")
    lines = text.splitlines()

    energy = _last_float(
        r"free\s+energy\s+TOTEN\s*=\s*([-+0-9.Ee]+)\s*eV", text
    )
    if energy is None:
        energy = _last_float(
            r"energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)", text
        )

    nelect = _last_float(r"number of electron\s+([-+0-9.Ee]+)", text)

    stress_matches = re.findall(
        r"^\s*in kB\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
        text,
        re.MULTILINE,
    )
    stress = [float(x) for x in stress_matches[-1]] if stress_matches else None
    pressure = _last_float(r"external pressure\s*=\s*([-+0-9.Ee]+)\s*kB", text)

    version_match = re.search(
        r"^\s*vasp\.(\S+)\s+(\S+)\s+\(build ([^)]+)\)\s+(\S+)",
        text,
        re.MULTILINE,
    )
    software = {
        "name": "vasp",
        "version": version_match.group(1) if version_match else None,
        "tag": version_match.group(2) if version_match else None,
        "build": version_match.group(3) if version_match else None,
        "variant": version_match.group(4) if version_match else None,
    }

    return {
        "energy_eV": energy,
        "forces_eV_per_Angstrom": _last_force_block(lines),
        "stress_kBar": stress,
        "external_pressure_kBar": pressure,
        "nelect": nelect,
        "software": software,
    }


def _stress_matrix(values: list[float]) -> list[list[float]]:
    # VASP order: xx yy zz xy yz zx
    xx, yy, zz, xy, yz, zx = values
    return [[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]]


# ---------------------------------------------------------------------------
# ORCA
# ---------------------------------------------------------------------------


def _section_values(lines: list[str], header: str, count: int | None = None) -> list[float]:
    start = None
    for index, line in enumerate(lines):
        if header in line:
            start = index
            break
    if start is None:
        return []
    values: list[float] = []
    for line in lines[start + 1:]:
        stripped = line.strip()
        if stripped.startswith("#"):
            if values:
                break
            continue
        if not stripped:
            if values:
                break
            continue
        try:
            values.extend(float(x) for x in stripped.split())
        except ValueError:
            break
        if count is not None and len(values) >= count:
            break
    return values


def parse_orca_engrad(path: str | Path) -> dict:
    lines = Path(path).read_text(errors="replace").splitlines()
    natoms_values = _section_values(lines, "# Number of atoms", 1)
    energy_values = _section_values(lines, "# The current total energy in Eh", 1)
    natoms = int(natoms_values[0]) if natoms_values else None
    expected = 3 * natoms if natoms else None
    gradient = _section_values(lines, "# The current gradient in Eh/bohr", expected)
    forces = None
    if natoms and len(gradient) == expected:
        forces = [
            [-gradient[3 * i + j] * ORCA_GRAD_TO_EV_PER_ANGSTROM for j in range(3)]
            for i in range(natoms)
        ]
    return {
        "natoms": natoms,
        "energy_hartree": energy_values[0] if energy_values else None,
        "gradient_hartree_per_bohr": (
            gradient if expected is None or len(gradient) == expected else None
        ),
        "forces_eV_per_Angstrom": forces,
    }


def parse_orca_out(path: str | Path) -> dict:
    text = Path(path).read_text(errors="replace")
    energy = _last_float(r"FINAL SINGLE POINT ENERGY\s+([-+0-9.]+)", text)
    charge = _last_float(
        r"Total Charge\s+Charge\s+\.+\s+([-+0-9.]+)", text
    )
    multiplicity = _last_float(
        r"Multiplicity\s+Mult\s+\.+\s+(\d+)", text
    )
    if multiplicity is None:
        multiplicity = _last_float(r"Multiplicity\s*:\s*(\d+)", text)
    version_match = re.search(r"Program Version\s+(\S+)", text)
    libxc_match = re.search(r"libXC version:\s*(\S+)", text)
    software = {
        "name": "orca",
        "version": version_match.group(1) if version_match else None,
        "libxc_version": libxc_match.group(1) if libxc_match else None,
    }
    return {
        "energy_hartree": energy,
        "charge": int(charge) if charge is not None else None,
        "multiplicity": int(multiplicity) if multiplicity is not None else None,
        "terminated_normally": "ORCA TERMINATED NORMALLY" in text,
        "scf_converged": "SCF CONVERGED" in text,
        "software": software,
    }


def _block(text: str, name: str) -> str | None:
    match = re.search(rf"%{name}\b(.*?)\bend\b", text, re.IGNORECASE | re.DOTALL)
    return match.group(1) if match else None


def parse_orca_inp(path: str | Path) -> dict:
    text = Path(path).read_text(errors="replace")
    simple_lines = [
        line.strip()[1:].strip()
        for line in text.splitlines()
        if line.strip().startswith("!")
    ]
    simple = " ".join(simple_lines)
    tokens = simple.split()

    charge = multiplicity = None
    match = re.search(r"^\*\s*(?:xyz|xyzfile)\s+(-?\d+)\s+(\d+)", text, re.MULTILINE)
    if match:
        charge = int(match.group(1))
        multiplicity = int(match.group(2))

    method_block = _block(text, "method") or ""
    scf_block = _block(text, "scf") or ""
    cpcm_block = _block(text, "cpcm") or ""
    pal_block = _block(text, "pal") or ""

    nprocs = _last_float(r"nprocs\s+(\d+)", pal_block)
    maxcore = _last_float(r"%maxcore\s+(\d+)", text)
    grid = _last_float(r"\bGrid\s+(\d+)", method_block)
    efield = _last_float(r"\bEField\s+([-+0-9.Ee]+)", scf_block)
    epsilon = _last_float(r"\bepsilon\s+([-+0-9.Ee]+)", cpcm_block)
    refrac = _last_float(r"\brefrac\s+([-+0-9.Ee]+)", cpcm_block)

    point_charge_match = re.search(
        r"%pointcharges\s+\"?([^\"\n]+)\"?",
        text,
        re.IGNORECASE,
    )
    point_charge_file = point_charge_match.group(1).strip() if point_charge_match else None

    def is_basis(token: str) -> bool:
        prefixes = (
            "def2-", "def-", "cc-", "aug-", "ma-", "SARC",
            "STO-", "3-21G", "4-22G", "6-31", "6-311",
            "LANL", "SDD", "TZVP", "TZVPP", "SVP", "QZVP",
        )
        return token.startswith(prefixes) and not any(
            marker in token for marker in ("/J", "/C", "/JK")
        )

    def is_aux(token: str) -> bool:
        return any(marker in token for marker in ("/J", "/C", "/JK"))

    method_keywords = {
        "RIJCOSX", "RI-JK", "RI-J", "RI-MP2", "RI-MP2-F12",
        "DLPNO-CCSD(T)", "DLPNO-MP2", "VERYTIGHTSCF", "TIGHTSCF",
        "NORMALSCF", "ENGRAD", "CPCM", "NOPOP", "DECONTRACT",
    }
    basis = next((token for token in tokens if is_basis(token)), None)
    autoaux = "autoaux" in text.lower()
    aux_basis = (
        "AutoAux"
        if autoaux
        else next((token for token in tokens if is_aux(token)), None)
    )
    approximations = [
        token for token in tokens if token.upper() in method_keywords
    ]
    functional = next(
        (
            token
            for token in tokens
            if not is_basis(token)
            and not is_aux(token)
            and token.upper() not in method_keywords
            and token.upper() != "AUTOAUX"
        ),
        None,
    )

    return {
        "simple_input": simple,
        "tokens": tokens,
        "functional": functional,
        "basis": basis,
        "aux_basis": aux_basis,
        "approximations": approximations,
        "charge": charge,
        "multiplicity": multiplicity,
        "nprocs": int(nprocs) if nprocs is not None else None,
        "maxcore_mb": int(maxcore) if maxcore is not None else None,
        "grid": int(grid) if grid is not None else None,
        "efield_au": efield,
        "point_charges_file": point_charge_file,
        "cpcm": (
            {"epsilon": epsilon, "refraction_index": refrac}
            if cpcm_block
            else None
        ),
    }


# ---------------------------------------------------------------------------
# Record assembly
# ---------------------------------------------------------------------------


def read_esp_report(path: str | Path) -> dict:
    report = json.loads(Path(path).read_text())
    if not isinstance(report, dict) or "esp" not in report:
        raise ValueError(f"not an esp_tools ESP report: {path}")
    esp = report["esp"]
    if "positions_Angstrom" not in esp or "phi_eV_per_e" not in esp:
        raise ValueError(
            f"{path} has no phi_eV_per_e; run extraction with a wavefunction/GBW"
        )
    positions = np.asarray(esp["positions_Angstrom"], dtype=float)
    phi = np.asarray(esp["phi_eV_per_e"], dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("ESP positions must be Nx3")
    if phi.shape != (len(positions),):
        raise ValueError("ESP positions and phi lengths differ")
    if not np.isfinite(positions).all() or not np.isfinite(phi).all():
        raise ValueError("ESP report contains non-finite values")
    return report


def _structure_from_esp(report: dict) -> dict:
    structure = report.get("structure")
    if not isinstance(structure, dict) or "positions_Angstrom" not in structure:
        raise ValueError("ESP report has no structure; pass --poscar or --xyz")
    return structure


def _check_positions(expected: np.ndarray, actual: np.ndarray, label: str) -> None:
    if expected.shape != actual.shape or not np.allclose(
        expected, actual, atol=1e-3, rtol=0
    ):
        raise ValueError(f"{label} positions do not match the ESP report")


def _max_norm(values: list[list[float]] | None) -> float | None:
    if not values:
        return None
    return float(np.linalg.norm(np.asarray(values, dtype=float), axis=1).max())


def _counts_from_symbols(symbols: list[str]) -> list[int]:
    counts: dict[str, int] = {}
    order: list[str] = []
    for symbol in symbols:
        if symbol not in counts:
            counts[symbol] = 0
            order.append(symbol)
        counts[symbol] += 1
    return [counts[symbol] for symbol in order]


def _esp_block(report: dict, source_path: Path) -> dict:
    esp = report["esp"]
    positions = np.asarray(esp["positions_Angstrom"], dtype=float)
    phi = np.asarray(esp["phi_eV_per_e"], dtype=float)
    points = [
        [float(x), float(y), float(z), float(q)]
        for (x, y, z), q in zip(positions, phi)
    ]
    code = report.get("format", "")
    if "vasp" in code:
        reference = report.get("locpot", {}).get("reference")
        components = {
            "total_eV_per_e": phi.tolist(),
            "locpot_electron_potential_energy_eV": esp.get("phi_locpot_eV"),
        }
    else:
        reference = {
            "mode": "infinity",
            "value_eV": 0.0,
            "description": "orca_vpot references the potential to infinity",
        }
        components = {
            "total_eV_per_e": phi.tolist(),
            "molecular_eV_per_e": esp.get("phi_molecular_eV_per_e"),
            "background_eV_per_e": esp.get("phi_background_eV_per_e"),
            "uniform_field_eV_per_e": esp.get("phi_uniform_field_eV_per_e"),
        }
    return {
        "points": points,
        "point_columns": [
            "x_Angstrom",
            "y_Angstrom",
            "z_Angstrom",
            "phi_eV_per_e",
        ],
        "point_format": "[x, y, z, q]; q = total electrostatic potential in eV/e (= V)",
        "units": {"positions": "Angstrom", "phi": "eV/e (= V)"},
        "reference": reference,
        "sampler": report.get("sampler"),
        "components_eV_per_e": components,
        "source_report": {
            "format": report.get("format"),
            "path": str(source_path.resolve()),
            "sha256": _sha(source_path),
        },
    }


def build_vasp_record(args: argparse.Namespace) -> dict:
    esp_path = Path(args.esp).expanduser()
    report = read_esp_report(esp_path)
    esp_structure = _structure_from_esp(report)
    natoms = int(esp_structure["natoms"])
    esp_positions = np.asarray(esp_structure["positions_Angstrom"], dtype=float)

    if args.poscar:
        geometry = read_poscar(args.poscar)
        _check_positions(esp_positions, geometry.positions, "POSCAR")
        lattice = geometry.lattice.tolist()
        symbols = geometry.symbols
        positions = geometry.positions.tolist()
        counts = geometry.counts
    else:
        lattice = esp_structure.get("lattice_Angstrom")
        symbols = esp_structure.get("symbols", [])
        positions = esp_positions.tolist()
        counts = esp_structure.get("counts") or _counts_from_symbols(symbols)

    outcar = parse_outcar(args.outcar)
    incar = parse_incar(args.incar) if args.incar else {}
    potcar = parse_potcar(args.potcar) if args.potcar else {"titles": [], "zvals": [], "functional": None}
    kpoints_text = Path(args.kpoints).read_text(errors="replace") if args.kpoints else None

    if outcar["forces_eV_per_Angstrom"] is not None and len(outcar["forces_eV_per_Angstrom"]) != natoms:
        raise ValueError(
            f"OUTCAR has {len(outcar['forces_eV_per_Angstrom'])} force rows, "
            f"ESP structure has {natoms} atoms"
        )

    nelect = outcar["nelect"]
    if nelect is None and incar.get("NELECT") not in (None, True):
        try:
            nelect = float(incar["NELECT"])
        except (TypeError, ValueError):
            nelect = None
    charge = None
    if potcar["zvals"] and nelect is not None:
        if len(counts) == len(potcar["zvals"]):
            valence = sum(
                int(count) * zval for count, zval in zip(counts, potcar["zvals"])
            )
            charge = float(valence - nelect)
            if abs(charge - round(charge)) < 1e-4:
                charge = float(round(charge))
    ispin = incar.get("ISPIN")
    try:
        ispin = int(str(ispin)) if ispin is not None else None
    except ValueError:
        ispin = None
    magmom = _parse_magmom(incar.get("MAGMOM"))

    stress = None
    if outcar["stress_kBar"] is not None:
        kbar = _stress_matrix(outcar["stress_kBar"])
        ev_a3 = [
            [value * KBAR_TO_EV_PER_ANGSTROM3 for value in row] for row in kbar
        ]
        pressure = outcar["external_pressure_kBar"]
        if pressure is None:
            pressure = sum(outcar["stress_kBar"][:3]) / 3.0
        stress = {
            "available": True,
            "values_kBar": kbar,
            "values_eV_per_Angstrom3": ev_a3,
            "pressure_kBar": pressure,
            "unit": "kBar and eV/Angstrom^3",
            "sign_convention": (
                "VASP OUTCAR values as printed; a positive trace corresponds "
                "to external pressure in the OUTCAR convention"
            ),
        }
    else:
        stress = {"available": False, "reason": "OUTCAR has no 'in kB' stress block"}

    energy = None
    if outcar["energy_eV"] is not None:
        energy = {
            "total_eV": outcar["energy_eV"],
            "total_hartree": outcar["energy_eV"] / HARTREE_EV,
            "per_atom_eV": outcar["energy_eV"] / natoms,
            "unit": "eV",
            "source": "OUTCAR free energy TOTEN (last ionic step)",
        }

    files = {"esp": _file_entry(esp_path), "outcar": _file_entry(args.outcar)}
    for name in ("incar", "kpoints", "potcar"):
        value = getattr(args, name)
        if value:
            files[name] = _file_entry(value)

    pseudopotentials = []
    for index, title in enumerate(potcar["titles"]):
        pseudopotentials.append(
            {
                "symbol": symbols[index] if index < len(symbols) else None,
                "title": title,
                "zval": potcar["zvals"][index] if index < len(potcar["zvals"]) else None,
            }
        )

    return {
        "format": RECORD_FORMAT,
        "source": {"code": "vasp", "files": files},
        "structure": {
            "natoms": natoms,
            "symbols": symbols,
            "counts": counts,
            "formula": esp_structure.get("formula"),
            "positions_Angstrom": positions,
            "lattice_Angstrom": lattice,
            "pbc": [True, True, True] if lattice is not None else [False, False, False],
            "charge": charge,
            "multiplicity": None,
            "spin": {"ispin": ispin, "magmom": magmom},
        },
        "method": {
            "code": "vasp",
            "functional": potcar["functional"] or incar.get("GGA"),
            "basis": "PAW",
            "software": outcar["software"],
            "pseudopotentials": pseudopotentials,
            "potcar_titles": potcar["titles"],
            "settings": incar,
            "kpoints": kpoints_text,
        },
        "energy": energy,
        "forces": (
            {
                "values_eV_per_Angstrom": outcar["forces_eV_per_Angstrom"],
                "max_norm_eV_per_Angstrom": _max_norm(outcar["forces_eV_per_Angstrom"]),
                "unit": "eV/Angstrom",
                "source": "OUTCAR TOTAL-FORCE (last ionic step)",
            }
            if outcar["forces_eV_per_Angstrom"] is not None
            else None
        ),
        "stress": stress,
        "electrostatic_potential": _esp_block(report, esp_path),
        "units": {
            "positions": "Angstrom",
            "energy": "eV",
            "forces": "eV/Angstrom",
            "stress": "kBar (raw) and eV/Angstrom^3",
            "phi": "eV/e (= V)",
        },
        "provenance": _provenance("vasp", esp_path, args),
    }


def build_orca_record(args: argparse.Namespace) -> dict:
    esp_path = Path(args.esp).expanduser()
    report = read_esp_report(esp_path)
    esp_structure = _structure_from_esp(report)
    natoms = int(esp_structure["natoms"])
    esp_positions = np.asarray(esp_structure["positions_Angstrom"], dtype=float)

    if args.xyz:
        symbols, xyz_positions = read_xyz(args.xyz)
        _check_positions(esp_positions, xyz_positions, "XYZ")
        positions = xyz_positions.tolist()
    else:
        symbols = esp_structure.get("symbols", [])
        positions = esp_positions.tolist()
    counts = esp_structure.get("counts") or _counts_from_symbols(symbols)

    out = parse_orca_out(args.out)
    engrad = parse_orca_engrad(args.engrad) if args.engrad else {
        "natoms": None,
        "energy_hartree": None,
        "gradient_hartree_per_bohr": None,
        "forces_eV_per_Angstrom": None,
    }
    inp = parse_orca_inp(args.inp) if args.inp else {
        "simple_input": None,
        "tokens": [],
        "functional": None,
        "basis": None,
        "aux_basis": None,
        "approximations": [],
        "charge": None,
        "multiplicity": None,
        "nprocs": None,
        "maxcore_mb": None,
        "grid": None,
        "efield_au": None,
        "point_charges_file": None,
        "cpcm": None,
    }

    if engrad["natoms"] is not None and int(engrad["natoms"]) != natoms:
        raise ValueError(
            f"engrad has {engrad['natoms']} atoms, ESP structure has {natoms}"
        )
    if engrad["forces_eV_per_Angstrom"] is not None and len(
        engrad["forces_eV_per_Angstrom"]
    ) != natoms:
        raise ValueError("engrad force count does not match the structure")

    energy_hartree = out["energy_hartree"]
    if energy_hartree is None:
        energy_hartree = engrad["energy_hartree"]
    energy = None
    if energy_hartree is not None:
        energy = {
            "total_hartree": energy_hartree,
            "total_eV": energy_hartree * HARTREE_EV,
            "per_atom_eV": energy_hartree * HARTREE_EV / natoms,
            "unit": "eV",
            "source": "ORCA FINAL SINGLE POINT ENERGY (last value)",
        }

    charge = inp["charge"] if inp["charge"] is not None else out["charge"]
    multiplicity = (
        inp["multiplicity"]
        if inp["multiplicity"] is not None
        else out["multiplicity"]
    )

    files = {"esp": _file_entry(esp_path), "out": _file_entry(args.out)}
    for name in ("engrad", "inp", "xyz"):
        value = getattr(args, name)
        if value:
            files[name] = _file_entry(value)

    software = dict(out["software"])
    orca_provenance = report.get("provenance", {})
    if isinstance(orca_provenance, dict):
        orca_runtime = orca_provenance.get("orca", {})
        if isinstance(orca_runtime, dict) and orca_runtime.get("container"):
            software["container"] = orca_runtime["container"]

    return {
        "format": RECORD_FORMAT,
        "source": {"code": "orca", "files": files},
        "structure": {
            "natoms": natoms,
            "symbols": symbols,
            "counts": counts,
            "formula": esp_structure.get("formula"),
            "positions_Angstrom": positions,
            "lattice_Angstrom": None,
            "pbc": [False, False, False],
            "charge": charge,
            "multiplicity": multiplicity,
            "spin": {"multiplicity": multiplicity},
        },
        "method": {
            "code": "orca",
            "functional": inp["functional"],
            "basis": inp["basis"],
            "aux_basis": inp["aux_basis"],
            "software": software,
            "approximations": inp["approximations"],
            "solvent": "CPCM" if inp["cpcm"] else None,
            "settings": {
                "simple_input": inp["simple_input"],
                "tokens": inp["tokens"],
                "nprocs": inp["nprocs"],
                "maxcore_mb": inp["maxcore_mb"],
                "grid": inp["grid"],
                "efield_au": inp["efield_au"],
                "point_charges_file": inp["point_charges_file"],
                "cpcm": inp["cpcm"],
                "scf_converged": out["scf_converged"],
                "terminated_normally": out["terminated_normally"],
            },
        },
        "energy": energy,
        "forces": (
            {
                "values_eV_per_Angstrom": engrad["forces_eV_per_Angstrom"],
                "max_norm_eV_per_Angstrom": _max_norm(
                    engrad["forces_eV_per_Angstrom"]
                ),
                "unit": "eV/Angstrom",
                "source": (
                    "ORCA engrad: force = -gradient(Eh/bohr) * "
                    f"{ORCA_GRAD_TO_EV_PER_ANGSTROM:.9f}"
                ),
            }
            if engrad["forces_eV_per_Angstrom"] is not None
            else None
        ),
        "stress": {
            "available": False,
            "reason": (
                "ORCA non-periodic molecular calculations do not define a "
                "stress tensor"
            ),
        },
        "electrostatic_potential": _esp_block(report, esp_path),
        "units": {
            "positions": "Angstrom",
            "energy": "eV",
            "forces": "eV/Angstrom",
            "stress": "not available for non-periodic ORCA calculations",
            "phi": "eV/e (= V)",
        },
        "provenance": _provenance("orca", esp_path, args),
    }


def _esp_tools_info() -> dict:
    try:
        commit = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        status = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        return {
            "commit": commit.stdout.strip() if commit.returncode == 0 else None,
            "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
        }
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None}


def _provenance(code: str, esp_path: Path, args: argparse.Namespace) -> dict:
    return {
        "tool": "esp_tools/unify.py",
        "format": RECORD_FORMAT,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_code": code,
        "esp_report": str(esp_path.resolve()),
        "esp_tools": _esp_tools_info(),
        "arguments": {
            key: str(value)
            for key, value in vars(args).items()
            if key not in {"output", "force", "func"} and value is not None
        },
    }


def _write_json(path: str | Path, record: dict, force: bool) -> Path:
    output = Path(path).expanduser()
    if output.exists() and not force:
        raise FileExistsError(
            f"output already exists: {output}; choose another name or pass --force"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w" if force else "x") as handle:
        json.dump(record, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="code", required=True)

    vasp = subparsers.add_parser("vasp", help="Merge a VASP ESP report with OUTCAR data")
    vasp.add_argument("--esp", required=True, help="ESP JSON from vasp/extract_esp.py")
    vasp.add_argument("--outcar", required=True, help="OUTCAR with final energy/forces/stress")
    vasp.add_argument("--incar", help="INCAR (method settings)")
    vasp.add_argument("--kpoints", help="KPOINTS (raw text is stored)")
    vasp.add_argument("--potcar", help="POTCAR (titles, ZVAL, functional)")
    vasp.add_argument("--poscar", help="POSCAR (structure cross-check)")
    vasp.add_argument("--output", required=True, help="Unified record JSON")
    vasp.add_argument("--force", action="store_true")
    vasp.set_defaults(func=build_vasp_record)

    orca = subparsers.add_parser("orca", help="Merge an ORCA ESP report with out/engrad data")
    orca.add_argument("--esp", required=True, help="ESP JSON from orca/extract_esp.py")
    orca.add_argument("--out", required=True, help="ORCA .out file")
    orca.add_argument("--engrad", help="ORCA .engrad file (forces)")
    orca.add_argument("--inp", help="ORCA .inp file (functional/basis/charge/spin)")
    orca.add_argument("--xyz", help="XYZ geometry (structure cross-check)")
    orca.add_argument("--output", required=True, help="Unified record JSON")
    orca.add_argument("--force", action="store_true")
    orca.set_defaults(func=build_orca_record)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        record = args.func(args)
        output = _write_json(args.output, record, args.force)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(
        f"Wrote {record['source']['code']} record with "
        f"{record['structure']['natoms']} atoms and "
        f"{len(record['electrostatic_potential']['points'])} potential points "
        f"to {output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
