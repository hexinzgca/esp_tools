from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import unify  # noqa: E402

try:
    import jsonschema
except ImportError:  # validation is optional; numpy is the only runtime dependency
    jsonschema = None

SCHEMA = json.loads((ROOT / "schema" / "record-v1.schema.json").read_text())


def validate_record(record: dict) -> None:
    if jsonschema is not None:
        jsonschema.validate(record, SCHEMA)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload))


def esp_payload(structure: dict, positions: list, phi: list, vasp: bool) -> dict:
    report = {
        "format": "esp_tools/vasp-esp-v1" if vasp else "esp_tools/orca-esp-v1",
        "structure": structure,
        "sampler": {"npoints": len(positions)},
        "esp": {"positions_Angstrom": positions, "phi_eV_per_e": phi},
    }
    if vasp:
        report["esp"]["phi_locpot_eV"] = [-value for value in phi]
        report["locpot"] = {"reference": {"mode": "vacuum", "value_eV": 0.0}}
    else:
        report["esp"]["phi_molecular_eV_per_e"] = phi
        report["esp"]["phi_background_eV_per_e"] = [0.0] * len(phi)
        report["esp"]["phi_uniform_field_eV_per_e"] = [0.0] * len(phi)
    return report


class UnifyTests(unittest.TestCase):
    def test_vasp_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            structure = {
                "natoms": 2,
                "symbols": ["O", "H"],
                "counts": [1, 1],
                "formula": "OH",
                "lattice_Angstrom": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "positions_Angstrom": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            }
            esp_path = tmp_path / "esp.json"
            write_json(
                esp_path,
                esp_payload(
                    structure,
                    [[2.0, 2.0, 2.0], [3.0, 3.0, 3.0]],
                    [1.0, 2.0],
                    vasp=True,
                ),
            )
            outcar = tmp_path / "OUTCAR"
            outcar.write_text(
                " vasp.6.5.1 10Mar25 (build Jun 15 2026 06:12:10) gamma-only\n"
                " free  energy   TOTEN  =      -10.00000000 eV\n"
                " number of electron       7.0000000 magnetization\n"
                " POSITION                                       TOTAL-FORCE (eV/Angst)\n"
                " -----------------------------------------------------------------------------------\n"
                "      0.00000      0.00000      0.00000        0.100000      0.200000      0.300000\n"
                "      1.00000      0.00000      0.00000       -0.100000     -0.200000     -0.300000\n"
                " -----------------------------------------------------------------------------------\n"
                "  in kB     100.00000   200.00000   300.00000    10.00000    20.00000    30.00000\n"
                "  external pressure =      200.00 kB  Pullay stress =        0.00 kB\n"
            )
            incar = tmp_path / "INCAR"
            incar.write_text("ISPIN = 2\nMAGMOM = 1.0 -1.0\nENCUT = 520\n")
            potcar = tmp_path / "POTCAR"
            potcar.write_text(
                "TITEL = PAW_PBE O 08Apr2002\nPOMASS = 16.0; ZVAL = 6.000\n"
                "TITEL = PAW_PBE H 15Jun2001\nPOMASS = 1.0; ZVAL = 1.000\n"
            )
            output = tmp_path / "record.json"
            rc = unify.main([
                "vasp",
                "--esp", str(esp_path),
                "--outcar", str(outcar),
                "--incar", str(incar),
                "--potcar", str(potcar),
                "--output", str(output),
            ])
            self.assertEqual(rc, 0)
            record = json.loads(output.read_text())
            validate_record(record)
            self.assertEqual(record["format"], unify.RECORD_FORMAT)
            self.assertEqual(record["structure"]["charge"], 0.0)
            self.assertEqual(record["structure"]["spin"], {"ispin": 2, "magmom": [1.0, -1.0]})
            self.assertEqual(record["method"]["functional"], "PBE")
            self.assertEqual(record["method"]["software"]["version"], "6.5.1")
            self.assertEqual(
                record["method"]["pseudopotentials"],
                [
                    {"symbol": "O", "title": "PAW_PBE O 08Apr2002", "zval": 6.0},
                    {"symbol": "H", "title": "PAW_PBE H 15Jun2001", "zval": 1.0},
                ],
            )
            self.assertIn("esp_tools", record["provenance"])
            self.assertAlmostEqual(record["energy"]["total_eV"], -10.0)
            self.assertAlmostEqual(record["energy"]["per_atom_eV"], -5.0)
            self.assertAlmostEqual(record["forces"]["max_norm_eV_per_Angstrom"], 0.3741657387)
            self.assertEqual(record["stress"]["pressure_kBar"], 200.0)
            self.assertAlmostEqual(
                record["stress"]["values_kBar"][0][2], 30.0
            )
            self.assertAlmostEqual(
                record["stress"]["values_eV_per_Angstrom3"][0][0],
                100.0 * unify.KBAR_TO_EV_PER_ANGSTROM3,
            )
            self.assertEqual(
                record["electrostatic_potential"]["points"],
                [[2.0, 2.0, 2.0, 1.0], [3.0, 3.0, 3.0, 2.0]],
            )

    def test_orca_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            structure = {
                "natoms": 1,
                "symbols": ["O"],
                "formula": "O",
                "positions_Angstrom": [[0.0, 0.0, 0.0]],
            }
            esp_path = tmp_path / "esp.json"
            write_json(
                esp_path,
                esp_payload(structure, [[0.0, 0.0, 1.0]], [0.5], vasp=False),
            )
            out = tmp_path / "water.out"
            out.write_text(
                " Program Version 6.1.0-f.0  -   RELEASE   -\n"
                " Your ORCA version has been built with support for libXC version: 7.0.0\n"
                " SCF CONVERGED AFTER 10 CYCLES\n"
                "FINAL SINGLE POINT ENERGY       -1.000000000000\n"
                " Total Charge           Charge          ....    0\n"
                " Multiplicity           Mult            ....    1\n"
                "                             ****ORCA TERMINATED NORMALLY****\n"
            )
            engrad = tmp_path / "water.engrad"
            engrad.write_text(
                "#\n# Number of atoms\n#\n 1\n"
                "#\n# The current total energy in Eh\n#\n -1.000000000000\n"
                "#\n# The current gradient in Eh/bohr\n#\n"
                " 0.010000000000\n 0.000000000000\n 0.000000000000\n"
            )
            inp = tmp_path / "water.inp"
            inp.write_text(
                "! wB97M-V def2-TZVPD AutoAux RIJCOSX VeryTightSCF\n"
                "%pal nprocs 4 end\n%maxcore 2000\n"
                "%method\n  Grid 5\nend\n"
                "* xyz 0 1\nO 0 0 0\n*\n"
            )
            output = tmp_path / "record.json"
            rc = unify.main([
                "orca",
                "--esp", str(esp_path),
                "--out", str(out),
                "--engrad", str(engrad),
                "--inp", str(inp),
                "--output", str(output),
            ])
            self.assertEqual(rc, 0)
            record = json.loads(output.read_text())
            validate_record(record)
            self.assertEqual(record["structure"]["charge"], 0)
            self.assertEqual(record["structure"]["multiplicity"], 1)
            self.assertEqual(record["method"]["functional"], "wB97M-V")
            self.assertEqual(record["method"]["basis"], "def2-TZVPD")
            self.assertEqual(record["method"]["aux_basis"], "AutoAux")
            self.assertEqual(record["method"]["software"]["version"], "6.1.0-f.0")
            self.assertEqual(record["method"]["software"]["libxc_version"], "7.0.0")
            self.assertIn("RIJCOSX", record["method"]["approximations"])
            self.assertAlmostEqual(record["energy"]["total_hartree"], -1.0)
            self.assertAlmostEqual(
                record["energy"]["total_eV"], -unify.HARTREE_EV
            )
            self.assertAlmostEqual(
                record["forces"]["values_eV_per_Angstrom"][0][0],
                -0.01 * unify.ORCA_GRAD_TO_EV_PER_ANGSTROM,
            )
            self.assertFalse(record["stress"]["available"])
            self.assertEqual(
                record["electrostatic_potential"]["points"],
                [[0.0, 0.0, 1.0, 0.5]],
            )


if __name__ == "__main__":
    unittest.main()
