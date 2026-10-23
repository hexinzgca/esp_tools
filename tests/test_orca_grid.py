from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from orca.extract_esp import main as extract_main  # noqa: E402


class OrcaGridTests(unittest.TestCase):
    def test_grid_only_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "grid.json"
            rc = extract_main([
                "--xyz", str(ROOT / "orca" / "water.xyz"),
                "--output", str(output),
            ])
            self.assertEqual(rc, 0)
            report = json.loads(output.read_text())
            self.assertEqual(report["npoints"], 30)
            self.assertEqual(report["sampler"]["requested_bin_counts"], [3, 12, 12, 3])
            self.assertEqual(report["sampler"]["actual_bin_counts"], [3, 12, 12, 3])
            self.assertEqual(len(report["esp"]["positions_Angstrom"]), 30)
            self.assertNotIn("phi_molecular_eV_per_e", report["esp"])
            self.assertNotIn("provenance", report)

    def test_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "grid.json"
            output.write_text("{}")
            with self.assertRaises(SystemExit):
                extract_main([
                    "--xyz", str(ROOT / "orca" / "water.xyz"),
                    "--output", str(output),
                ])


if __name__ == "__main__":
    unittest.main()
