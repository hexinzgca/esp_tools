from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vasp.extract_esp import main as extract_main  # noqa: E402
from vasp.locpot import interpolate_periodic, read_locpot  # noqa: E402


SHAPE = (4, 5, 6)
LATTICE = np.diag([4.0, 5.0, 6.0])


def analytic(frac: np.ndarray) -> np.ndarray:
    return 1.0 + 2.0 * frac[:, 0] + 3.0 * frac[:, 1] + 4.0 * frac[:, 2]


def write_locpot(path: Path, values: np.ndarray, blocks: int = 1) -> None:
    """Write a synthetic LOCPOT; values has shape (nx, ny, nz) in eV."""
    nx, ny, nz = values.shape
    lines = [
        "synthetic LOCPOT",
        "1.0",
        "4.0 0.0 0.0",
        "0.0 5.0 0.0",
        "0.0 0.0 6.0",
        "H",
        "1",
        "Cartesian",
        "0.0 0.0 0.0",
        "",
    ]
    for _ in range(blocks):
        lines.append(f"{nx} {ny} {nz}")
        # VASP ordering: x fastest.
        flat = values.transpose(2, 1, 0).ravel()
        for start in range(0, len(flat), 5):
            lines.append(" ".join(f"{v:.12e}" for v in flat[start:start + 5]))
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


class LocpotTests(unittest.TestCase):
    def test_ordering_and_interpolation(self):
        frac_grid = np.stack(
            np.meshgrid(
                np.arange(SHAPE[0]) / SHAPE[0],
                np.arange(SHAPE[1]) / SHAPE[1],
                np.arange(SHAPE[2]) / SHAPE[2],
                indexing="ij",
            ),
            axis=-1,
        )
        values = analytic(frac_grid.reshape(-1, 3)).reshape(SHAPE)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "LOCPOT"
            write_locpot(path, values)
            locpot = read_locpot(path)
            self.assertEqual(locpot.grid_shape, SHAPE)
            self.assertEqual(locpot.n_blocks, 1)
            np.testing.assert_allclose(locpot.potential, values)

            frac = np.array([
                [0.12, 0.22, 0.32],
                [0.31, 0.41, 0.53],
                [0.58, 0.11, 0.62],
            ])
            points = frac @ LATTICE
            got = interpolate_periodic(points, LATTICE, locpot.potential)
            np.testing.assert_allclose(got, analytic(frac), atol=1e-12)

    def test_delta_peak_ordering(self):
        values = np.zeros(SHAPE)
        values[1, 2, 3] = 7.0
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "LOCPOT"
            write_locpot(path, values)
            locpot = read_locpot(path)
            frac = np.array([[1 / 4, 2 / 5, 3 / 6]])
            got = interpolate_periodic(frac @ LATTICE, LATTICE, locpot.potential)
            np.testing.assert_allclose(got, [7.0], atol=1e-12)

    def test_two_blocks_are_averaged(self):
        values = np.ones(SHAPE)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "LOCPOT"
            write_locpot(path, values, blocks=2)
            locpot = read_locpot(path)
            self.assertEqual(locpot.n_blocks, 2)
            self.assertAlmostEqual(locpot.spin_block_max_abs_diff_eV, 0.0)
            np.testing.assert_allclose(locpot.potential, 1.0)

    def test_cli_dense_cell_reallocates_and_falls_back(self):
        # 4 x 5 x 6 A cell: no point is >=5 A from the single atom under the
        # minimum-image convention.  The far-field quota must be reallocated
        # to [1.5, 5), and the vacuum reference must fall back to LOCPOT zero
        # instead of aborting.
        values = np.zeros(SHAPE)
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            locpot_path = tmp_path / "LOCPOT"
            write_locpot(locpot_path, values)
            output = tmp_path / "esp.json"
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                rc = extract_main([
                    "--locpot", str(locpot_path),
                    "--sample-mode", "interpolate",
                    "--output", str(output),
                ])
            self.assertEqual(rc, 0)
            text = stderr.getvalue()
            self.assertIn("WARNING", text)
            self.assertIn("[1.5, 5)", text)
            report = json.loads(output.read_text())
            self.assertEqual(report["npoints"], 10)
            self.assertEqual(report["sampler"]["shortfall_bins"], [3])
            self.assertEqual(report["sampler"]["reallocated_bin_counts"], [0, 0, 1, 0])
            self.assertEqual(report["sampler"]["actual_bin_counts"], [1, 4, 5, 0])
            reference = report["locpot"]["reference"]
            self.assertEqual(reference["mode"], "none")
            self.assertEqual(reference["requested_mode"], "vacuum")
            self.assertTrue(reference["fallback"])

    def test_cli_dft_grid_mode(self):
        # 64^3 over 16 A => 0.25 A grid spacing.  Direct mode must return
        # exact grid points and values, with no interpolation.
        lattice = np.eye(3) * 16.0
        shape = (64, 64, 64)
        positions = np.array([
            [8.0, 8.0, 8.119262],
            [8.0, 8.763239, 7.522953],
            [8.0, 7.236761, 7.522953],
        ])
        i, j, k = np.indices(shape)
        values = (i + 2.0 * j + 3.0 * k).astype(float)
        lines = [
            "synthetic water LOCPOT for direct grid mode",
            "1.0",
            "16.0 0.0 0.0",
            "0.0 16.0 0.0",
            "0.0 0.0 16.0",
            "O H",
            "1 2",
            "Cartesian",
        ]
        lines += [f"{row[0]:.10f} {row[1]:.10f} {row[2]:.10f}" for row in positions]
        lines.append("")
        lines.append("64 64 64")
        flat = values.transpose(2, 1, 0).ravel()
        for start in range(0, len(flat), 5):
            lines.append(" ".join(f"{v:.12e}" for v in flat[start:start + 5]))
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            locpot_path = tmp_path / "LOCPOT"
            locpot_path.write_text("\n".join(lines) + "\n")
            output = tmp_path / "esp.json"
            rc = extract_main([
                "--locpot", str(locpot_path),
                "--output", str(output),
            ])
            self.assertEqual(rc, 0)
            report = json.loads(output.read_text())
            self.assertEqual(report["sampler"]["sampling_mode"], "dft-grid")
            self.assertEqual(report["npoints"], 30)
            self.assertEqual(report["sampler"]["requested_bin_counts"], [3, 12, 12, 3])
            self.assertEqual(report["sampler"]["actual_bin_counts"], [3, 12, 12, 3])
            self.assertEqual(report["locpot"]["interpolation"], "none (direct LOCPOT grid values)")
            points = np.asarray(report["esp"]["positions_Angstrom"])
            frac = points @ np.linalg.inv(lattice)
            np.testing.assert_allclose(frac * shape, np.round(frac * shape), atol=1e-10)
            flat_indices = report["esp"]["grid_indices_ijk"]
            raw = np.asarray(report["esp"]["phi_locpot_eV"])
            for (ii, jj, kk), q in zip(flat_indices, raw):
                self.assertAlmostEqual(q, values[ii, jj, kk])

    def test_cli_end_to_end(self):
        # Water in a 16 A box; a smooth synthetic potential is enough to test
        # sampling, interpolation, reference selection, and JSON writing.
        lattice = np.eye(3) * 16.0
        positions = np.array([
            [8.0, 8.0, 8.119262],
            [8.0, 8.763239, 7.522953],
            [8.0, 7.236761, 7.522953],
        ])
        frac = positions @ np.linalg.inv(lattice)
        shape = (16, 16, 16)
        grid = np.stack(
            np.meshgrid(
                np.arange(shape[0]) / shape[0],
                np.arange(shape[1]) / shape[1],
                np.arange(shape[2]) / shape[2],
                indexing="ij",
            ),
            axis=-1,
        )
        values = np.sin(2 * np.pi * grid[..., 0]) + np.cos(2 * np.pi * grid[..., 1])
        lines = [
            "synthetic water LOCPOT",
            "1.0",
            "16.0 0.0 0.0",
            "0.0 16.0 0.0",
            "0.0 0.0 16.0",
            "O H",
            "1 2",
            "Cartesian",
        ]
        lines += [f"{row[0]:.10f} {row[1]:.10f} {row[2]:.10f}" for row in positions]
        lines.append("")
        lines.append("16 16 16")
        flat = values.transpose(2, 1, 0).ravel()
        for start in range(0, len(flat), 5):
            lines.append(" ".join(f"{v:.12e}" for v in flat[start:start + 5]))
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            locpot_path = tmp_path / "LOCPOT"
            locpot_path.write_text("\n".join(lines) + "\n")
            output = tmp_path / "esp.json"
            rc = extract_main([
                "--locpot", str(locpot_path),
                "--sample-mode", "interpolate",
                "--output", str(output),
            ])
            self.assertEqual(rc, 0)
            report = json.loads(output.read_text())
            self.assertEqual(report["npoints"], 30)
            self.assertEqual(report["sampler"]["sampling_mode"], "interpolate")
            self.assertEqual(report["sampler"]["actual_bin_counts"], [3, 12, 12, 3])
            self.assertEqual(len(report["esp"]["positions_Angstrom"]), 30)
            self.assertEqual(len(report["esp"]["phi_eV_per_e"]), 30)
            self.assertTrue(np.isfinite(report["esp"]["phi_eV_per_e"]).all())
            raw = np.asarray(report["esp"]["phi_locpot_eV"])
            phi = np.asarray(report["esp"]["phi_eV_per_e"])
            locpot_mean = report["locpot"]["reference"]["locpot_mean_eV"]
            np.testing.assert_allclose(phi, locpot_mean - raw, atol=1e-12)
            self.assertAlmostEqual(
                float(report["locpot"]["reference"]["value_eV"]), -locpot_mean
            )


if __name__ == "__main__":
    unittest.main()
