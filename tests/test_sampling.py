from __future__ import annotations

from pathlib import Path
import sys
import unittest
import warnings

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common.sampling import (  # noqa: E402
    bin_indices,
    bin_quotas,
    nearest_atom_distances,
    sample_grid,
    structure_seed,
)

WATER = np.array([
    [0.000000, 0.000000, 0.119262],
    [0.000000, 0.763239, -0.477047],
    [0.000000, -0.763239, -0.477047],
])


class SamplingTests(unittest.TestCase):
    def test_ten_points_per_atom_and_quota(self):
        points, info = sample_grid(WATER)
        self.assertEqual(len(points), 30)
        self.assertEqual(info.npoints, 30)
        self.assertEqual(info.requested_bin_counts, (3, 12, 12, 3))
        self.assertEqual(info.actual_bin_counts, (3, 12, 12, 3))
        self.assertEqual(info.shortfall_bins, ())
        self.assertEqual(info.reallocated_points, 0)
        self.assertFalse(info.bin_reallocation)
        distances = nearest_atom_distances(points, WATER)
        np.testing.assert_array_equal(
            np.bincount(bin_indices(distances), minlength=4), [3, 12, 12, 3]
        )
        self.assertGreaterEqual(float(distances.min()), 0.25 - 1e-12)

    def test_reproducible_and_geometry_dependent(self):
        first, _ = sample_grid(WATER)
        second, _ = sample_grid(WATER)
        np.testing.assert_allclose(first, second)
        shifted = WATER + np.array([0.1, 0.0, 0.0])
        other, _ = sample_grid(shifted)
        self.assertFalse(np.allclose(first, other))
        self.assertNotEqual(structure_seed(WATER), structure_seed(shifted))

    def test_single_atom(self):
        points, info = sample_grid(np.zeros((1, 3)))
        self.assertEqual(len(points), 10)
        self.assertEqual(info.actual_bin_counts, (1, 4, 4, 1))

    def test_periodic_grid(self):
        lattice = np.eye(3) * 16.0
        points, info = sample_grid(WATER, lattice=lattice)
        self.assertEqual(len(points), 30)
        self.assertTrue(info.periodic)
        self.assertEqual(info.actual_bin_counts, (3, 12, 12, 3))
        self.assertTrue(np.all(points >= 0.0))
        self.assertTrue(np.all(points < 16.0))
        distances = nearest_atom_distances(points, WATER, lattice)
        np.testing.assert_array_equal(
            np.bincount(bin_indices(distances), minlength=4), [3, 12, 12, 3]
        )

    def test_dense_cell_shortfall_raises_by_default(self):
        lattice = np.eye(3) * 4.0
        with self.assertRaises(RuntimeError):
            sample_grid(np.zeros((1, 3)), lattice=lattice)

    def test_dense_cell_reallocates_with_warning(self):
        # A 4 A cell has no point >=5 A from the atom (minimum-image), so the
        # far-field quota (1 point) cannot be filled.  VASP asks for
        # reallocation: keep 10 points and warn instead of failing.
        lattice = np.eye(3) * 4.0
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            points, info = sample_grid(
                np.zeros((1, 3)),
                lattice=lattice,
                shortfall_policy="reallocate",
            )
        self.assertEqual(len(points), 10)
        self.assertEqual(sum(info.actual_bin_counts), 10)
        self.assertEqual(info.shortfall_bins, (3,))
        self.assertGreaterEqual(info.reallocated_points, 1)
        self.assertTrue(info.bin_reallocation)
        # A far-field shortfall must go to [1.5, 5) (bin 2) first.
        self.assertEqual(info.reallocated_bin_counts, (0, 0, 1, 0))
        self.assertEqual(info.actual_bin_counts, (1, 4, 5, 0))
        self.assertTrue(any(issubclass(w.category, RuntimeWarning) for w in caught))
        self.assertIn("[1.5, 5)", str(caught[0].message))
        distances = nearest_atom_distances(points, np.zeros((1, 3)), lattice)
        self.assertGreaterEqual(float(distances.min()), 0.25 - 1e-12)

    def test_quota_exact_sum(self):
        for total in (1, 3, 7, 10, 30, 101):
            quota = bin_quotas(total)
            self.assertEqual(sum(quota), total)
            self.assertTrue(all(q >= 0 for q in quota))


if __name__ == "__main__":
    unittest.main()
