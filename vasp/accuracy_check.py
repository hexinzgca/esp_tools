#!/usr/bin/env python3
"""Compare coarse LOCPOT trilinear interpolation with fine direct grid values.

Inputs:

- ``--coarse-locpot``: coarse LOCPOT (the interpolation source).
- ``--fine-locpot``: fine LOCPOT (the reference grid).
- ``--fine-esp``: a ``dft-grid`` ESP JSON produced from the fine LOCPOT; its
  points are exact fine-grid coordinates and its ``phi_locpot_eV`` values are
  exact fine-grid values.

The tool interpolates the coarse LOCPOT at the fine-grid points and reports the
error against the direct fine-grid values, overall and per distance bin.  It
also checks that trilinear interpolation at the fine grid's own nodes returns
the direct values exactly.

Example::

    python vasp/accuracy_check.py \
        --coarse-locpot coarse/LOCPOT \
        --fine-locpot fine/LOCPOT \
        --fine-esp fine/grid.json \
        --output accuracy.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common.sampling import bin_indices, nearest_atom_distances  # noqa: E402
from vasp.locpot import interpolate_periodic, read_locpot  # noqa: E402

BIN_LABELS = ("0.25-0.5", "0.5-1.5", "1.5-5.0", "5.0-inf")


def _stats(error: np.ndarray) -> dict:
    return {
        "n": int(len(error)),
        "mean_abs_eV": float(np.mean(np.abs(error))),
        "max_abs_eV": float(np.max(np.abs(error))),
        "rms_eV": float(np.sqrt(np.mean(error ** 2))),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse-locpot", required=True)
    parser.add_argument("--fine-locpot", required=True)
    parser.add_argument("--fine-esp", required=True, help="fine dft-grid ESP JSON")
    parser.add_argument("--output", help="optional JSON report")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    coarse = read_locpot(args.coarse_locpot)
    fine = read_locpot(args.fine_locpot)
    if coarse.positions.shape != fine.positions.shape or not np.allclose(
        coarse.positions, fine.positions, atol=1e-5, rtol=0
    ):
        parser.error("coarse and fine LOCPOT have different atomic positions")
    if not np.allclose(coarse.lattice, fine.lattice, atol=1e-5, rtol=0):
        parser.error("coarse and fine LOCPOT have different lattices")

    report = json.loads(Path(args.fine_esp).read_text())
    if report.get("sampler", {}).get("sampling_mode") != "dft-grid":
        parser.error("--fine-esp must be produced with --sample-mode dft-grid")
    points = np.asarray(report["esp"]["positions_Angstrom"], dtype=float)
    fine_direct = np.asarray(report["esp"]["phi_locpot_eV"], dtype=float)

    coarse_interp = interpolate_periodic(points, coarse.lattice, coarse.potential)
    fine_interp = interpolate_periodic(points, fine.lattice, fine.potential)
    self_check = float(np.max(np.abs(fine_interp - fine_direct)))

    bins = bin_indices(
        nearest_atom_distances(points, fine.positions, fine.lattice)
    )
    error = coarse_interp - fine_direct
    rows = [{"bin": BIN_LABELS[b], **_stats(error[bins == b])}
            for b in range(4) if np.any(bins == b)]

    result = {
        "format": "esp_tools/accuracy-check-v1",
        "units": "eV",
        "coarse_locpot": str(Path(args.coarse_locpot).resolve()),
        "fine_locpot": str(Path(args.fine_locpot).resolve()),
        "fine_esp": str(Path(args.fine_esp).resolve()),
        "grid_spacing_Angstrom": {
            "coarse": (
                np.linalg.norm(coarse.lattice, axis=1)
                / np.asarray(coarse.grid_shape, dtype=float)
            ).tolist(),
            "fine": (
                np.linalg.norm(fine.lattice, axis=1)
                / np.asarray(fine.grid_shape, dtype=float)
            ).tolist(),
        },
        "interpolation_at_fine_grid_nodes_max_abs_eV": self_check,
        "coarse_interpolation_error_vs_fine_direct": {
            "all": _stats(error),
            "per_bin": rows,
        },
    }

    print(
        f"fine-grid self-check (interp - direct) max abs: {self_check:.3e} eV\n"
        f"{'bin':>10} {'n':>4} {'mean_abs':>12} {'max_abs':>12} {'rms':>12}"
    )
    for row in rows:
        print(
            f"{row['bin']:>10} {row['n']:4d} {row['mean_abs_eV']:12.6f} "
            f"{row['max_abs_eV']:12.6f} {row['rms_eV']:12.6f}"
        )
    all_stats = result["coarse_interpolation_error_vs_fine_direct"]["all"]
    print(
        f"{'all':>10} {all_stats['n']:4d} {all_stats['mean_abs_eV']:12.6f} "
        f"{all_stats['max_abs_eV']:12.6f} {all_stats['rms_eV']:12.6f}"
    )

    if args.output:
        output = Path(args.output).expanduser()
        if output.exists() and not args.force:
            parser.error(f"output already exists: {output}; pass --force")
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w" if args.force else "x") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
