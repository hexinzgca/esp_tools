#!/usr/bin/env python3
"""Extract ESP on a 10 x natoms grid from a VASP LOCPOT file.

The grid is sampled with the shared 1:4:4:1 nearest-atom-distance scheme
(see ``common/sampling.py``).  Only those grid points are interpolated from
LOCPOT; the full NGXF x NGYF x NGZF grid is never written out.

Typical use::

    python vasp/extract_esp.py --locpot LOCPOT --output esp.json

Add ``--poscar POSCAR`` to verify that LOCPOT and POSCAR describe the same
geometry.  The default reference is the mean potential of the far-field bin
(>= 5 A from every atom), which estimates the vacuum level for a neutral
molecule in a sufficiently large box.  Use ``--reference none`` for the
LOCPOT zero or ``--reference <eV>`` for an explicit electrostatic offset.

SIGN
    LOCPOT (``LVHAR``) contains the potential energy of an electron, which
    is negative near nuclei.  The electrostatic potential reported by
    ``orca_vpot`` and by this tool is therefore ``phi = -LOCPOT`` (e = 1 in
    atomic units).  The JSON keeps the raw ``phi_locpot_eV`` for inspection.
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

from common.sampling import (  # noqa: E402
    ESP_BIN_EDGES,
    ESP_BIN_RATIO,
    ESP_FAR_SPAN,
    ESP_GRID_SEED,
    ESP_POINTS_PER_ATOM,
    bin_indices,
    nearest_atom_distances,
    sample_grid,
)
from vasp.locpot import (  # noqa: E402
    interpolate_periodic,
    read_locpot,
    read_poscar,
    sha256_file,
)


def parse_reference(value: str) -> tuple[str, float | None]:
    text = value.strip().lower()
    if text in {"vacuum", "far", "far-mean"}:
        return "vacuum", None
    if text in {"none", "raw"}:
        return "none", 0.0
    try:
        return "value", float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "reference must be 'vacuum', 'none', or an electrostatic potential in eV"
        ) from exc


def formula(symbols: list[str]) -> str | None:
    if not symbols:
        return None
    counts: dict[str, int] = {}
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
    return "".join(
        symbol + (str(count) if count > 1 else "")
        for symbol, count in counts.items()
    )


def compare_geometry(locpot, poscar) -> None:
    if not np.allclose(locpot.lattice, poscar.lattice, atol=1e-4, rtol=0):
        raise ValueError(
            "POSCAR lattice differs from the lattice in LOCPOT; "
            "refusing to mix geometry and potential"
        )
    if locpot.positions.shape != poscar.positions.shape or not np.allclose(
        locpot.positions, poscar.positions, atol=1e-4, rtol=0
    ):
        raise ValueError(
            "POSCAR atomic positions differ from those in LOCPOT; "
            "refusing to mix geometry and potential"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locpot", required=True, help="VASP LOCPOT file (LVHAR=.TRUE.)")
    parser.add_argument("--poscar", help="Optional POSCAR for a geometry cross-check")
    parser.add_argument("--output", required=True, help="Output JSON; never overwritten unless --force")
    parser.add_argument(
        "--points-per-atom",
        type=int,
        default=ESP_POINTS_PER_ATOM,
        help=f"Probe points per atom (default {ESP_POINTS_PER_ATOM}; total = this x natoms)",
    )
    parser.add_argument("--seed", type=int, default=ESP_GRID_SEED, help="Base sampling seed")
    parser.add_argument(
        "--reference",
        default="vacuum",
        help=(
            "'vacuum' (mean of the >=5 A bin, converted to electrostatic sign), "
            "'none' (LOCPOT zero), or an explicit electrostatic potential in eV"
        ),
    )
    parser.add_argument("--far-span", type=float, default=ESP_FAR_SPAN, help="Far-shell span in A (periodic mode ignores it)")
    parser.add_argument("--min-distance", type=float, default=ESP_BIN_EDGES[0], help="Minimum distance to any atom in A")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file")
    args = parser.parse_args(argv)

    output = Path(args.output).expanduser()
    if output.exists() and not args.force:
        parser.error(f"output already exists: {output}; choose another name or pass --force")
    reference_mode, reference_value_arg = parse_reference(args.reference)

    locpot = read_locpot(args.locpot)
    if locpot.n_blocks > 1:
        print(
            f"LOCPOT contains {locpot.n_blocks} potential blocks; "
            f"using their mean (max deviation {locpot.spin_block_max_abs_diff_eV:.3e} eV)"
        )
    if args.poscar:
        compare_geometry(locpot, read_poscar(args.poscar))

    points, sampler = sample_grid(
        locpot.positions,
        lattice=locpot.lattice,
        points_per_atom=args.points_per_atom,
        bin_ratio=ESP_BIN_RATIO,
        seed=args.seed,
        min_atom_distance=args.min_distance,
        far_span=args.far_span,
        shortfall_policy="reallocate",
        warn_on_reallocation=False,  # the CLI prints a structured warning below
    )
    if sampler.bin_reallocation:
        print(
            "WARNING: the requested bin distribution could not be filled in "
            f"this periodic cell; bins {list(sampler.shortfall_bins)} were short. "
            f"Reallocated {sampler.reallocated_points} point(s) "
            f"(per bin {list(sampler.reallocated_bin_counts)}; "
            "far-field deficit goes to [1.5, 5) first): "
            f"requested {list(sampler.requested_bin_counts)}, "
            f"actual {list(sampler.actual_bin_counts)}.",
            file=sys.stderr,
        )
    raw = interpolate_periodic(points, locpot.lattice, locpot.potential)

    distances = nearest_atom_distances(points, locpot.positions, locpot.lattice)
    bins = bin_indices(distances)
    # LOCPOT stores the potential energy of an electron.  The electrostatic
    # potential (the quantity orca_vpot reports and the one normally called
    # ESP) is phi = -LOCPOT/e, with e = 1 in atomic units.
    requested_reference_mode = reference_mode
    reference_fallback = False
    if reference_mode == "vacuum":
        far = raw[bins == 3]
        if not len(far):
            # No far-field region (dense/small cell).  Do not abort the
            # extraction; fall back to the LOCPOT zero and say so.
            print(
                "WARNING: --reference vacuum requested but the sampled grid has "
                "no point in the >=5 A far bin; falling back to --reference none "
                "(LOCPOT zero).",
                file=sys.stderr,
            )
            reference_mode = "none"
            reference_fallback = True
            locpot_reference = 0.0
            reference_value = 0.0
        else:
            locpot_reference = float(far.mean())
            reference_value = -locpot_reference
    elif reference_mode == "none":
        locpot_reference = 0.0
        reference_value = 0.0
    else:
        reference_value = float(reference_value_arg)
        locpot_reference = -reference_value

    phi = -raw - reference_value
    report = {
        "format": "esp_tools/vasp-esp-v1",
        "method": (
            "VASP LOCPOT (LVHAR: electron potential energy V_ionic + V_Hartree) "
            "trilinear interpolation; phi = -LOCPOT"
        ),
        "units": {
            "positions": "Angstrom",
            "potential": "eV/e (= V)",
            "locpot": "eV",
        },
        "structure": {
            "natoms": len(locpot.positions),
            "symbols": locpot.symbols,
            "counts": locpot.counts,
            "formula": formula(locpot.symbols),
            "lattice_Angstrom": locpot.lattice.tolist(),
            "positions_Angstrom": locpot.positions.tolist(),
        },
        "sampler": sampler.as_dict(),
        "locpot": {
            "path": str(Path(args.locpot).expanduser()),
            "sha256": sha256_file(args.locpot),
            "grid_shape": list(locpot.grid_shape),
            "n_blocks": locpot.n_blocks,
            "spin_block_max_abs_diff_eV": locpot.spin_block_max_abs_diff_eV,
            "reference": {
                "mode": reference_mode,
                "requested_mode": requested_reference_mode,
                "fallback": reference_fallback,
                "value_eV": reference_value,
                "locpot_mean_eV": locpot_reference,
                "description": (
                    "no far-field points; fell back to the LOCPOT zero"
                    if reference_fallback
                    else "mean of the far-field (>=5 A) sampled points, "
                    "converted from LOCPOT electron potential energy to ESP"
                    if reference_mode == "vacuum"
                    else "LOCPOT zero" if reference_mode == "none"
                    else "explicit electrostatic-potential offset"
                ),
            },
            "interpolation": "periodic trilinear on the LOCPOT grid",
            "sign_convention": (
                "LOCPOT = potential energy of an electron; "
                "phi_eV_per_e = -LOCPOT_eV - reference_eV"
            ),
        },
        "esp": {
            "positions_Angstrom": points.tolist(),
            "phi_locpot_eV": raw.tolist(),
            "phi_eV_per_e": phi.tolist(),
        },
        "npoints": int(len(points)),
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.force else "x"
    with output.open(mode) as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")

    print(
        f"Wrote {len(points)} points ({args.points_per_atom}/atom) to {output}\n"
        f"  bin counts (requested/actual): {list(sampler.requested_bin_counts)} / {list(sampler.actual_bin_counts)}\n"
        f"  reference (electrostatic): {reference_mode} = {reference_value:.6f} eV\n"
        f"  phi range: {phi.min():.6f} .. {phi.max():.6f} eV"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
