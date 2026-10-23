#!/usr/bin/env python3
"""Sample a 10 x natoms grid and extract ORCA ESP with ``orca_vpot``.

SAMPLING
    ``--xyz water.xyz`` samples the shared 1:4:4:1 nearest-atom-distance
    grid.  Only ``10 * natoms`` points are produced; no full cube is written.
    The innermost bin starts at 0.25 A, so points never sit on a nucleus.

EVALUATION
    ``--gbw water.gbw`` evaluates the molecular electrostatic potential
    (nuclei + SCF electrons) with the matching SCF density.  ORCA 6.1 stores
    densities in ``water.densities``; the tool discovers the ``.scfp`` name
    through ``orca_vpot`` instead of guessing it.  A native ORCA installation
    can be used with ``--orca-home``, or the calculation can run inside an
    Apptainer image with ``--orca-container``.

EXAMPLES
    # Native ORCA
    python orca/extract_esp.py --xyz orca/water.xyz --gbw water.gbw \
        --orca-home /opt/orca/bin --output water_esp.json

    # ORCA in the project's Apptainer image
    python orca/extract_esp.py --xyz orca/water.xyz --gbw water.gbw \
        --orca-container /mnt/e/apptainer_sifs/build/orca-6.1.sif \
        --output water_esp.json

    # Sample only (no SCF wavefunction needed)
    python orca/extract_esp.py --xyz orca/water.xyz --output water_grid.json

    # Reuse the same probe points for several charge states
    python orca/extract_esp.py --points water_grid.json --gbw anion.gbw \
        --orca-container /path/orca-6.1.sif --output anion_esp.json
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common.sampling import (  # noqa: E402
    BOHR_ANGSTROM,
    ESP_BIN_EDGES,
    ESP_BIN_RATIO,
    ESP_FAR_SPAN,
    ESP_GRID_SEED,
    ESP_POINTS_PER_ATOM,
    HARTREE_EV,
    read_xyz,
    sample_grid,
)


class OrcaBackend:
    """Run ``orca_vpot`` natively or inside an Apptainer image."""

    def __init__(
        self,
        home: str | None = None,
        container: str | None = None,
        apptainer: str = "apptainer",
    ):
        if home and container:
            raise ValueError("use either --orca-home or --orca-container, not both")
        self.home: Path | None = None
        self.container: Path | None = None
        self.apptainer = apptainer

        if container:
            self.container = Path(container).expanduser().resolve(strict=True)
            if not self.container.is_file():
                raise FileNotFoundError(f"container image not found: {self.container}")
        elif home:
            self.home = Path(home).expanduser().resolve(strict=True)
            if not self.home.is_dir():
                raise NotADirectoryError(f"ORCA home is not a directory: {self.home}")
        else:
            binary = shutil.which("orca_vpot") or shutil.which("orca")
            if not binary:
                raise FileNotFoundError(
                    "cannot locate orca_vpot; pass --orca-home or --orca-container"
                )
            self.home = Path(binary).resolve().parent

    def tool(self, name: str) -> str:
        if self.container:
            return name
        assert self.home is not None
        candidates = [self.home / name, self.home / "bin" / name]
        for candidate in candidates:
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
        raise FileNotFoundError(f"executable not found: {candidates[0]}")

    def env(self) -> dict[str, str]:
        env = os.environ.copy()
        if self.home is not None:
            env["PATH"] = str(self.home) + os.pathsep + env.get("PATH", "")
            env["ORCA_HOME"] = str(self.home)
        env.setdefault("OMP_NUM_THREADS", "1")
        return env

    def describe(self) -> dict:
        if self.container:
            return {"mode": "apptainer", "container": str(self.container)}
        return {"mode": "native", "home": str(self.home)}

    def run(
        self,
        argv: list[str],
        work: Path,
        *,
        timeout: float,
    ) -> tuple[subprocess.CompletedProcess, list[str]]:
        if self.container:
            command = [
                self.apptainer,
                "exec",
                "--bind",
                f"{work}:{work}",
                str(self.container),
            ] + argv
        else:
            command = [self.tool(argv[0])] + argv[1:]
        proc = subprocess.run(
            command,
            cwd=work,
            env=self.env(),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return proc, command


def _density_names_from_listing(text: str) -> list[str]:
    names = []
    for line in text.splitlines():
        match = re.match(r"^\s*\d+\s*:\s*(\S+)\s*$", line)
        if match:
            names.append(match.group(1))
    return names


def locate_density(
    backend: OrcaBackend,
    gbw: Path,
    work: Path,
    prefix: str,
) -> tuple[str, str, list[str]]:
    """Return (density_name, listing_or_path, command)."""
    density_container = gbw.with_suffix(".densities")
    if density_container.is_file():
        proc, command = backend.run(
            ["orca_vpot", density_container.name], work, timeout=60
        )
        log = work / f"{prefix}_density_list.log"
        log.write_text(proc.stdout + proc.stderr)
        names = _density_names_from_listing(proc.stdout)
        if proc.returncode != 0:
            raise RuntimeError(
                f"orca_vpot could not list {density_container}; see {log}"
            )
        scfp = [name for name in names if name.endswith(".scfp")]
        if len(scfp) > 1 and f"{prefix}.scfp" in scfp:
            scfp = [f"{prefix}.scfp"]
        if len(scfp) != 1:
            raise RuntimeError(
                f"cannot uniquely locate the SCF density in {density_container}: "
                f"{names}; see {log}"
            )
        return scfp[0], str(density_container), command

    # ORCA versions that keep a separate .scfp file.
    scfp = work / f"{prefix}.scfp"
    if scfp.is_file():
        return scfp.name, str(scfp), []
    raise FileNotFoundError(
        f"neither {density_container} nor {scfp} exists; "
        "run the ORCA SCF calculation first"
    )


def evaluate_points(
    backend: OrcaBackend,
    gbw: Path,
    points_angstrom: np.ndarray,
    work: Path,
    prefix: str,
) -> tuple[np.ndarray, dict]:
    """Evaluate nuclear + electronic ESP (Hartree/e) with this GBW density."""
    gbw = gbw.resolve(strict=True)
    work = work.resolve()
    points = np.asarray(points_angstrom, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("ESP points must have shape (n, 3) and be finite")

    density_name, density_source, list_command = locate_density(
        backend, gbw, work, prefix
    )
    stamp = time.time_ns()
    grid = work / f"{prefix}_grid_{stamp}_bohr.txt"
    output = work / f"{prefix}_esp_{stamp}_bohr.txt"
    np.savetxt(
        grid,
        points / BOHR_ANGSTROM,
        fmt="%.14f",
        header=str(len(points)),
        comments="",
    )
    if output.exists():
        output.unlink()

    proc, command = backend.run(
        ["orca_vpot", gbw.name, density_name, grid.name, output.name],
        work,
        timeout=300,
    )
    log = work / f"{prefix}_vpot_{stamp}.log"
    log.write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0 or not output.exists():
        raise RuntimeError(f"orca_vpot failed; see {log}")

    raw = np.loadtxt(output, skiprows=1, ndmin=2)
    if raw.shape != (len(points), 4) or not np.isfinite(raw).all():
        raise ValueError(f"invalid orca_vpot output shape/values: {raw.shape}")
    np.testing.assert_allclose(
        raw[:, :3] * BOHR_ANGSTROM, points, atol=1e-10, rtol=0
    )

    provenance = {
        "orca": backend.describe(),
        "gbw": str(gbw),
        "density": density_name,
        "density_source": density_source,
        "density_list_command": list_command,
        "vpot_command": command,
        "grid_bohr": str(grid),
        "raw_output_bohr_hartree": str(output),
        "log": str(log),
        "grid_input_units": "Bohr",
        "raw_output_units": "Bohr, Hartree/e",
    }
    return raw[:, 3], provenance


def external_components(
    points_angstrom: np.ndarray,
    background_charges: np.ndarray,
    electric_field: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Background point-charge and uniform-field potentials in eV/e."""
    points_bohr = np.asarray(points_angstrom, dtype=np.float64) / BOHR_ANGSTROM
    charges = np.asarray(background_charges, dtype=np.float64)
    phi_background_ha = np.zeros(len(points_bohr), dtype=np.float64)
    if len(charges):
        charge_positions_bohr = charges[:, :3] / BOHR_ANGSTROM
        distances = np.linalg.norm(
            points_bohr[:, None, :] - charge_positions_bohr[None, :, :], axis=2
        )
        if np.any(distances < 1e-10):
            raise ValueError("ESP grid intersects a background point charge")
        phi_background_ha = (charges[None, :, 3] / distances).sum(axis=1)

    field = np.asarray(electric_field, dtype=np.float64)
    phi_field_ha = -points_bohr @ field  # gauge origin = (0, 0, 0)
    return phi_background_ha * HARTREE_EV, phi_field_ha * HARTREE_EV


def read_points(path: str | Path) -> tuple[np.ndarray, dict | None]:
    path = Path(path)
    metadata: dict | None = None
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text())
        grid = data.get("esp", data)
        if not isinstance(grid, dict) or "positions_Angstrom" not in grid:
            raise ValueError("JSON must contain esp.positions_Angstrom")
        points = np.asarray(grid["positions_Angstrom"], dtype=np.float64)
        metadata = data.get("sampler")
    else:
        points = np.loadtxt(path, ndmin=2)
    if (
        points.ndim != 2
        or points.shape[1] != 3
        or not len(points)
        or not np.isfinite(points).all()
    ):
        raise ValueError("expected a nonempty finite Nx3 point set in Angstrom")
    return points, metadata


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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--xyz", help="One XYZ frame in Angstrom; sample a new grid")
    source.add_argument("--points", help="Reuse a grid JSON or plain Nx3 Angstrom text")

    parser.add_argument("--gbw", help="Existing ORCA .gbw; omit to only export the grid")
    parser.add_argument("--orca-home", help="Directory containing orca_vpot (native ORCA)")
    parser.add_argument("--orca-container", help="Apptainer/Singularity image containing ORCA")
    parser.add_argument("--apptainer", default="apptainer", help="Apptainer executable (default: apptainer)")
    parser.add_argument("--output", required=True, help="Output JSON; never overwritten unless --force")
    parser.add_argument("--background", help="JSON list of [x_A, y_A, z_A, q_e] rows")
    parser.add_argument(
        "--field",
        nargs=3,
        type=float,
        default=[0.0, 0.0, 0.0],
        metavar=("EX", "EY", "EZ"),
        help="Uniform electric field in atomic units; gauge origin (0,0,0)",
    )
    parser.add_argument(
        "--points-per-atom",
        type=int,
        default=ESP_POINTS_PER_ATOM,
        help=f"Probe points per atom (default {ESP_POINTS_PER_ATOM}; total = this x natoms)",
    )
    parser.add_argument("--seed", type=int, default=ESP_GRID_SEED, help="Base sampling seed")
    parser.add_argument("--far-span", type=float, default=ESP_FAR_SPAN, help="Far-shell span in A")
    parser.add_argument("--min-distance", type=float, default=ESP_BIN_EDGES[0], help="Minimum distance to any atom in A")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file")
    args = parser.parse_args(argv)

    output = Path(args.output).expanduser()
    if output.exists() and not args.force:
        parser.error(f"output already exists: {output}; choose another name or pass --force")
    if args.points and not args.gbw:
        print("note: --points without --gbw only re-exports the grid", file=sys.stderr)

    background = np.zeros((0, 4), dtype=np.float64)
    if args.background:
        background = np.asarray(json.loads(Path(args.background).read_text()), dtype=np.float64)
        if background.size == 0:
            background = np.zeros((0, 4), dtype=np.float64)
        if background.ndim != 2 or background.shape[1] != 4 or not np.isfinite(background).all():
            parser.error("--background must be a JSON array of [x, y, z, q] rows")
    field = np.asarray(args.field, dtype=np.float64)
    if not np.isfinite(field).all():
        parser.error("--field values must be finite")

    symbols: list[str] = []
    if args.xyz:
        symbols, source_positions = read_xyz(args.xyz)
        points, sampler = sample_grid(
            source_positions,
            lattice=None,
            points_per_atom=args.points_per_atom,
            bin_ratio=ESP_BIN_RATIO,
            seed=args.seed,
            min_atom_distance=args.min_distance,
            far_span=args.far_span,
            avoid_charges=background,
        )
        sampler_dict = sampler.as_dict()
        structure = {
            "natoms": len(source_positions),
            "symbols": symbols,
            "formula": formula(symbols),
            "positions_Angstrom": source_positions.tolist(),
            "xyz": str(Path(args.xyz).expanduser()),
        }
    else:
        points, reused_sampler = read_points(args.points)
        sampler_dict = {
            "kind": "reused grid",
            "npoints": int(len(points)),
            "points_per_atom": args.points_per_atom,
            "original_sampler": reused_sampler,
        }
        structure = None
        source_positions = None

    report: dict = {
        "format": "esp_tools/orca-esp-v1",
        "method": "ORCA orca_vpot on the SCF density",
        "units": {
            "positions": "Angstrom",
            "potential": "eV/e (= V)",
            "raw_orca": "Hartree/e at coordinates in Bohr",
        },
        "structure": structure,
        "sampler": sampler_dict,
        "background_charges_A_e": background.tolist(),
        "electric_field_au": field.tolist(),
        "gauge_origin_Angstrom": [0.0, 0.0, 0.0],
        "npoints": int(len(points)),
    }

    if args.gbw:
        gbw = Path(args.gbw).expanduser().resolve(strict=True)
        work = gbw.parent
        if not gbw.with_suffix(".densities").is_file() and not (work / f"{gbw.stem}.scfp").is_file():
            raise FileNotFoundError(
                f"matching density data not found for {gbw}: expected "
                f"{gbw.with_suffix('.densities')} or {work / (gbw.stem + '.scfp')}"
            )
        backend = OrcaBackend(
            home=args.orca_home,
            container=args.orca_container,
            apptainer=args.apptainer,
        )
        molecular_ha, provenance = evaluate_points(
            backend, gbw, points, work, prefix=f"esp_{gbw.stem}_{time.time_ns()}"
        )
        background_ev, field_ev = external_components(points, background, field)
        molecular_ev = molecular_ha * HARTREE_EV
        total_ev = molecular_ev + background_ev + field_ev
        report["esp"] = {
            "positions_Angstrom": points.tolist(),
            "phi_molecular_eV_per_e": molecular_ev.tolist(),
            "phi_background_eV_per_e": background_ev.tolist(),
            "phi_uniform_field_eV_per_e": field_ev.tolist(),
            "phi_eV_per_e": total_ev.tolist(),
        }
        report["provenance"] = provenance
        report["definition"] = (
            "phi_total = phi_molecular (nuclei + SCF electrons) "
            "+ phi_background + phi_uniform_field"
        )
    else:
        report["esp"] = {"positions_Angstrom": points.tolist()}
        report["definition"] = "grid only; no SCF density was evaluated"

    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.force else "x"
    with output.open(mode) as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")

    summary = f"Wrote {len(points)} points to {output}"
    if args.xyz:
        summary += (
            f"\n  bin counts (requested/actual): {list(sampler.requested_bin_counts)}"
            f" / {list(sampler.actual_bin_counts)}"
        )
    if args.gbw:
        phi = report["esp"]["phi_eV_per_e"]
        summary += f"\n  phi range: {min(phi):.6f} .. {max(phi):.6f} eV/e"
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
