"""ESP grid sampling shared by the VASP and ORCA extractors.

The sampler produces exactly ``points_per_atom * natoms`` probe points.
Points are distributed over four bins defined by the distance to the
nearest atom:

    [0.25, 0.5) : [0.5, 1.5) : [1.5, 5.0) : [5.0, inf)  =  1 : 4 : 4 : 1

For a neutral molecule these weights put only 10 % of the points in the
innermost shell, 80 % in the chemically relevant 0.5-5 A region, and 10 %
in the far field.  The innermost shell starts at 0.25 A, so no point sits
on a nucleus.  This is the sampling scheme of the reference ``esp_grid``
tool, with the total count changed from a per-atom target rounded to a
multiple of 14 to the exact ``10 * natoms`` requested here.

Two modes are supported:

``lattice is None``
    Non-periodic molecule (ORCA).  Near-field candidates are generated
    around randomly chosen atoms; far-field candidates lie on a shell
    around the centroid with radius ``r_max + 5 ... r_max + 10`` A.

``lattice is not None``
    Periodic cell (VASP).  All points are wrapped into the cell and all
    distances use the minimum-image convention.  Far-field candidates are
    uniform random points in the cell accepted only if their minimum-image
    distance to every atom is at least 5 A.  The cell must therefore
    contain enough vacuum; a 16 A cubic box is comfortable for water.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Sequence
import warnings

import numpy as np

BOHR_ANGSTROM = 0.529177210903
HARTREE_EV = 27.211386245988

ESP_BIN_EDGES = (0.25, 0.5, 1.5, 5.0)
ESP_BIN_RATIO = (1, 4, 4, 1)
ESP_GRID_SEED = 20260916
ESP_FAR_SPAN = 5.0
ESP_POINTS_PER_ATOM = 10
ESP_CHARGE_MIN_DISTANCE = 0.30


def structure_seed(positions: np.ndarray, seed: int | None = ESP_GRID_SEED) -> int | None:
    """Mix the base seed with a geometry hash.

    The same structure always gives the same grid; different structures do
    not share the same (atom, shell, radius, orientation) sequence.
    """
    if seed is None:
        return None
    digest = hashlib.sha256(
        np.round(np.asarray(positions, dtype=np.float64), 6).tobytes()
    )
    return (int(seed) + int.from_bytes(digest.digest()[:8], "little")) % (2 ** 63)


def fibonacci_directions(n: int, rng: np.random.Generator, start: int = 0) -> np.ndarray:
    """Approximately uniform sphere directions (golden-angle spiral + random rotation)."""
    idx = np.arange(start, start + n, dtype=np.float64) + 0.5 + rng.random()
    z = 1.0 - 2.0 * idx / (start + n + 1.0)
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    theta = np.pi * (3.0 - np.sqrt(5.0)) * idx
    vecs = np.stack([np.cos(theta) * r, z, np.sin(theta) * r], axis=1)

    quat = rng.normal(size=4)
    quat /= np.linalg.norm(quat)
    w, x, y, zq = quat
    rot = np.array([
        [1 - 2 * (y * y + zq * zq), 2 * (x * y - zq * w), 2 * (x * zq + y * w)],
        [2 * (x * y + zq * w), 1 - 2 * (x * x + zq * zq), 2 * (y * zq - x * w)],
        [2 * (x * zq - y * w), 2 * (y * zq + x * w), 1 - 2 * (x * x + y * y)],
    ])
    return vecs @ rot.T


def wrap_points(points: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """Wrap Cartesian points into the primary cell defined by row-vector lattice."""
    inv = np.linalg.inv(np.asarray(lattice, dtype=np.float64))
    frac = np.asarray(points, dtype=np.float64) @ inv
    frac -= np.floor(frac)
    frac[np.isclose(frac, 1.0, atol=1e-12)] = 0.0
    return frac @ np.asarray(lattice, dtype=np.float64)


def nearest_atom_distances(
    points: np.ndarray,
    positions: np.ndarray,
    lattice: np.ndarray | None = None,
    chunk: int = 8192,
) -> np.ndarray:
    """Distance from each point to the nearest atom (A), optionally minimum-image."""
    points = np.asarray(points, dtype=np.float64)
    positions = np.asarray(positions, dtype=np.float64)
    if lattice is not None:
        lattice = np.asarray(lattice, dtype=np.float64)
        if lattice.shape != (3, 3):
            raise ValueError("lattice must have shape (3, 3)")
        inv = np.linalg.inv(lattice)

    out = np.empty(len(points), dtype=np.float64)
    for start in range(0, len(points), chunk):
        block = points[start:start + chunk]
        delta = block[:, None, :] - positions[None, :, :]
        if lattice is not None:
            frac = delta @ inv
            frac -= np.round(frac)
            delta = frac @ lattice
        out[start:start + chunk] = np.sqrt((delta ** 2).sum(axis=-1)).min(axis=1)
    return out


def bin_indices(distances: np.ndarray) -> np.ndarray:
    """0..3 for the four bins; -1 for distances below the first bin."""
    distances = np.asarray(distances, dtype=np.float64)
    idx = np.digitize(distances, ESP_BIN_EDGES) - 1
    idx[distances < ESP_BIN_EDGES[0]] = -1
    return idx


def bin_quotas(total: int, ratio: Sequence[float] = ESP_BIN_RATIO) -> list[int]:
    """Integer quotas that sum exactly to ``total`` (largest remainder method)."""
    total = int(total)
    weights = np.asarray(ratio, dtype=np.float64)
    if weights.size != 4 or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError(f"bin ratio must contain four non-negative weights: {ratio!r}")
    raw = total * weights / weights.sum()
    quotas = np.floor(raw).astype(np.int64)
    remainder = total - int(quotas.sum())
    if remainder:
        order = np.argsort(-(raw - quotas), kind="stable")
        quotas[order[:remainder]] += 1
    if int(quotas.sum()) != total:
        raise RuntimeError("internal error: bin quotas do not sum to the requested total")
    return quotas.tolist()


@dataclass(frozen=True)
class SamplingInfo:
    natoms: int
    npoints: int
    points_per_atom: int
    periodic: bool
    bin_edges_Angstrom: tuple[float, ...]
    bin_ratio: tuple[float, ...]
    requested_bin_counts: tuple[int, ...]
    actual_bin_counts: tuple[int, ...]
    base_seed: int | None
    structure_seed: int | None
    min_atom_distance_Angstrom: float
    far_span_Angstrom: float
    charge_min_distance_Angstrom: float
    shortfall_policy: str
    shortfall_bins: tuple[int, ...]
    reallocated_points: int
    reallocated_bin_counts: tuple[int, ...]

    @property
    def bin_reallocation(self) -> bool:
        return self.reallocated_points > 0

    def as_dict(self) -> dict:
        return {
            "natoms": self.natoms,
            "npoints": self.npoints,
            "points_per_atom": self.points_per_atom,
            "periodic": self.periodic,
            "bin_edges_Angstrom": list(self.bin_edges_Angstrom) + [None],
            "bin_ratio": list(self.bin_ratio),
            "requested_bin_counts": list(self.requested_bin_counts),
            "actual_bin_counts": list(self.actual_bin_counts),
            "base_seed": self.base_seed,
            "structure_seed": self.structure_seed,
            "min_atom_distance_Angstrom": self.min_atom_distance_Angstrom,
            "far_span_Angstrom": self.far_span_Angstrom,
            "charge_min_distance_Angstrom": self.charge_min_distance_Angstrom,
            "shortfall_policy": self.shortfall_policy,
            "shortfall_bins": list(self.shortfall_bins),
            "reallocated_points": self.reallocated_points,
            "reallocated_bin_counts": list(self.reallocated_bin_counts),
            "bin_reallocation": self.bin_reallocation,
        }


def _propose_candidates(
    n_near: int,
    n_far: int,
    *,
    positions: np.ndarray,
    lattice: np.ndarray | None,
    rng: np.random.Generator,
    k: int,
    centroid: np.ndarray,
    r_max: float,
    far_span: float,
    near_bands: Sequence[int] | None = None,
) -> tuple[np.ndarray, int]:
    """Generate one batch of near- and far-field candidate points.

    ``near_bands`` restricts the near-field proposals to selected bins
    (0, 1, 2).  Reallocation uses this to send a far-field shortfall to
    ``[1.5, 5)`` first instead of distributing it randomly.
    """
    if n_near > 0:
        dirs = fibonacci_directions(n_near, rng, k)
        k += n_near
        atom_idx = rng.integers(0, len(positions), n_near)
        if near_bands is None:
            band = rng.integers(0, 3, n_near)
        else:
            band = rng.choice(np.asarray(tuple(near_bands), dtype=np.int64), n_near)
        lo = np.asarray(ESP_BIN_EDGES[:3])[band]
        hi = np.asarray(ESP_BIN_EDGES[1:])[band]
        radius = lo + rng.random(n_near) * (hi - lo)
        near = positions[atom_idx] + dirs * radius[:, None]
        if lattice is not None:
            near = wrap_points(near, lattice)
    else:
        near = np.empty((0, 3))

    if n_far > 0:
        if lattice is not None:
            frac = rng.random((n_far, 3))
            far = frac @ lattice
        else:
            f_dirs = fibonacci_directions(n_far, rng, k)
            k += n_far
            f_radius = r_max + ESP_BIN_EDGES[3] + rng.random(n_far) * far_span
            far = centroid + f_dirs * f_radius[:, None]
    else:
        far = np.empty((0, 3))

    return np.concatenate([near, far], axis=0), k


def sample_grid(
    positions: np.ndarray,
    *,
    lattice: np.ndarray | None = None,
    points_per_atom: int = ESP_POINTS_PER_ATOM,
    bin_ratio: Sequence[float] = ESP_BIN_RATIO,
    seed: int | None = ESP_GRID_SEED,
    min_atom_distance: float = ESP_BIN_EDGES[0],
    far_span: float = ESP_FAR_SPAN,
    avoid_charges: np.ndarray | None = None,
    charge_min_distance: float = ESP_CHARGE_MIN_DISTANCE,
    max_rounds: int = 64,
    shortfall_policy: str = "error",
    warn_on_reallocation: bool = True,
) -> tuple[np.ndarray, SamplingInfo]:
    """Sample a 10 x natoms ESP grid with a 1:4:4:1 nearest-atom-distance split.

    ``shortfall_policy`` controls what happens when a bin cannot be filled:

    ``"error"`` (default)
        Raise.  Used by the ORCA path, where the far shell always exists.
    ``"reallocate"``
        Keep the total ``points_per_atom * natoms`` but move the unfillable
        far-field deficit to ``[1.5, 5)`` first (then to lower bins if that
        also cannot absorb it), and emit a ``RuntimeWarning``.  Used by the
        VASP path, where a periodic cell may have no region >= 5 A from every
        atom (dense solids, small boxes).  Set ``warn_on_reallocation=False``
        when the caller prints its own structured warning.

    Parameters are intentionally explicit so the caller can record them in
    the output provenance.  ``avoid_charges`` has shape (n, 4) with
    ``[x_A, y_A, z_A, q_e]``; it only affects the exclusion of probe points,
    not the potential evaluation.
    """
    if shortfall_policy not in {"error", "reallocate"}:
        raise ValueError(f"unknown shortfall_policy: {shortfall_policy!r}")
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3 or not np.isfinite(positions).all():
        raise ValueError("positions must be a finite (natoms, 3) array in Angstrom")
    natoms = len(positions)
    if natoms == 0:
        raise ValueError("positions must contain at least one atom")
    if int(points_per_atom) < 1:
        raise ValueError("points_per_atom must be positive")
    total = int(points_per_atom) * natoms

    if lattice is not None:
        lattice = np.asarray(lattice, dtype=np.float64)
        if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
            raise ValueError("lattice must be a finite (3, 3) array")
        if abs(np.linalg.det(lattice)) < 1e-12:
            raise ValueError("lattice is singular")

    charges = None
    if avoid_charges is not None and len(np.asarray(avoid_charges)):
        charges = np.asarray(avoid_charges, dtype=np.float64)
        if charges.ndim != 2 or charges.shape[1] != 4 or not np.isfinite(charges).all():
            raise ValueError("avoid_charges must have shape (n, 4) with finite values")
        charges = charges[:, :3]

    want = bin_quotas(total, bin_ratio)
    have = [0, 0, 0, 0]
    pool: list[list[np.ndarray]] = [[], [], [], []]

    rng = np.random.default_rng(structure_seed(positions, seed))
    centroid = positions.mean(axis=0)
    r_max = float(np.linalg.norm(positions - centroid, axis=1).max())
    k = 0

    for _ in range(max_rounds):
        if have == want:
            break
        need = [max(0, want[b] - have[b]) for b in range(4)]

        # Near-field candidates: atom-centred points in the first three bins.
        near_need = need[0] + need[1] + need[2]
        n_near = min(8 * near_need + 512, 200000) if near_need else 0
        # Far-field candidates: shell around the centroid (molecule) or uniform
        # points in the cell (periodic).
        n_far = min(8 * need[3] + 512, 200000) if need[3] else 0

        candidates, k = _propose_candidates(
            n_near,
            n_far,
            positions=positions,
            lattice=lattice,
            rng=rng,
            k=k,
            centroid=centroid,
            r_max=r_max,
            far_span=far_span,
        )
        if not len(candidates):
            break

        dmin = nearest_atom_distances(candidates, positions, lattice)
        keep = dmin >= float(min_atom_distance)
        if charges is not None:
            dchg = np.sqrt(
                ((candidates[:, None, :] - charges[None, :, :]) ** 2).sum(axis=-1)
            ).min(axis=1)
            keep &= dchg >= float(charge_min_distance)
        idx = bin_indices(dmin)

        for b in range(4):
            if need[b] <= 0:
                continue
            selected = np.flatnonzero(keep & (idx == b))
            if not len(selected):
                continue
            selected = selected[:need[b]]
            pool[b].append(candidates[selected])
            have[b] += len(selected)

    shortfall = total - sum(have)
    shortfall_bins = tuple(b for b in range(4) if have[b] < want[b])
    reallocated = 0
    reallocated_bin_counts = [0, 0, 0, 0]

    if shortfall:
        if shortfall_policy == "error":
            raise RuntimeError(
                f"sampler could not fill requested bins {list(shortfall_bins)}: "
                f"{sum(have)} / {total} points; requested {want}, actual {have}. "
                "Pass shortfall_policy='reallocate' to keep the total and "
                "redistribute the deficit."
            )

        # A far-field shortfall goes to [1.5, 5) first, then to the two
        # remaining near bins.  This is the VASP case for a cell without a
        # >=5 A region: keep the requested total, but put the spare points in
        # the outermost chemically relevant shell instead of spreading them
        # randomly.
        if 3 in shortfall_bins:
            target_bins = [2, 1, 0]
        else:
            target_bins = [b for b in (2, 1, 0, 3) if b not in shortfall_bins]

        remaining = shortfall
        for target in target_bins:
            if remaining <= 0:
                break
            for _ in range(max_rounds):
                if remaining <= 0:
                    break
                n_near = min(8 * remaining + 512, 200000) if target in (0, 1, 2) else 0
                n_far = min(8 * remaining + 512, 200000) if target == 3 else 0
                candidates, k = _propose_candidates(
                    n_near,
                    n_far,
                    positions=positions,
                    lattice=lattice,
                    rng=rng,
                    k=k,
                    centroid=centroid,
                    r_max=r_max,
                    far_span=far_span,
                    near_bands=(target,) if target in (0, 1, 2) else None,
                )
                if not len(candidates):
                    break

                dmin = nearest_atom_distances(candidates, positions, lattice)
                keep = dmin >= float(min_atom_distance)
                if charges is not None:
                    dchg = np.sqrt(
                        ((candidates[:, None, :] - charges[None, :, :]) ** 2).sum(axis=-1)
                    ).min(axis=1)
                    keep &= dchg >= float(charge_min_distance)
                idx = bin_indices(dmin)
                selected = np.flatnonzero(keep & (idx == target))
                if not len(selected):
                    break
                selected = selected[:remaining]
                pool[target].append(candidates[selected])
                have[target] += len(selected)
                reallocated += len(selected)
                reallocated_bin_counts[target] += len(selected)
                remaining -= len(selected)

        if remaining > 0:
            raise RuntimeError(
                "sampler could not place all requested points even after "
                f"reallocation: {total - remaining} / {total}; "
                f"requested {want}, actual {have}, reallocated {reallocated}"
            )

    blocks = [
        np.concatenate(pool[b], axis=0) if pool[b] else np.empty((0, 3))
        for b in range(4)
    ]
    points = np.concatenate(blocks, axis=0)
    if len(points) != total:
        raise RuntimeError(
            "sampler could not fill the requested grid: "
            f"{len(points)} / {total} points; actual bin counts {have}, requested {want}"
        )
    if len(points) > 1:
        points = points[rng.permutation(len(points))]
    points = np.ascontiguousarray(points, dtype=np.float64)

    dmin = nearest_atom_distances(points, positions, lattice)
    actual = np.bincount(bin_indices(dmin), minlength=4).tolist()
    if reallocated and warn_on_reallocation:
        warnings.warn(
            f"ESP sampler could not fill requested bins {list(shortfall_bins)}; "
            f"reallocated {reallocated} point(s) to bins {reallocated_bin_counts} "
            f"(far-field deficit goes to [1.5, 5) first). "
            f"requested={want}, actual={actual}",
            RuntimeWarning,
            stacklevel=2,
        )
    info = SamplingInfo(
        natoms=natoms,
        npoints=total,
        points_per_atom=int(points_per_atom),
        periodic=lattice is not None,
        bin_edges_Angstrom=ESP_BIN_EDGES,
        bin_ratio=tuple(float(x) for x in bin_ratio),
        requested_bin_counts=tuple(want),
        actual_bin_counts=tuple(actual),
        base_seed=seed,
        structure_seed=structure_seed(positions, seed),
        min_atom_distance_Angstrom=float(min_atom_distance),
        far_span_Angstrom=float(far_span),
        charge_min_distance_Angstrom=float(charge_min_distance),
        shortfall_policy=shortfall_policy,
        shortfall_bins=shortfall_bins,
        reallocated_points=reallocated,
        reallocated_bin_counts=tuple(reallocated_bin_counts),
    )
    return points, info


def read_xyz(path: str | Path) -> tuple[list[str], np.ndarray]:
    """Read exactly one standard XYZ frame; coordinates are Angstrom."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"invalid XYZ file: {path}")
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"first XYZ line must be an atom count: {path}") from exc
    if natoms <= 0 or len(lines) < natoms + 2:
        raise ValueError(f"XYZ atom count/payload mismatch: {path}")
    symbols: list[str] = []
    positions: list[list[float]] = []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"invalid XYZ coordinate line: {line!r}")
        symbols.append(parts[0])
        positions.append([float(x) for x in parts[1:4]])
    if len(lines) > natoms + 2 and any(line.strip() for line in lines[natoms + 2:]):
        raise ValueError("XYZ must contain exactly one frame")
    positions_array = np.asarray(positions, dtype=np.float64)
    if not np.isfinite(positions_array).all():
        raise ValueError("XYZ coordinates must be finite")
    return symbols, positions_array
