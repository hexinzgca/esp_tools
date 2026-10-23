"""Sample actual LOCPOT grid points instead of interpolating.

The four distance bins and their ``1:4:4:1`` target are identical to the
continuous sampler.  The difference is that every selected point is an exact
``LOCPOT`` grid point ``(i, j, k)``; its potential is read directly, so there
is no interpolation error.

The FFT/LOCPOT grid is finite: if a bin cannot supply its quota (far field in a
small/dense cell), the VASP policy reallocates the deficit to ``[1.5, 5)``
first, then lower bins, and records the warning metadata.
"""

from __future__ import annotations

import numpy as np

from common.sampling import (
    ESP_BIN_EDGES,
    ESP_BIN_RATIO,
    ESP_GRID_SEED,
    SamplingInfo,
    bin_indices,
    bin_quotas,
    nearest_atom_distances,
    structure_seed,
)


def grid_indices_to_points(
    indices: np.ndarray,
    lattice: np.ndarray,
    grid_shape: tuple[int, int, int],
) -> np.ndarray:
    """Convert C-order LOCPOT indices to Cartesian Angstrom coordinates."""
    nx, ny, nz = grid_shape
    indices = np.asarray(indices, dtype=np.int64)
    i = indices // (ny * nz)
    j = (indices // nz) % ny
    k = indices % nz
    frac = np.column_stack((i / nx, j / ny, k / nz))
    return frac @ np.asarray(lattice, dtype=np.float64)


def sample_lopot_grid(
    positions: np.ndarray,
    lattice: np.ndarray,
    grid_shape: tuple[int, int, int],
    *,
    points_per_atom: int = 10,
    bin_ratio=ESP_BIN_RATIO,
    seed: int | None = ESP_GRID_SEED,
    min_atom_distance: float = ESP_BIN_EDGES[0],
    shortfall_policy: str = "reallocate",
    chunk: int = 1_000_000,
) -> tuple[np.ndarray, np.ndarray, SamplingInfo]:
    """Return ``(flat_indices, points_Angstrom, info)`` from exact grid points.

    ``flat_indices`` index ``potential.ravel(order="C")``.  The caller reads
    ``potential.ravel(order="C")[flat_indices]`` directly.
    """
    positions = np.asarray(positions, dtype=np.float64)
    lattice = np.asarray(lattice, dtype=np.float64)
    nx, ny, nz = (int(x) for x in grid_shape)
    nxyz = nx * ny * nz
    natoms = len(positions)
    total = int(points_per_atom) * natoms
    if total <= 0:
        raise ValueError("points_per_atom * natoms must be positive")
    if shortfall_policy not in {"error", "reallocate"}:
        raise ValueError(f"unknown shortfall_policy: {shortfall_policy!r}")

    want = bin_quotas(total, bin_ratio)
    pools: list[list[np.ndarray]] = [[], [], [], []]

    # Classify every C-order grid point by its nearest-atom distance.
    for start in range(0, nxyz, chunk):
        stop = min(start + chunk, nxyz)
        flat = np.arange(start, stop, dtype=np.int64)
        i = flat // (ny * nz)
        j = (flat // nz) % ny
        k = flat % nz
        frac = np.column_stack((i / nx, j / ny, k / nz))
        cart = frac @ lattice
        distances = nearest_atom_distances(cart, positions, lattice)
        bins = bin_indices(distances)
        bins[distances < float(min_atom_distance)] = -1
        for b in range(4):
            selected = flat[bins == b]
            if len(selected):
                pools[b].append(selected)

    pool_arrays = [
        np.concatenate(pool, axis=0) if pool else np.empty(0, dtype=np.int64)
        for pool in pools
    ]
    rng = np.random.default_rng(structure_seed(positions, seed))
    chosen: list[np.ndarray] = []
    for b in range(4):
        rng.shuffle(pool_arrays[b])
        take = min(want[b], len(pool_arrays[b]))
        chosen.append(pool_arrays[b][:take])
        pool_arrays[b] = pool_arrays[b][take:]

    shortfall = total - sum(len(part) for part in chosen)
    shortfall_bins = tuple(b for b in range(4) if len(chosen[b]) < want[b])
    reallocated = 0
    reallocated_bin_counts = [0, 0, 0, 0]

    if shortfall:
        if shortfall_policy == "error":
            raise RuntimeError(
                f"LOCPOT grid cannot fill requested bins {list(shortfall_bins)}: "
                f"{total - shortfall} / {total} points; requested {want}, "
                f"actual {[len(part) for part in chosen]}. "
                "Use shortfall_policy='reallocate' or a finer FFT grid."
            )
        target_bins = [2, 1, 0] if 3 in shortfall_bins else [
            b for b in (2, 1, 0, 3) if b not in shortfall_bins
        ]
        for target in target_bins:
            if reallocated >= shortfall:
                break
            need = shortfall - reallocated
            take = min(need, len(pool_arrays[target]))
            if take:
                chosen[target] = np.concatenate([chosen[target], pool_arrays[target][:take]])
                pool_arrays[target] = pool_arrays[target][take:]
                reallocated += take
                reallocated_bin_counts[target] += take
        if reallocated < shortfall:
            raise RuntimeError(
                "LOCPOT grid cannot place all requested points even after "
                f"reallocation: {total - shortfall + reallocated} / {total}"
            )

    indices = np.concatenate(chosen, axis=0)
    rng.shuffle(indices)
    points = grid_indices_to_points(indices, lattice, grid_shape)

    dmin = nearest_atom_distances(points, positions, lattice)
    actual = np.bincount(bin_indices(dmin), minlength=4).tolist()
    if sum(actual) != total:
        raise RuntimeError(
            f"internal grid-sampling error: {sum(actual)} / {total} points"
        )
    info = SamplingInfo(
        natoms=natoms,
        npoints=total,
        points_per_atom=int(points_per_atom),
        periodic=True,
        bin_edges_Angstrom=ESP_BIN_EDGES,
        bin_ratio=tuple(float(x) for x in bin_ratio),
        requested_bin_counts=tuple(want),
        actual_bin_counts=tuple(actual),
        base_seed=seed,
        structure_seed=structure_seed(positions, seed),
        min_atom_distance_Angstrom=float(min_atom_distance),
        far_span_Angstrom=0.0,
        charge_min_distance_Angstrom=0.0,
        shortfall_policy=shortfall_policy,
        shortfall_bins=shortfall_bins,
        reallocated_points=reallocated,
        reallocated_bin_counts=tuple(reallocated_bin_counts),
        sampling_mode="dft-grid",
    )
    return indices, points, info
