"""Minimal VASP POSCAR/LOCPOT reader and periodic trilinear interpolation.

Only the data written by ``LVHAR = .TRUE.`` is interpreted: the local
electrostatic potential ``V_ionic + V_Hartree`` in eV.  The file layout is
the POSCAR header followed by one or more ``NGXF NGYF NGZF`` blocks.  VASP
stores the values with x varying fastest; the array exposed here has shape
``(NGXF, NGYF, NGZF)`` and is indexed as ``potential[i, j, k]``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

import numpy as np


def _is_int_token(token: str) -> bool:
    try:
        int(token)
    except ValueError:
        return False
    return True


def _scale_from_line(tokens: list[str]) -> float:
    if len(tokens) != 1:
        raise ValueError(
            "only a single positive POSCAR scale factor is supported; "
            f"found {len(tokens)} values"
        )
    scale = float(tokens[0])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"POSCAR scale must be positive and finite, got {scale}")
    return scale


@dataclass(frozen=True)
class PoscarGeometry:
    lattice: np.ndarray
    positions: np.ndarray
    symbols: list[str]
    counts: list[int]


@dataclass(frozen=True)
class Locpot:
    lattice: np.ndarray
    positions: np.ndarray
    symbols: list[str]
    counts: list[int]
    grid_shape: tuple[int, int, int]
    potential: np.ndarray
    n_blocks: int
    spin_block_max_abs_diff_eV: float
    source: str


def _parse_geometry(lines: list[str], start: int = 0) -> tuple[PoscarGeometry, int]:
    """Parse the POSCAR-like header; return geometry and the next line index."""
    if len(lines) < start + 8:
        raise ValueError("file is too short to contain a POSCAR header")
    cursor = start + 1  # skip the comment line
    scale = _scale_from_line(lines[cursor].split())
    cursor += 1

    lattice = []
    for _ in range(3):
        tokens = lines[cursor].split()
        if len(tokens) < 3:
            raise ValueError(f"invalid lattice line: {lines[cursor]!r}")
        lattice.append([float(x) for x in tokens[:3]])
        cursor += 1
    lattice_array = np.asarray(lattice, dtype=np.float64) * scale

    tokens = lines[cursor].split()
    cursor += 1
    if tokens and all(_is_int_token(token) for token in tokens):
        symbols: list[str] = []
        counts = [int(x) for x in tokens]
    else:
        symbols = tokens
        if cursor >= len(lines):
            raise ValueError("missing atom-count line")
        counts = [int(x) for x in lines[cursor].split()]
        cursor += 1

    if not counts or any(count <= 0 for count in counts):
        raise ValueError(f"invalid atom counts: {counts}")
    if symbols and len(symbols) != len(counts):
        raise ValueError("species and atom-count lines have different lengths")

    if cursor < len(lines) and lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines):
        raise ValueError("missing coordinate-mode line")
    mode = lines[cursor].strip().lower()
    cursor += 1
    if not mode:
        raise ValueError("empty coordinate-mode line")
    direct = mode[0] == "d"
    cartesian = mode[0] == "c" or mode[0] == "k"
    if not (direct or cartesian):
        raise ValueError(f"unknown coordinate mode: {lines[cursor - 1]!r}")

    natoms = int(sum(counts))
    coordinates = []
    for _ in range(natoms):
        if cursor >= len(lines):
            raise ValueError("truncated coordinate block")
        tokens = lines[cursor].split()
        if len(tokens) < 3:
            raise ValueError(f"invalid coordinate line: {lines[cursor]!r}")
        coordinates.append([float(x) for x in tokens[:3]])
        cursor += 1
    coordinates_array = np.asarray(coordinates, dtype=np.float64)
    if not np.isfinite(coordinates_array).all() or not np.isfinite(lattice_array).all():
        raise ValueError("non-finite geometry")

    if direct:
        positions = coordinates_array @ lattice_array
    else:
        positions = coordinates_array * scale

    return PoscarGeometry(
        lattice=lattice_array,
        positions=positions,
        symbols=symbols,
        counts=counts,
    ), cursor


def read_poscar(path: str | Path) -> PoscarGeometry:
    lines = Path(path).read_text().splitlines()
    geometry, _ = _parse_geometry(lines)
    return geometry


def read_locpot(path: str | Path) -> Locpot:
    path = Path(path)
    lines = path.read_text().splitlines()
    geometry, cursor = _parse_geometry(lines)

    # Skip blank lines up to the first grid-dimension line.
    while cursor < len(lines) and not lines[cursor].strip():
        cursor += 1
    if cursor >= len(lines):
        raise ValueError("LOCPOT has no grid-dimension line")

    nblocks = 0
    raw_blocks: list[np.ndarray] = []
    grid_shape: tuple[int, int, int] | None = None

    while cursor < len(lines):
        line = lines[cursor].strip()
        cursor += 1
        if not line:
            continue
        tokens = line.split()
        if len(tokens) != 3 or not all(_is_int_token(token) for token in tokens):
            if raw_blocks:
                break  # tolerate trailing annotations
            raise ValueError(f"expected 'NGXF NGYF NGZF', got {line!r}")
        shape = tuple(int(x) for x in tokens)
        if any(n <= 0 for n in shape):
            raise ValueError(f"invalid grid shape: {shape}")
        if grid_shape is None:
            grid_shape = shape
        elif shape != grid_shape:
            raise ValueError(f"inconsistent grid shape {shape} vs {grid_shape}")
        nxyz = int(np.prod(shape))

        values: list[float] = []
        while len(values) < nxyz and cursor < len(lines):
            data_line = lines[cursor].strip()
            cursor += 1
            if not data_line:
                continue
            try:
                values.extend(float(x) for x in data_line.split())
            except ValueError as exc:
                raise ValueError(
                    f"invalid potential value near line {cursor}: {data_line!r}"
                ) from exc
        if len(values) != nxyz:
            raise ValueError(f"grid block has {len(values)} values, expected {nxyz}")

        # VASP writes x fastest, then y, then z.
        array = np.asarray(values, dtype=np.float64).reshape(shape[::-1]).transpose(2, 1, 0)
        raw_blocks.append(array)
        nblocks += 1

    if grid_shape is None or not raw_blocks:
        raise ValueError("LOCPOT contains no potential grid")

    stack = np.stack(raw_blocks, axis=0)
    potential = stack.mean(axis=0)
    max_diff = float(np.max(np.abs(stack - potential))) if nblocks > 1 else 0.0
    if not np.isfinite(potential).all():
        raise ValueError("LOCPOT contains non-finite potential values")

    return Locpot(
        lattice=geometry.lattice,
        positions=geometry.positions,
        symbols=geometry.symbols,
        counts=geometry.counts,
        grid_shape=tuple(int(x) for x in grid_shape),
        potential=potential,
        n_blocks=nblocks,
        spin_block_max_abs_diff_eV=max_diff,
        source=str(path),
    )


def interpolate_periodic(
    points_angstrom: np.ndarray,
    lattice: np.ndarray,
    potential: np.ndarray,
) -> np.ndarray:
    """Trilinear interpolation on a periodic VASP grid.

    ``potential`` has shape ``(NGXF, NGYF, NGZF)``.  Points outside the cell
    are folded by periodicity, which is the correct interpretation of a
    periodic LOCPOT.
    """
    points = np.asarray(points_angstrom, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    shape = np.asarray(potential.shape, dtype=np.int64)
    inv_lattice = np.linalg.inv(np.asarray(lattice, dtype=np.float64))
    frac = points @ inv_lattice
    u = frac * shape
    base = np.floor(u).astype(np.int64)
    weight = u - base
    i0 = np.mod(base, shape)
    i1 = np.mod(base + 1, shape)

    wx, wy, wz = weight[:, 0], weight[:, 1], weight[:, 2]
    i0x, i0y, i0z = i0[:, 0], i0[:, 1], i0[:, 2]
    i1x, i1y, i1z = i1[:, 0], i1[:, 1], i1[:, 2]

    c000 = potential[i0x, i0y, i0z]
    c100 = potential[i1x, i0y, i0z]
    c010 = potential[i0x, i1y, i0z]
    c110 = potential[i1x, i1y, i0z]
    c001 = potential[i0x, i0y, i1z]
    c101 = potential[i1x, i0y, i1z]
    c011 = potential[i0x, i1y, i1z]
    c111 = potential[i1x, i1y, i1z]

    c00 = c000 * (1 - wx) + c100 * wx
    c10 = c010 * (1 - wx) + c110 * wx
    c01 = c001 * (1 - wx) + c101 * wx
    c11 = c011 * (1 - wx) + c111 * wx
    c0 = c00 * (1 - wy) + c10 * wy
    c1 = c01 * (1 - wy) + c11 * wy
    return c0 * (1 - wz) + c1 * wz


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
