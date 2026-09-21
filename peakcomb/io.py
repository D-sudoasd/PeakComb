"""Load two-column XRD spectra from common ASCII formats."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .models import AXIS_D_A, DEFAULT_ENERGY_KEV, DEFAULT_WAVELENGTH_A, Spectrum, X_RAY_KEV_ANGSTROM

TEXT_ENCODINGS = ("utf-8-sig", "utf-8", "gb18030", "cp936", "latin1")
SUPPORTED_SUFFIXES = {".xy", ".xye", ".chi", ".csv", ".txt", ".dat"}


class SpectrumLoadError(ValueError):
    """Raised when a spectrum file cannot be parsed."""


def wavelength_from_energy_kev(energy_kev: float) -> float:
    return X_RAY_KEV_ANGSTROM / float(energy_kev)


def energy_from_wavelength_a(wavelength_a: float) -> float:
    return X_RAY_KEV_ANGSTROM / float(wavelength_a)


def _decode(path: Path) -> str:
    data = path.read_bytes()
    last_error: Exception | None = None
    for encoding in TEXT_ENCODINGS:
        try:
            return data.decode(encoding)
        except UnicodeDecodeError as exc:
            last_error = exc
    raise SpectrumLoadError(f"Cannot decode {path}: {last_error}")


def _parse_row(line: str) -> list[float] | None:
    stripped = line.strip()
    if not stripped or stripped[0] in "#@!;/%":
        return None
    if "," in stripped and "\t" not in stripped:
        parts = stripped.replace(";", " ").replace(",", " ").split()
    else:
        parts = stripped.replace(",", " ").split()
    values: list[float] = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            if values:
                break
            return None
    if len(values) < 2:
        return None
    return values


def load_spectrum(
    path: str | Path,
    *,
    axis: str = AXIS_D_A,
    wavelength_a: float | None = None,
    energy_kev: float | None = None,
) -> Spectrum:
    file_path = Path(path)
    if not file_path.is_file():
        raise SpectrumLoadError(f"File not found: {file_path}")
    text = _decode(file_path)
    rows: list[list[float]] = []
    for line in text.splitlines():
        parsed = _parse_row(line)
        if parsed is not None:
            rows.append(parsed)
    if len(rows) >= 3 and len(rows[0]) == 1 and len(rows[1]) == 1:
        # GSAS .chi style: column-count, npoints, then data.
        rows = [row for row in rows[2:] if len(row) >= 2]
    if len(rows) < 3:
        raise SpectrumLoadError(f"{file_path.name}: need at least 3 numeric XY rows.")
    array = np.asarray(rows, dtype=float)
    x = array[:, 0]
    y = array[:, 1]
    sigma = array[:, 2] if array.shape[1] >= 3 else None
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    if sigma is not None:
        sigma = sigma[order]
        if not np.all(np.isfinite(sigma)) or np.any(sigma <= 0.0):
            sigma = None
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if sigma is not None:
        sigma = sigma[finite]
    if x.size < 3:
        raise SpectrumLoadError(f"{file_path.name}: too few finite points.")
    resolved_wavelength = DEFAULT_WAVELENGTH_A
    resolved_energy = energy_kev
    if wavelength_a is not None:
        resolved_wavelength = float(wavelength_a)
        resolved_energy = energy_from_wavelength_a(resolved_wavelength)
    elif energy_kev is not None:
        resolved_energy = float(energy_kev)
        resolved_wavelength = wavelength_from_energy_kev(resolved_energy)
    elif axis != AXIS_D_A:
        resolved_energy = DEFAULT_ENERGY_KEV
        resolved_wavelength = wavelength_from_energy_kev(resolved_energy)
    return Spectrum(
        x=x,
        y=y,
        sigma=sigma,
        path=file_path,
        axis=axis,
        wavelength_a=resolved_wavelength,
        energy_kev=resolved_energy,
    )


def write_xy(path: str | Path, x: np.ndarray, y: np.ndarray) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# PeakComb XY\n")
        for xv, yv in zip(np.asarray(x, dtype=float), np.asarray(y, dtype=float), strict=True):
            handle.write(f"{xv:.12g}\t{yv:.12g}\n")
    return out
