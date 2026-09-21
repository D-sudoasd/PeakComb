"""Fityk-compatible height-parameterized PseudoVoigt kernel."""

from __future__ import annotations

import math

import numpy as np

LN2 = math.log(2.0)
GAUSS_AREA_FACTOR = math.sqrt(math.pi / LN2)
LORENTZ_AREA_FACTOR = math.pi


def hwhm_from_fwhm(fwhm: float) -> float:
    return 0.5 * float(fwhm)


def fwhm_from_hwhm(hwhm: float) -> float:
    return 2.0 * float(hwhm)


def unit_area_factor(shape: float) -> float:
    """Integral of a unit-height, unit-hwhm PseudoVoigt."""
    eta = float(shape)
    return (1.0 - eta) * GAUSS_AREA_FACTOR + eta * LORENTZ_AREA_FACTOR


def pseudo_voigt_area(height: float, hwhm: float, shape: float) -> float:
    return float(height) * abs(float(hwhm)) * unit_area_factor(shape)


def pseudo_voigt(
    x: np.ndarray,
    *,
    height: float,
    center: float,
    hwhm: float,
    shape: float,
) -> np.ndarray:
    """Fityk PseudoVoigt(height, center, hwhm, shape).

    shape=0 is Gaussian, shape=1 is Lorentzian. Both components share hwhm.
    """
    width = abs(float(hwhm))
    if width == 0.0:
        raise ValueError("hwhm must be non-zero.")
    eta = float(shape)
    if eta < 0.0 or eta > 1.0:
        raise ValueError("shape must be in [0, 1].")
    z = (np.asarray(x, dtype=float) - float(center)) / width
    gaussian = np.exp(-LN2 * z * z)
    lorentzian = 1.0 / (1.0 + z * z)
    return float(height) * ((1.0 - eta) * gaussian + eta * lorentzian)


def background_columns(
    x: np.ndarray,
    kind: str,
    *,
    x_mid: float | None = None,
) -> np.ndarray | None:
    """Return design columns for a local polynomial background, or None."""
    if kind in {"none", None}:
        return None
    xv = np.asarray(x, dtype=float)
    mid = float(xv.mean()) if x_mid is None else float(x_mid)
    dx = xv - mid
    if kind == "constant":
        return np.ones((xv.size, 1), dtype=float)
    if kind == "linear":
        return np.column_stack([np.ones(xv.size), dx])
    if kind == "quadratic":
        return np.column_stack([np.ones(xv.size), dx, dx * dx])
    raise ValueError(f"Unknown background kind: {kind!r}.")


def peak_design_matrix(
    x: np.ndarray,
    centers: np.ndarray,
    *,
    hwhm: float,
    shape: float,
) -> np.ndarray:
    centers = np.asarray(centers, dtype=float).reshape(-1)
    if centers.size == 0:
        raise ValueError("At least one peak center is required.")
    width = abs(float(hwhm))
    if width == 0.0:
        raise ValueError("hwhm must be non-zero.")
    eta = float(shape)
    z = (np.asarray(x, dtype=float)[:, None] - centers[None, :]) / width
    gaussian = np.exp(-LN2 * z * z)
    lorentzian = 1.0 / (1.0 + z * z)
    return (1.0 - eta) * gaussian + eta * lorentzian
