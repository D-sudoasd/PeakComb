"""Peak-center seeding for Domain and Distribution modes."""

from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks

from .limits import MAX_COMB_PEAKS
from .models import FitMode, PeakSeed, SeedMethod, Session, Spectrum
from .analysis import slice_window


def equal_centers(xmin: float, xmax: float, n_peaks: int) -> np.ndarray:
    if n_peaks < 1:
        raise ValueError("n_peaks must be >= 1.")
    if n_peaks == 1:
        return np.array([(xmin + xmax) / 2.0], dtype=float)
    span = xmax - xmin
    inset = span / (2.0 * n_peaks)
    return np.linspace(xmin + inset, xmax - inset, n_peaks)


def grid_centers(xmin: float, xmax: float, fwhm: float, spacing_factor: float) -> np.ndarray:
    step = float(spacing_factor) * float(fwhm)
    if step <= 0.0:
        raise ValueError("spacing_factor * fwhm must be positive.")
    inset = 0.5 * float(fwhm)
    left = xmin + inset
    right = xmax - inset
    if right <= left:
        return np.array([(xmin + xmax) / 2.0], dtype=float)
    n = int(np.floor((right - left) / step)) + 1
    n = min(max(n, 2), MAX_COMB_PEAKS)
    return np.linspace(left, right, n)


def symmetric_centers(center: float, n_peaks: int, step: float) -> np.ndarray:
    """Odd or even comb mirrored about `center`: ..., c-2Δ, c-Δ, c, c+Δ, c+2Δ."""
    n = int(n_peaks)
    if n < 1:
        raise ValueError("n_peaks must be >= 1.")
    if step <= 0.0:
        raise ValueError("symmetric comb step must be positive.")
    if n == 1:
        return np.array([float(center)], dtype=float)
    offsets = (np.arange(n, dtype=float) - 0.5 * (n - 1)) * float(step)
    return float(center) + offsets


def peak_center(x: np.ndarray, y: np.ndarray) -> float:
    return float(x[int(np.argmax(y))])


def maxima_centers(x: np.ndarray, y: np.ndarray, n_peaks: int, fwhm: float) -> np.ndarray:
    distance = max(int(fwhm / max(np.median(np.diff(x)), 1e-12)), 1)
    prominence = 0.05 * (float(np.max(y)) - float(np.min(y)))
    peaks, _ = find_peaks(y, distance=distance, prominence=max(prominence, 0.0))
    if peaks.size == 0:
        return equal_centers(float(x.min()), float(x.max()), n_peaks)
    order = np.argsort(y[peaks])[::-1]
    chosen = np.sort(peaks[order[:n_peaks]])
    centers = x[chosen]
    if centers.size < n_peaks:
        extra = equal_centers(float(x.min()), float(x.max()), n_peaks)
        merged = np.unique(np.concatenate([centers, extra]))
        return merged[:n_peaks]
    return centers


def seed_peaks(session: Session, spectrum: Spectrum) -> list[PeakSeed]:
    xmin, xmax = session.window(spectrum)
    comb_xmin, comb_xmax = session.comb_window(spectrum)
    x, y, _ = slice_window(spectrum, xmin, xmax)
    if x.size < 3:
        raise ValueError("Fitting window contains too few points.")
    if session.mode is FitMode.DISTRIBUTION:
        if session.seed_method is SeedMethod.GRID:
            centers = grid_centers(comb_xmin, comb_xmax, session.fwhm, session.spacing_factor)
        elif session.seed_method is SeedMethod.SYMMETRIC:
            n_peaks = min(max(int(session.n_peaks), 1), MAX_COMB_PEAKS)
            step = float(session.spacing_factor) * float(session.fwhm)
            center = peak_center(x, y)
            centers = symmetric_centers(center, n_peaks, step)
            centers = centers[(centers >= comb_xmin) & (centers <= comb_xmax)]
            if centers.size == 0:
                centers = np.array([center], dtype=float)
        else:
            n_peaks = min(max(int(session.n_peaks), 2), MAX_COMB_PEAKS)
            centers = equal_centers(comb_xmin, comb_xmax, n_peaks)
    elif session.seed_method is SeedMethod.GRID:
        centers = grid_centers(comb_xmin, comb_xmax, session.fwhm, session.spacing_factor)
    elif session.seed_method is SeedMethod.MANUAL and session.peaks:
        centers = np.array(
            [peak.center for peak in session.peaks if xmin <= peak.center <= xmax],
            dtype=float,
        )
        if centers.size == 0:
            centers = equal_centers(comb_xmin, comb_xmax, session.n_peaks)
    elif session.seed_method is SeedMethod.MAXIMA:
        centers = maxima_centers(x, y, session.n_peaks, session.fwhm)
    else:
        centers = equal_centers(comb_xmin, comb_xmax, session.n_peaks)
    y_span = float(np.max(y) - np.min(y)) if y.size else 1.0
    seeds: list[PeakSeed] = []
    for index, center in enumerate(centers):
        nearest = int(np.argmin(np.abs(x - center)))
        height = max(float(y[nearest] - np.min(y)), 0.05 * y_span)
        seeds.append(PeakSeed(center=float(center), height=height, label=f"P{index + 1}"))
    session.n_peaks = len(seeds)
    session.peaks = seeds
    return seeds
