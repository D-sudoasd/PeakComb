"""Axis conversion, windowing, and component analysis."""

from __future__ import annotations

import math

import numpy as np

from .kernel import pseudo_voigt, pseudo_voigt_area
from .models import (
    AXIS_D_A,
    AXIS_Q_INV_A,
    AXIS_TWO_THETA_DEG,
    Component,
    PeakSeed,
    Session,
    Spectrum,
)


def two_theta_deg_from_d(d_a: float, wavelength_a: float) -> float | None:
    arg = float(wavelength_a) / (2.0 * float(d_a))
    if not math.isfinite(arg) or abs(arg) > 1.0 or d_a == 0.0:
        return None
    return math.degrees(2.0 * math.asin(arg))


def d_from_two_theta_deg(two_theta_deg: float, wavelength_a: float) -> float | None:
    theta = math.radians(float(two_theta_deg) / 2.0)
    s = math.sin(theta)
    if s == 0.0 or not math.isfinite(s):
        return None
    return float(wavelength_a) / (2.0 * s)


def q_from_d(d_a: float) -> float | None:
    if d_a == 0.0 or not math.isfinite(d_a):
        return None
    return 2.0 * math.pi / float(d_a)


def d_from_q(q_inv_a: float) -> float | None:
    if q_inv_a == 0.0 or not math.isfinite(q_inv_a):
        return None
    return 2.0 * math.pi / float(q_inv_a)


def center_to_d(center: float, axis: str, wavelength_a: float) -> float | None:
    if axis == AXIS_D_A:
        return float(center)
    if axis == AXIS_Q_INV_A:
        return d_from_q(center)
    if axis == AXIS_TWO_THETA_DEG:
        return d_from_two_theta_deg(center, wavelength_a)
    raise ValueError(f"Unknown axis: {axis!r}.")


def d_to_center(d_a: float, axis: str, wavelength_a: float) -> float | None:
    if axis == AXIS_D_A:
        return float(d_a)
    if axis == AXIS_Q_INV_A:
        return q_from_d(d_a)
    if axis == AXIS_TWO_THETA_DEG:
        return two_theta_deg_from_d(d_a, wavelength_a)
    raise ValueError(f"Unknown axis: {axis!r}.")


def window_mask(spectrum: Spectrum, xmin: float, xmax: float) -> np.ndarray:
    return (spectrum.x >= xmin) & (spectrum.x <= xmax)


def slice_window(spectrum: Spectrum, xmin: float, xmax: float) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    mask = window_mask(spectrum, xmin, xmax)
    sigma = None if spectrum.sigma is None else spectrum.sigma[mask]
    return spectrum.x[mask], spectrum.y[mask], sigma


def half_max_fwhm(x: np.ndarray, y: np.ndarray) -> float:
    """Estimate FWHM of the strongest peak in a window from half-maximum width."""
    if x.size < 5:
        raise ValueError("Need at least 5 points to estimate FWHM.")
    baseline = float(np.nanmedian([y[0], y[-1]]))
    net = np.asarray(y, dtype=float) - baseline
    peak = float(np.nanmax(net))
    if not math.isfinite(peak) or peak <= 0.0:
        raise ValueError("No positive peak above the endpoint baseline.")
    above = net >= 0.5 * peak
    if not np.any(above):
        raise ValueError("Cannot find half-maximum crossing.")
    idx = np.flatnonzero(above)
    return float(x[idx[-1]] - x[idx[0]])


def evaluate_model(
    x: np.ndarray,
    seeds: list[PeakSeed],
    heights: np.ndarray,
    *,
    hwhm: float,
    shape: float,
    background: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    peak_curves = np.zeros((len(seeds), x.size), dtype=float)
    for index, (seed, height) in enumerate(zip(seeds, heights, strict=True)):
        peak_curves[index] = pseudo_voigt(
            x,
            height=float(height),
            center=float(seed.center),
            hwhm=hwhm,
            shape=shape,
        )
    y_peaks = peak_curves.sum(axis=0) if len(seeds) else np.zeros_like(x, dtype=float)
    y_bkg = np.zeros_like(x, dtype=float) if background is None else np.asarray(background, dtype=float)
    return y_bkg, y_peaks + y_bkg, peak_curves


def components_from_arrays(
    session: Session,
    seeds: list[PeakSeed],
    heights: np.ndarray,
    *,
    fwhm: float | None = None,
    shape: float | None = None,
) -> list[Component]:
    used_fwhm = float(session.fwhm if fwhm is None else fwhm)
    used_shape = float(session.shape if shape is None else shape)
    hwhm = 0.5 * used_fwhm
    areas = np.array(
        [pseudo_voigt_area(float(height), hwhm, used_shape) for height in heights],
        dtype=float,
    )
    total = float(np.sum(np.clip(areas, 0.0, None)))
    components: list[Component] = []
    for index, (seed, height, area) in enumerate(zip(seeds, heights, areas, strict=True)):
        d_a = center_to_d(seed.center, session.axis, session.wavelength_a)
        strain = None
        if d_a is not None and session.d0 not in (None, 0.0):
            strain = (d_a - float(session.d0)) / float(session.d0)
        two_theta = None if d_a is None else two_theta_deg_from_d(d_a, session.wavelength_a)
        q_val = None if d_a is None else q_from_d(d_a)
        area_frac = 0.0 if total <= 0.0 else float(max(area, 0.0) / total)
        label = seed.label or f"P{index + 1}"
        components.append(
            Component(
                index=index,
                label=label,
                center=float(seed.center),
                height=float(height),
                area=float(area),
                area_frac=area_frac,
                fwhm=used_fwhm,
                hwhm=hwhm,
                shape=used_shape,
                d_a=d_a,
                q_inv_a=q_val,
                two_theta_deg=two_theta,
                strain=strain,
                visible=seed.visible,
            )
        )
    return components


def r_metrics(y_obs: np.ndarray, y_fit: np.ndarray, sigma: np.ndarray | None = None) -> dict[str, float]:
    residual = np.asarray(y_obs, dtype=float) - np.asarray(y_fit, dtype=float)
    denom = np.sum(np.square(y_obs))
    rwp = math.sqrt(np.sum(np.square(residual)) / denom) if denom > 0.0 else math.nan
    abs_obs = np.sum(np.abs(y_obs))
    r_factor = float(np.sum(np.abs(residual)) / abs_obs) if abs_obs > 0.0 else math.nan
    wssr = None
    if sigma is not None and np.all(sigma > 0.0):
        wssr = float(np.sum(np.square(residual / sigma)))
    return {
        "rwp": float(rwp),
        "r_factor": float(r_factor),
        "rms_residual": float(np.sqrt(np.mean(np.square(residual)))),
        "wssr": wssr if wssr is not None else float(np.sum(np.square(residual))),
    }
