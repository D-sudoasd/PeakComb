"""Fit quality checks for shared-width comb models."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .models import Component, QCReport


def fwhm_relative_span(components: list[Component]) -> float:
    values = [item.fwhm for item in components if math.isfinite(item.fwhm) and item.fwhm > 0.0]
    if not values:
        return math.nan
    mid = 0.5 * (min(values) + max(values))
    if mid == 0.0:
        return math.nan
    return (max(values) - min(values)) / mid


def build_qc(
    components: list[Component],
    metrics: dict[str, Any],
    *,
    fitter: str,
    design_cond: float | None,
    empty_area_frac: float,
    expected_fwhm: float,
    extra_warnings: list[str] | None = None,
) -> QCReport:
    warnings = list(extra_warnings or [])
    n_empty = sum(1 for item in components if item.area_frac < empty_area_frac or item.height <= 0.0)
    span = fwhm_relative_span(components)
    shared_ok = bool(np.isfinite(span) and span < 1e-8)
    if components:
        mean_fwhm = float(np.mean([item.fwhm for item in components]))
        if abs(mean_fwhm - expected_fwhm) > 1e-8 * max(expected_fwhm, 1e-12):
            if fitter == "preview_nnls":
                shared_ok = True
            else:
                shared_ok = shared_ok and abs(mean_fwhm - expected_fwhm) <= 1e-6 * max(expected_fwhm, 1.0)
                if not shared_ok:
                    warnings.append("Fitted FWHM is not identical across components.")
    if design_cond is not None and design_cond > 1e8:
        warnings.append(f"Design matrix is ill-conditioned (cond={design_cond:.3g}).")
    if n_empty and components:
        warnings.append(f"{n_empty}/{len(components)} components are empty or near-zero.")
    if fitter == "nnls_fallback":
        warnings.append("Fityk fit was not used; result is the NNLS fallback.")
    return QCReport(
        rwp=float(metrics.get("rwp", math.nan)),
        r_factor=float(metrics.get("r_factor", math.nan)),
        rms_residual=float(metrics.get("rms_residual", math.nan)),
        wssr=metrics.get("wssr"),
        shared_fwhm_ok=shared_ok,
        fwhm_relative_span=float(span) if np.isfinite(span) else math.nan,
        empty_fraction=(n_empty / len(components)) if components else math.nan,
        design_cond=design_cond,
        n_peaks=len(components),
        n_empty=n_empty,
        fitter=fitter,
        warnings=warnings,
    )
