"""Non-negative least-squares preview with a frozen shared FWHM."""

from __future__ import annotations

import numpy as np
from scipy.optimize import nnls

from .analysis import components_from_arrays, r_metrics, slice_window
from .kernel import background_columns, peak_design_matrix
from .limits import MAX_COND_ELEMENTS, MAX_PREVIEW_POINTS
from .models import FitResult, PeakSeed, Session, Spectrum
from .mechanics import analyze_components
from .qc import build_qc


def _background_kind(session: Session) -> str:
    return session.background.value if hasattr(session.background, "value") else str(session.background)


def _downsample(x: np.ndarray, y: np.ndarray, sigma: np.ndarray | None) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    if x.size <= MAX_PREVIEW_POINTS:
        return x, y, sigma
    stride = int(np.ceil(x.size / MAX_PREVIEW_POINTS))
    sl = slice(None, None, stride)
    return x[sl], y[sl], None if sigma is None else sigma[sl]


def nnls_preview(
    session: Session,
    spectrum: Spectrum,
    seeds: list[PeakSeed] | None = None,
) -> FitResult:
    used_seeds = seeds if seeds is not None else list(session.peaks)
    if not used_seeds:
        raise ValueError("No peak seeds. Run seeding first.")
    xmin, xmax = session.window(spectrum)
    x, y, sigma = slice_window(spectrum, xmin, xmax)
    x, y, sigma = _downsample(x, y, sigma)
    centers = np.array([seed.center for seed in used_seeds], dtype=float)
    peak_g = peak_design_matrix(x, centers, hwhm=session.hwhm, shape=session.shape)
    bkg_g = background_columns(x, _background_kind(session), x_mid=0.5 * (xmin + xmax))
    if bkg_g is None:
        design = peak_g
        n_bkg = 0
    else:
        design = np.column_stack([peak_g, bkg_g])
        n_bkg = bkg_g.shape[1]
    coeff, _residual = nnls(design, y)
    heights = coeff[: len(used_seeds)]
    bkg_coeff = coeff[len(used_seeds) :]
    y_bkg = np.zeros_like(x) if n_bkg == 0 else bkg_g @ bkg_coeff
    peak_curves = (peak_g * heights).T
    y_fit = peak_curves.sum(axis=0) + y_bkg
    metrics = r_metrics(y, y_fit, sigma)
    components = components_from_arrays(session, used_seeds, heights)
    cond = None
    if design.size <= MAX_COND_ELEMENTS:
        try:
            cond = float(np.linalg.cond(design))
        except np.linalg.LinAlgError:
            cond = None
    qc = build_qc(
        components,
        metrics,
        fitter="preview_nnls",
        design_cond=cond,
        empty_area_frac=session.empty_area_frac,
        expected_fwhm=session.fwhm,
    )
    return FitResult(
        components=components,
        x=x,
        y_obs=y,
        y_bkg=y_bkg,
        y_fit=y_fit,
        residual=y - y_fit,
        peak_curves=peak_curves,
        qc=qc,
        fitter="preview_nnls",
        session=session,
        mechanics=analyze_components(components, session),
    )
