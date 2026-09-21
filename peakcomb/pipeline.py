"""Orchestrate seeding, NNLS preview, optional Fityk fit, and analysis."""

from __future__ import annotations

from pathlib import Path

from .analysis import evaluate_model, r_metrics, slice_window
from .fityk_runtime import (
    FitykRuntimeError,
    components_from_peaks,
    detect_cfityk,
    parse_peaks_text,
    run_cfityk,
)
from .io import write_xy
from .lua import render_fit_lua
from .models import FitMode, FitResult, PeakSeed, Session, Spectrum
from .preview import nnls_preview
from .mechanics import analyze_components
from .qc import build_qc
from .seed import seed_peaks


DISTRIBUTION_FITYK_PEAK_LIMIT = 30


def ensure_seeds(session: Session, spectrum: Spectrum) -> list[PeakSeed]:
    if session.peaks:
        xmin, xmax = session.window(spectrum)
        kept = [peak for peak in session.peaks if xmin <= peak.center <= xmax]
        if kept:
            session.peaks = kept
            session.n_peaks = len(kept)
            return kept
    return seed_peaks(session, spectrum)


def preview_fit(session: Session, spectrum: Spectrum) -> FitResult:
    seeds = ensure_seeds(session, spectrum)
    return nnls_preview(session, spectrum, seeds)


def _should_try_fityk(session: Session, seeds: list[PeakSeed], cfityk_path: str | None) -> bool:
    if detect_cfityk(cfityk_path or session.cfityk_path or None) is None:
        return False
    if session.mode is FitMode.DISTRIBUTION and len(seeds) > DISTRIBUTION_FITYK_PEAK_LIMIT:
        return False
    return True


def official_fit(
    session: Session,
    spectrum: Spectrum,
    work_dir: str | Path,
    *,
    force_nnls: bool = False,
) -> FitResult:
    seeds = ensure_seeds(session, spectrum)
    preview = nnls_preview(session, spectrum, seeds)
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    if force_nnls or not _should_try_fityk(session, seeds, session.cfityk_path):
        preview.fitter = "nnls_fallback" if force_nnls or detect_cfityk(session.cfityk_path or None) is None else preview.fitter
        if preview.fitter != "preview_nnls":
            preview.qc = build_qc(
                preview.components,
                {
                    "rwp": preview.qc.rwp,
                    "r_factor": preview.qc.r_factor,
                    "rms_residual": preview.qc.rms_residual,
                    "wssr": preview.qc.wssr,
                },
                fitter=preview.fitter,
                design_cond=preview.qc.design_cond,
                empty_area_frac=session.empty_area_frac,
                expected_fwhm=session.fwhm,
                extra_warnings=list(preview.qc.warnings),
            )
        preview.work_dir = work
        return preview

    xmin, xmax = session.window(spectrum)
    x, y, sigma = slice_window(spectrum, xmin, xmax)
    input_xy = write_xy(work / "input.xy", x, y)
    prefix = work / "fit"
    lua_text = render_fit_lua(
        session,
        seeds,
        input_xy=input_xy,
        output_prefix=prefix,
        xmin=xmin,
        xmax=xmax,
    )
    script_path = work / "fit.lua"
    script_path.write_text(lua_text, encoding="utf-8")
    try:
        run_cfityk(script_path, cfityk_path=session.cfityk_path or None)
        peaks_path = prefix.with_suffix(".peaks")
        if not peaks_path.is_file():
            raise FitykRuntimeError("cfityk finished without writing fit.peaks.")
        peaks_text = peaks_path.read_text(encoding="utf-8", errors="replace")
        parsed = parse_peaks_text(peaks_text)
        components = components_from_peaks(session, parsed, seeds)
        rebuilt_seeds = [
            PeakSeed(center=item.center, height=item.height, label=item.label)
            for item in components
        ]
        heights = [item.height for item in components]
        y_bkg, y_fit, peak_curves = evaluate_model(
            x,
            rebuilt_seeds,
            heights,
            hwhm=components[0].hwhm if components else session.hwhm,
            shape=components[0].shape if components else session.shape,
            background=preview.y_bkg,
        )
        # Recompute background-free peak sum; keep NNLS background as the local baseline.
        y_fit = peak_curves.sum(axis=0) + preview.y_bkg
        y_bkg = preview.y_bkg
        metrics = r_metrics(y, y_fit, sigma)
        log_path = prefix.with_suffix(".log")
        log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
        qc = build_qc(
            components,
            metrics,
            fitter="fityk",
            design_cond=preview.qc.design_cond,
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
            fitter="fityk",
            session=session,
            lua_text=lua_text,
            peaks_text=peaks_text,
            log_text=log_text,
            work_dir=work,
            mechanics=analyze_components(components, session),
        )
    except (FitykRuntimeError, OSError, TimeoutError) as exc:
        fallback = preview
        fallback.fitter = "nnls_fallback"
        fallback.lua_text = lua_text
        fallback.work_dir = work
        fallback.qc = build_qc(
            fallback.components,
            {
                "rwp": fallback.qc.rwp,
                "r_factor": fallback.qc.r_factor,
                "rms_residual": fallback.qc.rms_residual,
                "wssr": fallback.qc.wssr,
            },
            fitter="nnls_fallback",
            design_cond=fallback.qc.design_cond,
            empty_area_frac=session.empty_area_frac,
            expected_fwhm=session.fwhm,
            extra_warnings=[*fallback.qc.warnings, f"Fityk failed: {exc}"],
        )
        return fallback
