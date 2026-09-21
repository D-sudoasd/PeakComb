"""Write a reproducible PeakComb run directory."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import __version__
from .figures import plot_distribution, plot_overlay, plot_residual, plot_strain
from .models import FitResult, Session

INTERPRETATION_BOUNDARY = (
    "area_frac (w_i) is the intensity share under a shared-FWHM PseudoVoigt comb; "
    "it is a volume-fraction proxy only if structure factor, multiplicity and texture are comparable. "
    "epsilon_i=(d_i-d0)/d0. Mean strain sum(w_i epsilon_i)~0 is a physical check only when d0 is independent "
    "(not the area-weighted mean). Hydrostatic stress sigma_h=E/(1-2nu)*epsilon assumes isotropic elasticity "
    "and that peak shifts are local hydrostatic elastic strain."
)


@dataclass
class ExportBundle:
    output_dir: Path
    paths: list[Path] = field(default_factory=list)


def _components_frame(result: FitResult) -> pd.DataFrame:
    rows = []
    for item in result.components:
        rows.append(
            {
                "index": item.index,
                "label": item.label,
                "center": item.center,
                "d_a": item.d_a,
                "q_inv_a": item.q_inv_a,
                "two_theta_deg": item.two_theta_deg,
                "height": item.height,
                "area": item.area,
                "area_frac": item.area_frac,
                "fwhm": item.fwhm,
                "hwhm": item.hwhm,
                "shape": item.shape,
                "strain": item.strain,
                "sigma_h_mpa": item.sigma_h_mpa,
                "visible": item.visible,
            }
        )
    return pd.DataFrame(rows)


def _curves_frame(result: FitResult) -> pd.DataFrame:
    data = {
        "x": result.x,
        "y_obs": result.y_obs,
        "y_bkg": result.y_bkg,
        "y_fit": result.y_fit,
        "residual": result.residual,
    }
    for index, component in enumerate(result.components):
        data[f"peak_{index + 1:02d}"] = result.peak_curves[index]
    return pd.DataFrame(data)


def _distribution_frame(result: FitResult) -> pd.DataFrame:
    frame = _components_frame(result)[
        ["index", "label", "center", "d_a", "strain", "sigma_h_mpa", "area", "area_frac", "fwhm"]
    ]
    return frame


def export_run(result: FitResult, output_dir: str | Path, session: Session | None = None) -> ExportBundle:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    figures = out / "figures"
    figures.mkdir(exist_ok=True)
    used_session = session or result.session
    bundle = ExportBundle(output_dir=out)

    session_path = out / "session.yaml"
    payload = used_session.to_dict()
    payload["peakcomb_version"] = __version__
    payload["exported_at"] = datetime.now(timezone.utc).isoformat()
    payload["interpretation_boundary"] = INTERPRETATION_BOUNDARY
    payload["fitter"] = result.fitter
    session_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    bundle.paths.append(session_path)

    components_path = out / "components.csv"
    _components_frame(result).to_csv(components_path, index=False)
    bundle.paths.append(components_path)

    curves_path = out / "curves.csv"
    _curves_frame(result).to_csv(curves_path, index=False)
    bundle.paths.append(curves_path)

    distribution_path = out / "distribution.csv"
    _distribution_frame(result).to_csv(distribution_path, index=False)
    bundle.paths.append(distribution_path)

    qc_path = out / "qc.json"
    qc_payload = result.qc.to_dict()
    qc_payload["interpretation_boundary"] = INTERPRETATION_BOUNDARY
    qc_payload["peakcomb_version"] = __version__
    if result.mechanics is not None:
        qc_payload["mechanics"] = result.mechanics.to_dict()
    qc_path.write_text(json.dumps(qc_payload, indent=2), encoding="utf-8")
    bundle.paths.append(qc_path)

    if result.lua_text:
        lua_path = out / "fit.lua"
        lua_path.write_text(result.lua_text, encoding="utf-8")
        bundle.paths.append(lua_path)
    if result.peaks_text:
        peaks_path = out / "fit.peaks"
        peaks_path.write_text(result.peaks_text, encoding="utf-8")
        bundle.paths.append(peaks_path)
    if result.log_text:
        log_path = out / "fit.log"
        log_path.write_text(result.log_text, encoding="utf-8")
        bundle.paths.append(log_path)

    overlay = plot_overlay(result, figures / "overlay.png")
    residual = plot_residual(result, figures / "residual.png")
    distribution = plot_distribution(result, figures / "distribution.png")
    if result.mechanics is not None:
        bundle.paths.append(plot_strain(result, figures / "strain.png"))
        if result.mechanics.compute_stress:
            bundle.paths.append(plot_strain(result, figures / "stress.png", as_stress=True))
    plot_overlay(result, figures / "overlay.pdf")
    plot_overlay(result, figures / "overlay.svg")
    plot_distribution(result, figures / "distribution.pdf")
    bundle.paths.extend([overlay, residual, distribution])

    origin_dir = out / "origin"
    origin_dir.mkdir(exist_ok=True)
    xlsx_path = origin_dir / "peakcomb.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        _components_frame(result).to_excel(writer, sheet_name="components", index=False)
        _curves_frame(result).to_excel(writer, sheet_name="curves", index=False)
        _distribution_frame(result).to_excel(writer, sheet_name="distribution", index=False)
        pd.DataFrame([result.qc.to_dict()]).to_excel(writer, sheet_name="qc", index=False)
    bundle.paths.append(xlsx_path)
    return bundle
