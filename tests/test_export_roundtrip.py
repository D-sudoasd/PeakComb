from __future__ import annotations

from pathlib import Path

import pandas as pd

from peakcomb.export import export_run
from peakcomb.models import BackgroundKind, FitMode, Session
from peakcomb.pipeline import official_fit
from peakcomb.seed import seed_peaks
from peakcomb.synthesize import synthetic_domain_spectrum
from peakcomb.io import write_xy, load_spectrum


def test_export_bundle_has_identical_fwhm(tmp_path: Path):
    spectrum, truth = synthetic_domain_spectrum(noise=2.0)
    write_xy(tmp_path / "spec.xy", spectrum.x, spectrum.y)
    loaded = load_spectrum(tmp_path / "spec.xy", axis="d_a", energy_kev=83.0)
    session = Session(
        input_path=str(tmp_path / "spec.xy"),
        xmin=2.22,
        xmax=2.42,
        mode=FitMode.DOMAIN,
        fwhm=float(truth["fwhm"]),
        background=BackgroundKind.LINEAR,
        n_peaks=5,
        seed_method="equal",
    )
    seed_peaks(session, loaded)
    result = official_fit(session, loaded, tmp_path / "work", force_nnls=True)
    bundle = export_run(result, tmp_path / "run", session)
    components = pd.read_csv(tmp_path / "run" / "components.csv")
    assert components["fwhm"].nunique() == 1
    span = (components["fwhm"].max() - components["fwhm"].min()) / components["fwhm"].mean()
    assert span < 1e-8
    curves = pd.read_csv(tmp_path / "run" / "curves.csv")
    assert "y_fit" in curves.columns
    assert (tmp_path / "run" / "session.yaml").is_file()
    assert (tmp_path / "run" / "qc.json").is_file()
    assert (tmp_path / "run" / "figures" / "overlay.png").is_file()
    assert (tmp_path / "run" / "origin" / "peakcomb.xlsx").is_file()
    assert any(path.name == "components.csv" for path in bundle.paths)
    assert result.fitter in {"preview_nnls", "nnls_fallback", "fityk"}


def test_loader_reads_written_xy(tmp_path: Path):
    x = [2.2, 2.3, 2.4, 2.5]
    y = [1.0, 3.0, 2.0, 1.1]
    path = write_xy(tmp_path / "a.xy", x, y)
    spec = load_spectrum(path)
    assert spec.x.size == 4
    assert spec.y[1] == 3.0
