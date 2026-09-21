from __future__ import annotations

from pathlib import Path

from peakcomb.limits import MAX_COMB_PEAKS
from peakcomb.lua import render_fit_lua
from peakcomb.seed import grid_centers
from peakcomb.models import FitMode, PeakSeed, Session
from peakcomb.fityk_runtime import parse_peaks_text


def _session(mode: FitMode) -> Session:
    return Session(
        mode=mode,
        fwhm=0.012,
        shape=0.5,
        refine_shared_fwhm=False,
        n_peaks=3,
        peaks=[
            PeakSeed(center=2.30, height=100),
            PeakSeed(center=2.32, height=200),
            PeakSeed(center=2.34, height=150),
        ],
    )


def test_lua_uses_one_shared_hwhm_variable(tmp_path: Path):
    session = _session(FitMode.DOMAIN)
    lua = render_fit_lua(
        session,
        session.peaks,
        input_xy=tmp_path / "input.xy",
        output_prefix=tmp_path / "fit",
        xmin=2.25,
        xmax=2.40,
    )
    assert lua.count("$hwhm = ") == 1
    assert lua.count("$hwhm = ~") == 0
    assert lua.count("PseudoVoigt(") == 3
    assert lua.count(", $hwhm, $shape)") == 3
    assert "~2.3" in lua  # free centers in Domain
    assert "set fitting_method = levenberg_marquardt" in lua


def test_distribution_lua_locks_centers(tmp_path: Path):
    session = _session(FitMode.DISTRIBUTION)
    lua = render_fit_lua(
        session,
        session.peaks,
        input_xy=tmp_path / "input.xy",
        output_prefix=tmp_path / "fit",
        xmin=2.25,
        xmax=2.40,
    )
    assert "PseudoVoigt(~100" in lua or "PseudoVoigt(~100.0" in lua
    assert "~2.30" not in lua
    assert ", 2.3" in lua
    session.refine_shared_fwhm = True
    lua_free = render_fit_lua(
        session,
        session.peaks,
        input_xy=tmp_path / "input.xy",
        output_prefix=tmp_path / "fit",
        xmin=2.25,
        xmax=2.40,
    )
    assert "$hwhm = ~" in lua_free


def test_grid_centers_caps_wide_windows():
    centers = grid_centers(0.4, 12.0, fwhm=0.012, spacing_factor=0.4)
    assert 2 <= centers.size <= MAX_COMB_PEAKS


def test_parse_fityk_peaks_block():
    text = """# PeakType\tCenter\tHeight\tArea\tFWHM\tparameters...
%_1  PseudoVoigt\t2.30\t100\t2\t0.012\t 100 2.30 0.006 0.5
%_2  PseudoVoigt\t2.32\t200\t4\t0.012\t 200 2.32 0.006 0.5
"""
    parsed = parse_peaks_text(text)
    assert len(parsed) == 2
    assert parsed[0].center == 2.30
    assert parsed[1].fwhm == 0.012
