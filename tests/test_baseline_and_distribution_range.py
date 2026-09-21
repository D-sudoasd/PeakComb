from __future__ import annotations

import numpy as np

from peakcomb.baseline import estimate_baseline, pybaselines_available, subtract_in_window
from peakcomb.limits import MAX_COMB_PEAKS
from peakcomb.models import FitMode, SeedMethod, Session
from peakcomb.seed import seed_peaks
from peakcomb.synthesize import synthetic_distribution_spectrum, synthetic_domain_spectrum


def test_rubberband_or_asls_lowers_background():
    spectrum, _ = synthetic_domain_spectrum(noise=1.0)
    method = "asls" if pybaselines_available() else "rubberband"
    result = estimate_baseline(spectrum.x, spectrum.y, method)
    assert result.corrected.size == spectrum.y.size
    assert float(np.mean(result.corrected)) < float(np.mean(spectrum.y))
    peak = int(np.argmax(spectrum.y))
    assert result.corrected[peak] > np.percentile(result.corrected, 60)
    assert np.nanmax(result.corrected) > 0.0


def test_subtract_only_inside_window():
    spectrum, _ = synthetic_domain_spectrum(noise=1.0)
    method = "asls" if pybaselines_available() else "rubberband"
    corrected, baseline, _fitted = subtract_in_window(spectrum.x, spectrum.y, 2.26, 2.38, method)
    outside = (spectrum.x < 2.26) | (spectrum.x > 2.38)
    np.testing.assert_allclose(corrected[outside], spectrum.y[outside])
    inside = ~outside
    assert np.any(np.abs(corrected[inside] - spectrum.y[inside]) > 1.0)
    assert np.all(np.isnan(baseline[outside]))


def test_distribution_honors_comb_range_and_n():
    spectrum, _ = synthetic_distribution_spectrum(noise=1.0)
    session = Session(
        xmin=2.24,
        xmax=2.40,
        comb_xmin=2.28,
        comb_xmax=2.36,
        mode=FitMode.DISTRIBUTION,
        fwhm=0.010,
        n_peaks=11,
        seed_method=SeedMethod.EQUAL,
    )
    seeds = seed_peaks(session, spectrum)
    assert len(seeds) == 11
    centers = np.array([item.center for item in seeds])
    assert centers.min() >= 2.28 - 1e-9
    assert centers.max() <= 2.36 + 1e-9


def test_distribution_n_is_capped():
    spectrum, _ = synthetic_distribution_spectrum(noise=1.0)
    session = Session(
        xmin=2.24,
        xmax=2.40,
        mode=FitMode.DISTRIBUTION,
        fwhm=0.010,
        n_peaks=500,
        seed_method=SeedMethod.EQUAL,
    )
    seeds = seed_peaks(session, spectrum)
    assert len(seeds) == MAX_COMB_PEAKS
