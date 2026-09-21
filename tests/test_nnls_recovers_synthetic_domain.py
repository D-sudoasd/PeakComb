from __future__ import annotations

import numpy as np

from peakcomb.models import BackgroundKind, FitMode, PeakSeed, SeedMethod, Session
from peakcomb.preview import nnls_preview
from peakcomb.synthesize import synthetic_distribution_spectrum, synthetic_domain_spectrum


def test_nnls_recovers_five_equal_width_peaks():
    spectrum, truth = synthetic_domain_spectrum(noise=1.0)
    session = Session(
        axis="d_a",
        xmin=2.22,
        xmax=2.42,
        mode=FitMode.DOMAIN,
        fwhm=float(truth["fwhm"]),
        shape=0.5,
        background=BackgroundKind.LINEAR,
        n_peaks=5,
        seed_method=SeedMethod.MANUAL,
        peaks=[PeakSeed(center=float(c), height=float(h)) for c, h in zip(truth["centers"], truth["heights"])],
    )
    result = nnls_preview(session, spectrum)
    recovered = np.array([item.center for item in result.components])
    np.testing.assert_allclose(recovered, truth["centers"], atol=0.15 * truth["fwhm"])
    assert len(result.components) == 5
    height_rel = np.array([item.height for item in result.components]) / truth["heights"]
    np.testing.assert_allclose(height_rel, np.ones(5), rtol=0.12)
    true_frac = truth["heights"] / truth["heights"].sum()
    rec_frac = np.array([item.area_frac for item in result.components])
    np.testing.assert_allclose(rec_frac, true_frac, atol=0.08)
    assert result.qc.rwp < 0.05
    assert result.qc.shared_fwhm_ok


def test_distribution_comb_recovers_gaussian_d_weights():
    spectrum, truth = synthetic_distribution_spectrum(noise=1.0)
    session = Session(
        axis="d_a",
        xmin=2.24,
        xmax=2.40,
        mode=FitMode.DISTRIBUTION,
        fwhm=float(truth["fwhm"]),
        shape=0.5,
        background=BackgroundKind.CONSTANT,
        spacing_factor=0.45,
        seed_method="grid",
        d0=float(truth["d0"]),
    )
    from peakcomb.seed import seed_peaks

    seed_peaks(session, spectrum)
    result = nnls_preview(session, spectrum)
    centers = np.array([item.center for item in result.components])
    weights = np.array([item.area_frac for item in result.components])
    true = np.exp(-0.5 * ((centers - truth["d0"]) / truth["sigma"]) ** 2)
    true = true / true.sum()
    corr = np.corrcoef(weights, true)[0, 1]
    assert corr > 0.9
    assert result.qc.n_peaks >= 8
