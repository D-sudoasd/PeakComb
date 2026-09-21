from __future__ import annotations

import numpy as np

from peakcomb.mechanics import analyze_components, hydrostatic_modulus_gpa
from peakcomb.seed import symmetric_centers
from peakcomb.models import Component, Session


def _components(ds, weights, fwhm=0.01):
    items = []
    for i, (d, w) in enumerate(zip(ds, weights, strict=True)):
        items.append(
            Component(
                index=i,
                label=f"P{i+1}",
                center=d,
                height=w,
                area=w,
                area_frac=w,
                fwhm=fwhm,
                hwhm=0.5 * fwhm,
                shape=0.5,
                d_a=d,
            )
        )
    return items


def test_symmetric_centers_mirror_about_peak():
    centers = symmetric_centers(2.32907, 9, 0.004)
    assert centers.size == 9
    assert abs(centers[4] - 2.32907) < 1e-12
    np.testing.assert_allclose(centers + centers[::-1], 2.0 * 2.32907)


def test_area_weighted_d0_makes_mean_strain_zero():
    comps = _components([2.32, 2.33, 2.34], [0.2, 0.5, 0.3])
    session = Session(d0_source="area_weighted", e_gpa=80.0, nu=0.33)
    report = analyze_components(comps, session)
    assert abs(report.mean_strain) < 1e-12
    assert report.tautology
    assert report.balance_ok


def test_peak_max_d0_and_mean_strain():
    comps = _components([2.32, 2.33, 2.34], [0.2, 0.6, 0.2])
    session = Session(d0_source="peak_max", e_gpa=80.0, nu=0.33)
    report = analyze_components(comps, session)
    assert abs(report.d0 - 2.33) < 1e-12
    expected = 0.2 * (2.32 - 2.33) / 2.33 + 0.2 * (2.34 - 2.33) / 2.33
    assert abs(report.mean_strain - expected) < 1e-12


def test_hydrostatic_stress_example_0p1_percent():
    modulus = hydrostatic_modulus_gpa(80.0, 0.33)
    assert abs(modulus - 80.0 / 0.34) < 1e-9
    comps = _components([2.332329, 2.329], [1.0, 0.0])
    comps[0].area_frac = 1.0
    comps[1].area_frac = 0.0
    session = Session(d0=2.329, d0_source="manual", e_gpa=80.0, nu=0.33)
    # 2.332329/2.329 - 1 ≈ 0.00143, use exact 0.1%:
    comps = _components([2.329 * 1.001], [1.0])
    report = analyze_components(comps, session)
    assert abs(report.mean_strain - 0.001) < 1e-12
    assert abs(report.mean_stress_mpa - 1000.0 * modulus * 0.001) < 1e-6
    assert 230 < report.mean_stress_mpa < 240
