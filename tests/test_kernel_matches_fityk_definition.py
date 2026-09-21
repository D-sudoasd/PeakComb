from __future__ import annotations

import math

import numpy as np

from peakcomb.kernel import (
    GAUSS_AREA_FACTOR,
    LORENTZ_AREA_FACTOR,
    fwhm_from_hwhm,
    pseudo_voigt,
    pseudo_voigt_area,
)


def test_gaussian_half_maximum_is_hwhm():
    hwhm = 0.01
    x = np.array([0.0, hwhm])
    y = pseudo_voigt(x, height=10.0, center=0.0, hwhm=hwhm, shape=0.0)
    assert y[0] == 10.0
    np.testing.assert_allclose(y[1], 5.0, rtol=1e-12)


def test_lorentzian_half_maximum_is_hwhm():
    hwhm = 0.02
    y = pseudo_voigt(np.array([0.0, hwhm]), height=8.0, center=0.0, hwhm=hwhm, shape=1.0)
    np.testing.assert_allclose(y, [8.0, 4.0], rtol=1e-12)


def test_area_matches_analytic_factors():
    height, hwhm, shape = 12.0, 0.015, 0.4
    analytic = pseudo_voigt_area(height, hwhm, shape)
    expected = height * hwhm * ((1.0 - shape) * GAUSS_AREA_FACTOR + shape * LORENTZ_AREA_FACTOR)
    np.testing.assert_allclose(analytic, expected)
    x = np.linspace(-0.4, 0.4, 40001)
    y_g = pseudo_voigt(x, height=height, center=0.0, hwhm=hwhm, shape=0.0)
    trapz = getattr(np, "trapezoid", np.trapz)
    np.testing.assert_allclose(trapz(y_g, x), pseudo_voigt_area(height, hwhm, 0.0), rtol=2e-4)
    assert math.isclose(fwhm_from_hwhm(hwhm), 2.0 * hwhm)


def test_kernel_is_linear_in_height():
    x = np.linspace(-0.1, 0.1, 51)
    y1 = pseudo_voigt(x, height=1.0, center=0.0, hwhm=0.01, shape=0.5)
    y3 = pseudo_voigt(x, height=3.0, center=0.0, hwhm=0.01, shape=0.5)
    np.testing.assert_allclose(y3, 3.0 * y1)
