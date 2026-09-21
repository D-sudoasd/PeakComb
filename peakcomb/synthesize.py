"""Synthetic spectra used by tests and the example folder."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .kernel import background_columns, pseudo_voigt
from .models import AXIS_D_A, Spectrum


def synthetic_domain_spectrum(
    *,
    noise: float = 8.0,
    rng: np.random.Generator | None = None,
) -> tuple[Spectrum, dict]:
    """Five equal-FWHM PseudoVoigts on a linear background, d in 2.20–2.44 Å."""
    generator = rng or np.random.default_rng(0)
    x = np.linspace(2.18, 2.46, 561)
    truth = {
        "centers": np.array([2.28, 2.30, 2.32, 2.34, 2.36]),
        "heights": np.array([220.0, 820.0, 1500.0, 910.0, 280.0]),
        "fwhm": 0.012,
        "shape": 0.5,
        "background": (40.0, 15.0),
    }
    hwhm = 0.5 * truth["fwhm"]
    y = np.zeros_like(x)
    for center, height in zip(truth["centers"], truth["heights"], strict=True):
        y += pseudo_voigt(x, height=float(height), center=float(center), hwhm=hwhm, shape=truth["shape"])
    bkg = background_columns(x, "linear", x_mid=float(x.mean()))
    y = y + bkg @ np.array(truth["background"])
    y = y + generator.normal(0.0, noise, size=x.size)
    spectrum = Spectrum(x=x, y=np.clip(y, 0.0, None), axis=AXIS_D_A, wavelength_a=0.14938, energy_kev=83.0)
    return spectrum, truth


def synthetic_distribution_spectrum(
    *,
    noise: float = 6.0,
    rng: np.random.Generator | None = None,
) -> tuple[Spectrum, dict]:
    """A Gaussian d-distribution convolved with a fixed-width kernel."""
    generator = rng or np.random.default_rng(1)
    x = np.linspace(2.22, 2.42, 801)
    fwhm = 0.010
    hwhm = 0.5 * fwhm
    d0 = 2.32
    sigma = 0.012
    generating = np.linspace(2.24, 2.40, 81)
    weights = np.exp(-0.5 * ((generating - d0) / sigma) ** 2)
    weights = weights / weights.max() * 180.0
    y = np.zeros_like(x)
    for center, height in zip(generating, weights, strict=True):
        y += pseudo_voigt(x, height=float(height), center=float(center), hwhm=hwhm, shape=0.5)
    y = y + 25.0 + generator.normal(0.0, noise, size=x.size)
    spectrum = Spectrum(x=x, y=np.clip(y, 0.0, None), axis=AXIS_D_A, wavelength_a=0.14938, energy_kev=83.0)
    truth = {
        "d0": d0,
        "sigma": sigma,
        "fwhm": fwhm,
        "generating_centers": generating,
        "generating_heights": weights,
    }
    return spectrum, truth


def write_example_files(directory: str | Path | None = None) -> list[Path]:
    from .io import write_xy

    root = Path(directory) if directory is not None else Path(__file__).resolve().parent.parent / "examples"
    root.mkdir(parents=True, exist_ok=True)
    domain, _ = synthetic_domain_spectrum(noise=8.0)
    distribution, _ = synthetic_distribution_spectrum(noise=6.0)
    return [
        write_xy(root / "synthetic_domain.xy", domain.x, domain.y),
        write_xy(root / "synthetic_distribution.xy", distribution.x, distribution.y),
    ]


if __name__ == "__main__":
    for path in write_example_files():
        print(path)
