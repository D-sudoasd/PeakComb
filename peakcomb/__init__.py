"""PeakComb: shared-FWHM comb decomposition of broad XRD peaks."""

from __future__ import annotations

__version__ = "0.1.0"

from .models import (
    AXIS_D_A,
    AXIS_Q_INV_A,
    AXIS_TWO_THETA_DEG,
    BackgroundKind,
    Component,
    FitMode,
    FitResult,
    PeakSeed,
    Session,
    Spectrum,
)

__all__ = [
    "AXIS_D_A",
    "AXIS_Q_INV_A",
    "AXIS_TWO_THETA_DEG",
    "BackgroundKind",
    "Component",
    "FitMode",
    "FitResult",
    "PeakSeed",
    "Session",
    "Spectrum",
    "__version__",
]
