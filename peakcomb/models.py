"""Session, spectrum, and fit-result dataclasses."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np

AXIS_D_A = "d_a"
AXIS_Q_INV_A = "q_inv_a"
AXIS_TWO_THETA_DEG = "two_theta_deg"
VALID_AXES = (AXIS_D_A, AXIS_Q_INV_A, AXIS_TWO_THETA_DEG)

DEFAULT_WAVELENGTH_A = 1.5406
DEFAULT_ENERGY_KEV = 83.0
X_RAY_KEV_ANGSTROM = 12.398419843


class FitMode(str, Enum):
    DOMAIN = "domain"
    DISTRIBUTION = "distribution"


class BackgroundKind(str, Enum):
    NONE = "none"
    CONSTANT = "constant"
    LINEAR = "linear"
    QUADRATIC = "quadratic"


class SeedMethod(str, Enum):
    EQUAL = "equal"
    MAXIMA = "maxima"
    MANUAL = "manual"
    GRID = "grid"
    SYMMETRIC = "symmetric"


FITTERS = ("preview_nnls", "fityk", "nnls_fallback")


@dataclass
class Spectrum:
    x: np.ndarray
    y: np.ndarray
    sigma: np.ndarray | None = None
    path: Path | None = None
    axis: str = AXIS_D_A
    wavelength_a: float = DEFAULT_WAVELENGTH_A
    energy_kev: float | None = None

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        if self.sigma is not None:
            self.sigma = np.asarray(self.sigma, dtype=float)
        if self.x.ndim != 1 or self.y.ndim != 1 or self.x.size != self.y.size:
            raise ValueError("Spectrum x and y must be 1-D arrays of equal length.")
        if self.sigma is not None and self.sigma.shape != self.y.shape:
            raise ValueError("Spectrum sigma must match y.")


@dataclass
class PeakSeed:
    center: float
    height: float = 1.0
    visible: bool = True
    label: str = ""


@dataclass
class Session:
    input_path: str = ""
    axis: str = AXIS_D_A
    wavelength_a: float = DEFAULT_WAVELENGTH_A
    energy_kev: float | None = DEFAULT_ENERGY_KEV
    xmin: float | None = None
    xmax: float | None = None
    comb_xmin: float | None = None
    comb_xmax: float | None = None
    baseline_method: str = "none"
    mode: FitMode = FitMode.DOMAIN
    fwhm: float = 0.012
    refine_shared_fwhm: bool = False
    shape: float = 0.5
    refine_shared_shape: bool = False
    background: BackgroundKind = BackgroundKind.LINEAR
    n_peaks: int = 6
    seed_method: SeedMethod = SeedMethod.EQUAL
    spacing_factor: float = 0.4
    center_slack_factor: float = 0.35
    d0: float | None = None
    d0_source: str = "peak_max"
    e_gpa: float = 80.0
    nu: float = 0.33
    k_gpa: float | None = None
    compute_stress: bool = True
    mean_strain_tol: float = 5.0e-4
    fitting_method: str = "levenberg_marquardt"
    max_wssr_evaluations: int = 800
    max_fitting_time: float = 30.0
    empty_area_frac: float = 1.0e-4
    cfityk_path: str = ""
    peaks: list[PeakSeed] = field(default_factory=list)

    def __post_init__(self) -> None:
        if isinstance(self.mode, str):
            self.mode = FitMode(self.mode)
        if isinstance(self.background, str):
            self.background = BackgroundKind(self.background)
        if isinstance(self.seed_method, str):
            self.seed_method = SeedMethod(self.seed_method)
        if self.axis not in VALID_AXES:
            raise ValueError(f"axis must be one of {VALID_AXES}; got {self.axis!r}.")
        if not np.isfinite(self.fwhm) or self.fwhm <= 0.0:
            raise ValueError("fwhm must be a positive finite number.")
        if not 0.0 <= float(self.shape) <= 1.0:
            raise ValueError("shape must be in [0, 1].")
        if self.n_peaks < 1:
            raise ValueError("n_peaks must be >= 1.")

    @property
    def hwhm(self) -> float:
        return 0.5 * float(self.fwhm)

    def window(self, spectrum: Spectrum) -> tuple[float, float]:
        xmin = float(spectrum.x.min()) if self.xmin is None else float(self.xmin)
        xmax = float(spectrum.x.max()) if self.xmax is None else float(self.xmax)
        if xmax <= xmin:
            raise ValueError("xmax must be greater than xmin.")
        return xmin, xmax

    def comb_window(self, spectrum: Spectrum) -> tuple[float, float]:
        xmin, xmax = self.window(spectrum)
        cmin = xmin if self.comb_xmin is None else float(self.comb_xmin)
        cmax = xmax if self.comb_xmax is None else float(self.comb_xmax)
        cmin = max(cmin, xmin)
        cmax = min(cmax, xmax)
        if cmax <= cmin:
            return xmin, xmax
        return cmin, cmax

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["mode"] = self.mode.value
        payload["background"] = self.background.value
        payload["seed_method"] = self.seed_method.value
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> Session:
        allowed = {item.name for item in fields(cls)}
        data = {key: value for key, value in payload.items() if key in allowed}
        peaks = []
        for item in data.get("peaks") or []:
            if isinstance(item, PeakSeed):
                peaks.append(item)
            else:
                peaks.append(PeakSeed(**item))
        data["peaks"] = peaks
        return cls(**data)


@dataclass
class Component:
    index: int
    label: str
    center: float
    height: float
    area: float
    area_frac: float
    fwhm: float
    hwhm: float
    shape: float
    d_a: float | None = None
    q_inv_a: float | None = None
    two_theta_deg: float | None = None
    strain: float | None = None
    sigma_h_mpa: float | None = None
    visible: bool = True


@dataclass
class QCReport:
    rwp: float
    r_factor: float
    rms_residual: float
    wssr: float | None
    shared_fwhm_ok: bool
    fwhm_relative_span: float
    empty_fraction: float
    design_cond: float | None
    n_peaks: int
    n_empty: int
    fitter: str
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class FitResult:
    components: list[Component]
    x: np.ndarray
    y_obs: np.ndarray
    y_bkg: np.ndarray
    y_fit: np.ndarray
    residual: np.ndarray
    peak_curves: np.ndarray
    qc: QCReport
    fitter: str
    session: Session
    lua_text: str = ""
    peaks_text: str = ""
    log_text: str = ""
    work_dir: Path | None = None
    mechanics: Any = None
