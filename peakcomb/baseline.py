"""Windowed baseline subtraction. pybaselines when available, rubberband otherwise."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial import ConvexHull

BASELINE_METHODS = (
    "none",
    "asls",
    "arpls",
    "airpls",
    "snip",
    "modpoly",
    "rubberband",
)

_PYBASELINES_METHODS = {
    "asls": {"lam": 1.0e5, "p": 0.01},
    "arpls": {"lam": 1.0e5},
    "airpls": {"lam": 1.0e6},
    "snip": {},
    "modpoly": {"poly_order": 2},
}


class BaselineError(ValueError):
    """Raised when a baseline cannot be estimated."""


@dataclass
class BaselineResult:
    method: str
    baseline: np.ndarray
    corrected: np.ndarray
    backend: str
    params: dict


def pybaselines_available() -> bool:
    try:
        import pybaselines  # noqa: F401
    except ImportError:
        return False
    return True


def rubberband_baseline(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Lower convex-hull envelope; no pybaselines required."""
    if x.size < 3:
        raise BaselineError("Rubberband baseline needs at least 3 points.")
    points = np.column_stack([x, y])
    try:
        verts = np.asarray(ConvexHull(points).vertices, dtype=int)
    except Exception as exc:
        raise BaselineError(f"Rubberband hull failed: {exc}") from exc
    order = np.argsort(x[verts])
    vx = x[verts][order]
    vy = y[verts][order]
    chord = vy[0] + (vy[-1] - vy[0]) * (vx - vx[0]) / (vx[-1] - vx[0] + 1e-30)
    lower = np.flatnonzero(vy <= chord + 1e-12 * max(abs(float(np.max(y))), 1.0))
    if lower.size < 2:
        lower = np.array([0, vx.size - 1], dtype=int)
    interpolator = interp1d(vx[lower], vy[lower], kind="linear", bounds_error=False, fill_value="extrapolate")
    return np.asarray(interpolator(x), dtype=float)


def _pybaselines_fit(x: np.ndarray, y: np.ndarray, method: str) -> tuple[np.ndarray, dict, str]:
    from pybaselines import Baseline

    params = dict(_PYBASELINES_METHODS[method])
    fitter = Baseline(x_data=x)
    func = getattr(fitter, method, None)
    if func is None:
        raise BaselineError(f"pybaselines has no method {method!r}.")
    baseline, info = func(y, **params)
    version = getattr(__import__("pybaselines"), "__version__", "unknown")
    return np.asarray(baseline, dtype=float), {"pybaselines": version, **params}, "pybaselines"


def estimate_baseline(x: np.ndarray, y: np.ndarray, method: str = "asls") -> BaselineResult:
    method = str(method or "none").strip().lower()
    if method not in BASELINE_METHODS:
        raise BaselineError(f"Unknown baseline method: {method!r}.")
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    if xv.size != yv.size or xv.size < 3:
        raise BaselineError("Baseline needs matching x, y with at least 3 points.")
    if method == "none":
        baseline = np.zeros_like(yv)
        return BaselineResult(method, baseline, yv.copy(), "none", {})
    if method == "rubberband" or (method in _PYBASELINES_METHODS and not pybaselines_available()):
        if method != "rubberband" and method in _PYBASELINES_METHODS:
            baseline = rubberband_baseline(xv, yv)
            return BaselineResult("rubberband", baseline, yv - baseline, "rubberband_fallback", {"requested": method})
        baseline = rubberband_baseline(xv, yv)
        return BaselineResult("rubberband", baseline, yv - baseline, "rubberband", {})
    baseline, params, backend = _pybaselines_fit(xv, yv, method)
    return BaselineResult(method, baseline, yv - baseline, backend, params)


def subtract_in_window(
    x: np.ndarray,
    y: np.ndarray,
    xmin: float,
    xmax: float,
    method: str,
) -> tuple[np.ndarray, np.ndarray, BaselineResult]:
    """Subtract a baseline inside [xmin, xmax]; leave points outside unchanged."""
    mask = (x >= xmin) & (x <= xmax)
    if int(np.count_nonzero(mask)) < 8:
        raise BaselineError("扣背底窗口里点数太少。")
    fitted = estimate_baseline(x[mask], y[mask], method)
    corrected = np.asarray(y, dtype=float).copy()
    baseline_full = np.full_like(corrected, np.nan, dtype=float)
    baseline_full[mask] = fitted.baseline
    corrected[mask] = fitted.corrected
    return corrected, baseline_full, fitted
