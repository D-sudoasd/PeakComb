"""Publication-style overlay, residual, and distribution figures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from .kernel import pseudo_voigt
from .models import FitResult

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

PUBLICATION_PALETTE = (
    "#3B6EA5",
    "#E09F3E",
    "#9B5DE5",
    "#00A5A8",
    "#D64545",
    "#2A9D8F",
    "#F4A261",
    "#6C8EAD",
    "#B56576",
    "#4A4E69",
)


def _style(ax) -> None:
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def _fine_grid(result: FitResult, n: int = 1200) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[np.ndarray]]:
    x = np.linspace(float(result.x.min()), float(result.x.max()), n)
    y_bkg = np.interp(x, result.x, result.y_bkg)
    curves: list[np.ndarray] = []
    for component in result.components:
        curves.append(
            pseudo_voigt(
                x,
                height=float(component.height),
                center=float(component.center),
                hwhm=float(component.hwhm),
                shape=float(component.shape),
            )
        )
    y_fit = y_bkg + (np.sum(curves, axis=0) if curves else 0.0)
    return x, y_bkg, y_fit, curves


@dataclass
class ViewOptions:
    show_fill: bool = True
    show_outlines: bool = True
    show_bkg: bool = True
    show_residual: bool = True
    data_as_points: bool = True
    fill_alpha: float = 0.40
    hide_frac: float = 0.012
    marker_size: float = 3.6
    show_strain: bool = True
    strain_as_stress: bool = False


def render_split(ax, ax_res, result: FitResult, options: ViewOptions | None = None) -> None:
    """Draw the publication split on existing axes. GUI and export share this."""
    opts = options or ViewOptions()
    x_fine, y_bkg, y_fit, curves = _fine_grid(result)
    kept = [
        (component, curve)
        for component, curve in zip(result.components, curves, strict=True)
        if component.visible and component.area_frac >= opts.hide_frac
    ]
    kept.sort(key=lambda item: item[0].center)
    for index, (component, curve) in enumerate(kept):
        color = PUBLICATION_PALETTE[index % len(PUBLICATION_PALETTE)]
        if opts.show_fill:
            ax.fill_between(x_fine, y_bkg, curve + y_bkg, color=color, alpha=opts.fill_alpha, linewidth=0, zorder=2)
        if opts.show_outlines:
            ax.plot(x_fine, curve + y_bkg, color=color, lw=0.9, zorder=3)
    if opts.show_bkg:
        ax.plot(x_fine, y_bkg, color="#8A8A8A", lw=0.8, ls="--", label="bkg", zorder=1)
    ax.plot(x_fine, y_fit, color="#111111", lw=1.55, label="sum", zorder=4)
    if opts.data_as_points:
        ax.plot(result.x, result.y_obs, "o", color="#111111", ms=opts.marker_size, mew=0.0, zorder=5, label="data")
    else:
        ax.plot(result.x, result.y_obs, color="#111111", lw=1.0, zorder=5, label="data")
    ax.set_ylabel("Intensity (counts/s)")
    ax.legend(loc="upper right", fontsize=8, frameon=False, handlelength=1.4)
    y_min = float(np.min(result.y_obs))
    ax.set_ylim(bottom=min(0.0, y_min * 0.98) if y_min < 0 else None)
    _style(ax)
    if ax_res is None:
        return
    ax_res.set_visible(opts.show_residual)
    if not opts.show_residual:
        return
    ax_res.axhline(0.0, color="#B0B0B0", lw=0.7)
    ax_res.plot(result.x, result.residual, "o-", color="#333333", ms=2.4, lw=0.6)
    span = max(float(np.max(np.abs(result.residual))), 1e-6)
    ax_res.set_ylim(-1.25 * span, 1.25 * span)
    ax_res.set_ylabel("Res.")
    ax_res.set_xlabel(_x_label(result))
    _style(ax_res)


def plot_overlay(result: FitResult, path: Path, *, dpi: int = 300, options: ViewOptions | None = None) -> Path:
    opts = options or ViewOptions()
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(4.8, 4.6),
        sharex=True,
        gridspec_kw={"height_ratios": [3.4, 1.0]},
    )
    ax, ax_res = axes
    render_split(ax, ax_res, result, opts)
    fig.tight_layout(h_pad=0.25)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def render_strain(ax, result: FitResult, *, as_stress: bool = False) -> None:
    live = [item for item in result.components if item.strain is not None and item.area_frac >= 0.012]
    live.sort(key=lambda item: float(item.strain))
    if not live:
        ax.set_visible(False)
        return
    ax.set_visible(True)
    xs = [100.0 * float(item.strain) for item in live]
    ys = [float(item.area_frac) for item in live]
    if as_stress and all(item.sigma_h_mpa is not None for item in live):
        xs = [float(item.sigma_h_mpa) for item in live]
        xlabel = r"$\sigma_h$ (MPa)"
    else:
        xlabel = r"$\varepsilon$ (%)"
    colors = [PUBLICATION_PALETTE[i % len(PUBLICATION_PALETTE)] for i in range(len(live))]
    ax.axvline(0.0, color="#B0B0B0", lw=0.7)
    ax.vlines(xs, 0.0, ys, colors=colors, lw=1.4)
    ax.scatter(xs, ys, c=colors, s=22, zorder=3)
    report = getattr(result, "mechanics", None)
    if report is not None:
        x_mean = 100.0 * report.mean_strain if not as_stress else report.mean_stress_mpa
        ax.axvline(x_mean, color="#111111", lw=0.9, ls="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$w_i$")
    ax.set_ylim(bottom=0.0)
    _style(ax)


def plot_strain(result: FitResult, path: Path, *, dpi: int = 240, as_stress: bool = False) -> Path:
    fig, ax = plt.subplots(figsize=(4.8, 2.4))
    render_strain(ax, result, as_stress=as_stress)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_distribution(result: FitResult, path: Path, *, dpi: int = 200) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    if result.session.d0 not in (None, 0.0) and all(item.strain is not None for item in result.components):
        xs = [float(item.strain) for item in result.components]
        xlabel = r"$(d - d_0)/d_0$"
    elif all(item.d_a is not None for item in result.components):
        xs = [float(item.d_a) for item in result.components]
        xlabel = r"$d$ (Å)"
    else:
        xs = [item.center for item in result.components]
        xlabel = _x_label(result)
    ys = [item.area_frac for item in result.components]
    colors = [PUBLICATION_PALETTE[i % len(PUBLICATION_PALETTE)] for i in range(len(xs))]
    width = 0.6 * result.session.fwhm if result.session.d0 in (None, 0.0) else 0.6 * np.median(np.diff(np.sort(xs))) if len(xs) > 1 else 0.001
    ax.bar(xs, ys, width=max(width, 1e-6), color=colors, edgecolor="#222222", linewidth=0.4, alpha=0.85)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("area_frac (intensity share)")
    _style(ax)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def plot_residual(result: FitResult, path: Path, *, dpi: int = 200) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.plot(result.x, result.residual, color="#222222", lw=1.0)
    ax.set_xlabel(_x_label(result))
    ax.set_ylabel("Residual")
    _style(ax)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def _x_label(result: FitResult) -> str:
    axis = result.session.axis
    if axis == "d_a":
        return r"$d$ (Å)"
    if axis == "q_inv_a":
        return r"$q$ (Å$^{-1}$)"
    if axis == "two_theta_deg":
        return r"$2\theta$ (deg)"
    return "x"
