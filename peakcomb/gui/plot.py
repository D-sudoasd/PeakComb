"""Matplotlib canvas with drag-to-select range, like Fityk."""

from __future__ import annotations

from typing import Callable

import matplotlib
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
from PySide6.QtWidgets import QVBoxLayout, QWidget

from ..figures import ViewOptions, render_split, render_strain
from ..limits import MAX_PREVIEW_POINTS
from ..models import FitResult, Spectrum

AXIS_LABELS = {
    "d_a": r"$d$ (Å)",
    "q_inv_a": r"$q$ (Å$^{-1}$)",
    "two_theta_deg": r"$2\theta$ (deg)",
}


def _thin(x: np.ndarray, *ys: np.ndarray) -> tuple[np.ndarray, ...]:
    if x.size <= MAX_PREVIEW_POINTS:
        return (x, *ys)
    stride = int(np.ceil(x.size / MAX_PREVIEW_POINTS))
    sl = slice(None, None, stride)
    return (x[sl], *(y[sl] for y in ys))


class PlotCanvas(QWidget):
    def __init__(
        self,
        on_range: Callable[[float, float], None] | None = None,
        on_click: Callable[[float, float], None] | None = None,
    ) -> None:
        super().__init__()
        self.on_range = on_range
        self.on_click = on_click
        self.tool = "select"
        self.span: SpanSelector | None = None
        self.view = ViewOptions()
        self.figure = Figure(figsize=(6.4, 6.2), dpi=140, facecolor="white", layout=None)
        self.ax_main, self.ax_res, self.ax_strain = self.figure.subplots(
            3, 1, sharex=False, height_ratios=[3.2, 0.85, 1.15]
        )
        self.ax_res.sharex(self.ax_main)
        self.figure.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.08, hspace=0.28)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar.setVisible(False)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, 1)
        self._draw_empty()

    def set_tool(self, tool: str) -> None:
        self.tool = tool
        if tool == "browse":
            self.span and self.span.set_active(False)
            self.toolbar.setVisible(True)
        else:
            self.span and self.span.set_active(True)
            self.toolbar.setVisible(False)

    def _install_span(self) -> None:
        self.span = SpanSelector(
            self.ax_main,
            self._span_selected,
            "horizontal",
            useblit=True,
            interactive=False,
            button=1,
            minspan=0.0,
            props={"alpha": 0.22, "facecolor": "#0b6f71"},
        )
        self.span.set_active(self.tool != "browse")

    def _span_selected(self, xmin: float, xmax: float) -> None:
        if xmax < xmin:
            xmin, xmax = xmax, xmin
        x0, x1 = self.ax_main.get_xlim()
        minspan = max(abs(x1 - x0) * 0.004, 1e-12)
        if xmax - xmin < minspan:
            if self.tool == "add" and self.on_click is not None:
                self.on_click(0.5 * (xmin + xmax), 1.0)
            return
        if self.tool == "select" and self.on_range is not None:
            self.on_range(float(xmin), float(xmax))

    def _style_axes(self) -> None:
        for ax in (self.ax_main, self.ax_res, self.ax_strain):
            ax.tick_params(direction="in", top=True, right=True)
            ax.set_facecolor("white")

    def _draw_empty(self) -> None:
        for ax in (self.ax_main, self.ax_res, self.ax_strain):
            ax.clear()
        self._style_axes()
        self.ax_main.set_ylabel("Intensity (counts/s)")
        self.ax_res.set_ylabel("Res.")
        self.ax_res.set_xlabel(r"$d$ (Å)")
        self.ax_strain.set_ylabel(r"$w_i$")
        self.ax_strain.set_xlabel(r"$\varepsilon$ (%)")
        self._install_span()
        self.canvas.draw_idle()

    def draw_spectrum(
        self,
        spectrum: Spectrum,
        xmin: float | None,
        xmax: float | None,
        *,
        original: np.ndarray | None = None,
        baseline: np.ndarray | None = None,
        comb_xmin: float | None = None,
        comb_xmax: float | None = None,
        title: str | None = None,
    ) -> None:
        for ax in (self.ax_main, self.ax_res, self.ax_strain):
            ax.clear()
        self._style_axes()
        if original is not None:
            xo, yo = _thin(spectrum.x, np.asarray(original, dtype=float))
            self.ax_main.plot(xo, yo, color="#bbbbbb", lw=0.9, label="raw")
        x, y = _thin(spectrum.x, spectrum.y)
        self.ax_main.plot(x, y, "o-", color="#111111", ms=2.4, lw=0.7, label="data")
        if baseline is not None:
            xb, yb = _thin(spectrum.x, np.asarray(baseline, dtype=float))
            finite = np.isfinite(yb)
            if np.any(finite):
                self.ax_main.plot(xb[finite], yb[finite], color="#009E73", lw=1.2, ls="--", label="baseline")
        if xmin is not None and xmax is not None and abs(xmax - xmin) < 0.98 * (float(spectrum.x.max()) - float(spectrum.x.min()) + 1e-12):
            self.ax_main.axvspan(xmin, xmax, color="#0b6f71", alpha=0.10)
            self.ax_main.axvline(xmin, color="#0b6f71", lw=0.9, ls="--")
            self.ax_main.axvline(xmax, color="#0b6f71", lw=0.9, ls="--")
        if comb_xmin is not None and comb_xmax is not None:
            self.ax_main.axvline(comb_xmin, color="#E69F00", lw=1.1, ls=":")
            self.ax_main.axvline(comb_xmax, color="#E69F00", lw=1.1, ls=":")
        self.ax_main.set_ylabel("Intensity (counts/s)")
        self.ax_main.legend(loc="upper right", fontsize=8, frameon=False)
        xlabel = AXIS_LABELS.get(spectrum.axis, "x")
        self.ax_res.set_xlabel(xlabel)
        self.ax_res.set_ylabel("Res.")
        self._install_span()
        self.canvas.draw_idle()

    def draw_result(self, result: FitResult, spectrum: Spectrum | None = None, options: ViewOptions | None = None) -> None:
        for ax in (self.ax_main, self.ax_res, self.ax_strain):
            ax.clear()
        self._style_axes()
        view = options or self.view
        self.view = view
        render_split(self.ax_main, self.ax_res if view.show_residual else None, result, view)
        self.ax_res.set_visible(view.show_residual)
        if view.show_strain and getattr(result, "mechanics", None) is not None:
            render_strain(self.ax_strain, result, as_stress=view.strain_as_stress)
            self.ax_strain.set_visible(True)
        else:
            self.ax_strain.set_visible(False)
        self._install_span()
        self.canvas.draw_idle()
