"""PeakComb desktop application."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from PySide6.QtCore import QObject, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QAction, QCursor
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSlider,
    QSpinBox,
    QSplitter,
    QStatusBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .. import __version__
from ..analysis import half_max_fwhm, slice_window
from ..baseline import BASELINE_METHODS, pybaselines_available, subtract_in_window
from ..export import export_run
from ..fityk_runtime import detect_cfityk
from ..mechanics import D0_AREA_WEIGHTED, D0_MANUAL, D0_PEAK_MAX, analyze_components
from ..io import load_spectrum
from ..limits import MAX_COMB_PEAKS, MAX_TABLE_ROWS
from ..figures import ViewOptions, plot_overlay
from ..models import BackgroundKind, FitMode, FitResult, PeakSeed, SeedMethod, Session, Spectrum
from ..pipeline import official_fit, preview_fit
from ..seed import seed_peaks
from .plot import PlotCanvas
from .style import application_stylesheet

ABOUT_TEXT = (
    "PeakComb 把一条宽 XRD/SXRD 峰拆成一组共享 FWHM 的 PseudoVoigt 小峰。\n\n"
    "两个模式：\n"
    "• 密梳拆分（Distribution）：连续 d 分布，例如 Type II 内应变。"
    "等间距密梳，中心锁死，只拟合非负高度。\n"
    "• 少数分立峰（Domain）：少量 nanodomain / 劈裂亚峰。"
    "N 个等宽峰，高度自由，中心可在小窗口内动。\n\n"
    "共享 FWHM 代表单个 domain 的仪器+尺寸峰宽，不要让每个峰再自由变宽。\n"
    "表里的 w（area_frac）是衍射强度份额，不是自动的体积分数。\n\n"
    "预览用 NNLS；官方拟合用 Fityk（cfityk）。没有 Fityk 时仍可预览并导出。\n"
    "这是独立工具，不修改 PeakTrace。"
)


class _PreviewWorker(QObject):
    finished = Signal(object, int)
    failed = Signal(str, int)

    def __init__(self, session: Session, spectrum: Spectrum, generation: int) -> None:
        super().__init__()
        self._session = session
        self._spectrum = spectrum
        self._generation = generation

    def run(self) -> None:
        try:
            result = preview_fit(self._session, self._spectrum)
            self.finished.emit(result, self._generation)
        except Exception as exc:
            self.failed.emit(str(exc), self._generation)

AXIS_OPTIONS = [
    ("d-spacing, Å", "d_a"),
    ("q, Å⁻¹", "q_inv_a"),
    ("2θ, deg", "two_theta_deg"),
]


def _suggest_window(spectrum: Spectrum, fwhm: float) -> tuple[float, float]:
    peak_at = float(spectrum.x[int(spectrum.y.argmax())])
    half = max(8.0 * float(fwhm), 0.04)
    xmin = max(float(spectrum.x.min()), peak_at - half)
    xmax = min(float(spectrum.x.max()), peak_at + half)
    if xmax <= xmin:
        return float(spectrum.x.min()), float(spectrum.x.max())
    return xmin, xmax


class PeakCombWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"PeakComb {__version__}")
        self.resize(1400, 860)
        self.spectrum: Spectrum | None = None
        self.spectrum_raw: Spectrum | None = None
        self.baseline_y = None
        self.session = Session()
        self.result: FitResult | None = None
        self.click_adds = False
        self._preview_timer = QTimer(self)
        self._preview_timer.setSingleShot(True)
        self._preview_timer.setInterval(280)
        self._preview_timer.timeout.connect(self.refresh_preview)
        self._building = False
        self._preview_busy = False
        self._preview_pending = False
        self._preview_generation = 0
        self._thread: QThread | None = None
        self._worker: _PreviewWorker | None = None
        self._build()
        self._update_cfityk_status()

    def _build(self) -> None:
        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.plot = PlotCanvas(on_range=self._on_plot_range, on_click=self._on_plot_click)
        splitter.addWidget(self._build_controls())
        splitter.addWidget(self.plot)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([320, 1080])
        self.setCentralWidget(splitter)
        self.setStatusBar(QStatusBar())
        open_action = QAction("打开谱", self)
        open_action.triggered.connect(self.open_spectrum)
        export_action = QAction("导出", self)
        export_action.triggered.connect(self.export_current)
        menubar = self.menuBar()
        file_menu = menubar.addMenu("文件")
        file_menu.addAction(open_action)
        file_menu.addAction(export_action)
        help_menu = menubar.addMenu("帮助")
        about_action = QAction("这是什么", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    def _build_controls(self) -> QWidget:
        panel = QFrame()
        panel.setObjectName("control_panel")
        panel.setMinimumWidth(310)
        inner = QWidget()
        layout = QVBoxLayout(inner)
        layout.setContentsMargins(10, 10, 10, 10)

        self.step_label = QLabel("把宽峰拆成一组等宽 PseudoVoigt 小峰。打开谱 → 拖选峰区 → 调 N 和 FWHM。")
        self.step_label.setWordWrap(True)
        self.step_label.setObjectName("step_label")
        layout.addWidget(self.step_label)

        open_btn = QPushButton("打开谱…")
        open_btn.setObjectName("primary")
        open_btn.clicked.connect(self.open_spectrum)
        layout.addWidget(open_btn)
        self.path_label = QLabel("未打开文件")
        self.path_label.setWordWrap(True)
        layout.addWidget(self.path_label)

        self.axis_combo = QComboBox()
        for label, code in AXIS_OPTIONS:
            self.axis_combo.addItem(label, code)
        self.mode_combo = QComboBox()
        self.mode_combo.addItem("密梳拆分", FitMode.DISTRIBUTION.value)
        self.mode_combo.addItem("少数分立峰", FitMode.DOMAIN.value)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        form = QFormLayout()
        form.addRow("横轴", self.axis_combo)
        form.addRow("模式", self.mode_combo)
        layout.addLayout(form)

        self.n_label = QLabel("N  9")
        self.n_slider = QSlider(Qt.Orientation.Horizontal)
        self.n_slider.setRange(3, 16)
        self.n_slider.setValue(9)
        self.n_slider.valueChanged.connect(self._on_n_slider_moved)
        self.n_slider.sliderReleased.connect(self._on_n_slider_released)
        layout.addWidget(self.n_label)
        layout.addWidget(self.n_slider)

        self.fwhm_label = QLabel("FWHM  0.010")
        self.fwhm_slider = QSlider(Qt.Orientation.Horizontal)
        self.fwhm_slider.setRange(40, 200)
        self.fwhm_slider.setValue(100)
        self.fwhm_slider.valueChanged.connect(self._on_fwhm_slider_moved)
        self.fwhm_slider.sliderReleased.connect(self._on_fwhm_slider_released)
        layout.addWidget(self.fwhm_label)
        layout.addWidget(self.fwhm_slider)

        self.alpha_label = QLabel("填色  40%")
        self.alpha_slider = QSlider(Qt.Orientation.Horizontal)
        self.alpha_slider.setRange(10, 80)
        self.alpha_slider.setValue(40)
        self.alpha_slider.valueChanged.connect(self._redraw_view)
        layout.addWidget(self.alpha_label)
        layout.addWidget(self.alpha_slider)

        self.show_fill = QCheckBox("彩色填色")
        self.show_fill.setChecked(True)
        self.show_fill.toggled.connect(self._redraw_view)
        self.show_bkg = QCheckBox("背景线")
        self.show_bkg.setChecked(True)
        self.show_bkg.toggled.connect(self._redraw_view)
        self.show_res = QCheckBox("残差")
        self.show_res.setChecked(True)
        self.show_res.toggled.connect(self._redraw_view)
        self.data_points = QCheckBox("数据用点")
        self.data_points.setChecked(True)
        self.data_points.toggled.connect(self._redraw_view)
        for box in (self.show_fill, self.show_bkg, self.show_res, self.data_points):
            layout.addWidget(box)

        tools = QHBoxLayout()
        self.tool_select = QRadioButton("拖选")
        self.tool_add = QRadioButton("加点")
        self.tool_browse = QRadioButton("浏览")
        self.tool_select.setChecked(True)
        group = QButtonGroup(self)
        for button in (self.tool_select, self.tool_add, self.tool_browse):
            group.addButton(button)
            tools.addWidget(button)
        self.tool_select.toggled.connect(lambda on: on and self.plot.set_tool("select"))
        self.tool_add.toggled.connect(lambda on: on and self.plot.set_tool("add"))
        self.tool_browse.toggled.connect(lambda on: on and self.plot.set_tool("browse"))
        layout.addLayout(tools)

        self.auto_baseline = QCheckBox("拖选后自动扣背底")
        self.auto_baseline.setChecked(False)
        layout.addWidget(self.auto_baseline)
        self.baseline_status = QLabel("未扣背底")
        self.baseline_status.setWordWrap(True)
        layout.addWidget(self.baseline_status)

        export_btn = QPushButton("导出当前图")
        export_btn.setObjectName("primary")
        export_btn.clicked.connect(self.export_current)
        layout.addWidget(export_btn)
        fit_btn = QPushButton("Fityk 拟合")
        fit_btn.clicked.connect(self.run_fityk)
        layout.addWidget(fit_btn)
        self.cfityk_label = QLabel()
        self.cfityk_label.setWordWrap(True)
        layout.addWidget(self.cfityk_label)

        self.xmin_spin = self._float_spin(-50.0, 400.0, 2.20, 6)
        self.xmax_spin = self._float_spin(-50.0, 400.0, 2.40, 6)
        self.n_spin = QSpinBox()
        self.n_spin.setRange(1, MAX_COMB_PEAKS)
        self.n_spin.setValue(9)
        self.n_spin.editingFinished.connect(self._on_n_changed)
        self.fwhm_spin = self._float_spin(1e-5, 2.0, 0.010, 5)
        range_row = QHBoxLayout()
        range_row.addWidget(QLabel("范围"))
        range_row.addWidget(self.xmin_spin)
        range_row.addWidget(self.xmax_spin)
        layout.addLayout(range_row)

        mech = QGroupBox("内应变 / 静水内应力")
        mech_form = QFormLayout(mech)
        self.d0_source = QComboBox()
        self.d0_source.addItem("峰顶分量", D0_PEAK_MAX)
        self.d0_source.addItem("面积加权平均", D0_AREA_WEIGHTED)
        self.d0_source.addItem("手动输入", D0_MANUAL)
        self.d0_source.currentIndexChanged.connect(self._update_mechanics)
        mech_form.addRow("d0 来源", self.d0_source)
        self.d0_spin = QDoubleSpinBox()
        self.d0_spin.setRange(0.0, 20.0)
        self.d0_spin.setDecimals(5)
        self.d0_spin.setValue(0.0)
        self.d0_spin.editingFinished.connect(self._update_mechanics)
        mech_form.addRow("d0 / Å", self.d0_spin)
        self.e_spin = QDoubleSpinBox()
        self.e_spin.setRange(1.0, 400.0)
        self.e_spin.setDecimals(1)
        self.e_spin.setValue(80.0)
        self.e_spin.editingFinished.connect(self._update_mechanics)
        mech_form.addRow("E / GPa", self.e_spin)
        self.nu_spin = QDoubleSpinBox()
        self.nu_spin.setRange(0.0, 0.49)
        self.nu_spin.setDecimals(3)
        self.nu_spin.setValue(0.33)
        self.nu_spin.editingFinished.connect(self._update_mechanics)
        mech_form.addRow("ν", self.nu_spin)
        self.compute_stress = QCheckBox("换算 σ_h = E/(1−2ν) ε")
        self.compute_stress.setChecked(True)
        self.compute_stress.toggled.connect(self._update_mechanics)
        mech_form.addRow(self.compute_stress)
        self.show_strain = QCheckBox("显示 w–ε 图")
        self.show_strain.setChecked(True)
        self.show_strain.toggled.connect(self._redraw_view)
        self.strain_as_stress = QCheckBox("横轴用 σ_h")
        self.strain_as_stress.setChecked(False)
        self.strain_as_stress.toggled.connect(self._redraw_view)
        mech_form.addRow(self.show_strain)
        mech_form.addRow(self.strain_as_stress)
        self.mech_label = QLabel("拟合后给出 ε̄ 与平衡检查。")
        self.mech_label.setWordWrap(True)
        mech_form.addRow(self.mech_label)
        layout.addWidget(mech)

        adv_btn = QPushButton("高级参数")
        adv_btn.setCheckable(True)
        layout.addWidget(adv_btn)
        advanced = QWidget()
        advanced.setVisible(False)
        adv_btn.toggled.connect(advanced.setVisible)
        adv_form = QFormLayout(advanced)
        self.energy_spin = QDoubleSpinBox()
        self.energy_spin.setRange(0.1, 200.0)
        self.energy_spin.setDecimals(4)
        self.energy_spin.setValue(83.0)
        adv_form.addRow("能量 / keV", self.energy_spin)
        self.baseline_combo = QComboBox()
        for name in BASELINE_METHODS:
            self.baseline_combo.addItem(name, name)
        self.baseline_combo.setCurrentText("asls")
        adv_form.addRow("背底方法", self.baseline_combo)
        apply_base = QPushButton("手动扣背底")
        apply_base.clicked.connect(lambda: self.apply_baseline(preview=True))
        adv_form.addRow(apply_base)
        self.shape_spin = self._float_spin(0.0, 1.0, 0.5, 3)
        adv_form.addRow("shape", self.shape_spin)
        self.refine_fwhm = QCheckBox("Fityk 放开共享 FWHM")
        self.refine_shape = QCheckBox("Fityk 放开共享 shape")
        adv_form.addRow(self.refine_fwhm)
        adv_form.addRow(self.refine_shape)
        self.bkg_combo = QComboBox()
        for kind in BackgroundKind:
            self.bkg_combo.addItem(kind.value, kind.value)
        self.bkg_combo.setCurrentText("linear")
        adv_form.addRow("残差背景", self.bkg_combo)
        self.spacing_spin = self._float_spin(0.1, 2.0, 0.4, 3)
        adv_form.addRow("密梳间距/FWHM", self.spacing_spin)
        self.seed_combo = QComboBox()
        self.seed_combo.addItem("对称密梳", SeedMethod.SYMMETRIC.value)
        self.seed_combo.addItem("等间距", SeedMethod.EQUAL.value)
        self.seed_combo.addItem("局部极大", SeedMethod.MAXIMA.value)
        self.seed_combo.addItem("手动", SeedMethod.MANUAL.value)
        adv_form.addRow("铺峰", self.seed_combo)
        self.comb_same = QCheckBox("密梳范围 = 拟合范围")
        self.comb_same.setChecked(True)
        self.comb_same.toggled.connect(self._sync_comb_enabled)
        self.comb_xmin_spin = self._float_spin(-50.0, 400.0, 2.20, 6)
        self.comb_xmax_spin = self._float_spin(-50.0, 400.0, 2.40, 6)
        adv_form.addRow(self.comb_same)
        adv_form.addRow("密梳 xmin", self.comb_xmin_spin)
        adv_form.addRow("密梳 xmax", self.comb_xmax_spin)
        self.click_box = QCheckBox("加点模式")
        self.click_box.toggled.connect(self._on_click_box)
        adv_form.addRow(self.click_box)
        measure_btn = QPushButton("用当前范围估 FWHM")
        measure_btn.clicked.connect(self._measure_fwhm)
        adv_form.addRow(measure_btn)
        full_btn = QPushButton("范围铺满全谱")
        full_btn.clicked.connect(self._window_full_spectrum)
        adv_form.addRow(full_btn)
        layout.addWidget(advanced)
        self.xmin_spin.editingFinished.connect(self._on_fit_range_changed)
        self.xmax_spin.editingFinished.connect(self._on_fit_range_changed)
        self.spacing_spin.editingFinished.connect(self._on_spacing_changed)
        self._sync_comb_enabled()

        table_box = QGroupBox("分量")
        table_layout = QVBoxLayout(table_box)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["d", "w", "ε %", "σ_h MPa", "FWHM"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        del_btn = QPushButton("删除选中峰")
        del_btn.clicked.connect(self.delete_selected)
        table_layout.addWidget(self.table)
        table_layout.addWidget(del_btn)
        layout.addWidget(table_box, 1)

        layout.addStretch(0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(inner)
        wrap = QVBoxLayout(panel)
        wrap.setContentsMargins(0, 0, 0, 0)
        wrap.addWidget(scroll)
        return panel

    def _float_spin(self, lo: float, hi: float, value: float, decimals: int) -> QDoubleSpinBox:
        box = QDoubleSpinBox()
        box.setRange(lo, hi)
        box.setDecimals(decimals)
        box.setValue(value)
        box.editingFinished.connect(self._schedule_preview)
        return box

    def _show_about(self) -> None:
        QMessageBox.information(self, f"PeakComb {__version__}", ABOUT_TEXT)

    def _update_cfityk_status(self) -> None:
        found = detect_cfityk(self.session.cfityk_path or None)
        if found:
            self.cfityk_label.setText(f"cfityk: {found}")
        else:
            self.cfityk_label.setText("cfityk 未找到。仍可用 NNLS 预览与导出。")

    def open_spectrum(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "打开 XRD 谱",
            "",
            "Spectra (*.xy *.xye *.chi *.csv *.txt *.dat);;All files (*.*)",
        )
        if not path:
            return
        try:
            spectrum = load_spectrum(
                path,
                axis=self.axis_combo.currentData(),
                energy_kev=self.energy_spin.value(),
            )
        except Exception as exc:
            QMessageBox.warning(self, "读取失败", str(exc))
            return
        self.spectrum_raw = spectrum
        self.spectrum = Spectrum(
            x=spectrum.x.copy(),
            y=spectrum.y.copy(),
            sigma=None if spectrum.sigma is None else spectrum.sigma.copy(),
            path=spectrum.path,
            axis=spectrum.axis,
            wavelength_a=spectrum.wavelength_a,
            energy_kev=spectrum.energy_kev,
        )
        self.baseline_y = None
        self.session.baseline_method = "none"
        self.baseline_status.setText("未扣背底")
        self.path_label.setText(str(Path(path).name))
        self.session.input_path = path
        self.session.peaks = []
        self.result = None
        lo, hi = float(spectrum.x.min()), float(spectrum.x.max())
        pad = max(0.05 * (hi - lo), 1e-6)
        for box in (self.xmin_spin, self.xmax_spin, self.comb_xmin_spin, self.comb_xmax_spin):
            box.setRange(lo - pad, hi + pad)
        lo, hi = float(spectrum.x.min()), float(spectrum.x.max())
        self._set_range_spins(lo, hi)
        self._redraw_spectrum()
        self.step_label.setText("② 在图上按住左键，拖过要拆的那一团峰")
        self.statusBar().showMessage(f"已加载 {len(spectrum.x)} 点。拖选峰区即可拆峰。")

    def _collect_session(self) -> Session | None:
        if self.spectrum is None:
            return None
        d0 = self.d0_spin.value()
        self.session.d0_source = self.d0_source.currentData()
        self.session.e_gpa = self.e_spin.value()
        self.session.nu = self.nu_spin.value()
        self.session.compute_stress = self.compute_stress.isChecked()
        self.session.axis = self.axis_combo.currentData()
        self.session.energy_kev = self.energy_spin.value()
        self.session.wavelength_a = self.spectrum.wavelength_a
        self.session.xmin = self.xmin_spin.value()
        self.session.xmax = self.xmax_spin.value()
        if self.comb_same.isChecked():
            self.session.comb_xmin = None
            self.session.comb_xmax = None
        else:
            self.session.comb_xmin = self.comb_xmin_spin.value()
            self.session.comb_xmax = self.comb_xmax_spin.value()
        self.session.baseline_method = self.baseline_combo.currentData()
        self.session.mode = FitMode(self.mode_combo.currentData())
        self.session.fwhm = self.fwhm_spin.value()
        self.session.shape = self.shape_spin.value()
        self.session.refine_shared_fwhm = self.refine_fwhm.isChecked()
        self.session.refine_shared_shape = self.refine_shape.isChecked()
        self.session.background = BackgroundKind(self.bkg_combo.currentData())
        self.session.n_peaks = self.n_spin.value()
        self.session.spacing_factor = self.spacing_spin.value()
        self.session.seed_method = SeedMethod(self.seed_combo.currentData())
        self.session.d0 = None if d0 <= 0.0 else d0
        if self.session.d0_source == D0_MANUAL and self.session.d0 is None:
            self.session.d0_source = D0_PEAK_MAX
        self.spectrum.axis = self.session.axis
        return self.session

    def reseed(self) -> None:
        if self.spectrum is None:
            return
        session = self._collect_session()
        if session is None:
            return
        if session.mode is FitMode.DISTRIBUTION and session.seed_method not in {
            SeedMethod.EQUAL,
            SeedMethod.GRID,
            SeedMethod.SYMMETRIC,
        }:
            session.seed_method = SeedMethod.SYMMETRIC
        try:
            seed_peaks(session, self.spectrum)
        except Exception as exc:
            QMessageBox.warning(self, "铺峰失败", str(exc))
            return
        self.n_spin.blockSignals(True)
        self.n_spin.setValue(len(session.peaks))
        self.n_spin.blockSignals(False)
        self.refresh_preview()

    def _on_mode_changed(self) -> None:
        distribution = self.mode_combo.currentData() == FitMode.DISTRIBUTION.value
        self.click_box.setEnabled(not distribution)
        if distribution:
            self.click_box.setChecked(False)
            self.session.seed_method = SeedMethod.SYMMETRIC
            self.seed_combo.setCurrentText("对称密梳")
        self.session.peaks = []
        if self.spectrum is None:
            return
        lo, hi = float(self.spectrum.x.min()), float(self.spectrum.x.max())
        xmin, xmax = self.xmin_spin.value(), self.xmax_spin.value()
        if (xmax - xmin) > 0.95 * (hi - lo + 1e-12):
            self._redraw_spectrum()
            return
        self.reseed()

    def _on_fit_range_changed(self) -> None:
        if self.comb_same.isChecked():
            self._sync_comb_enabled()
        elif self.spectrum is not None:
            self._redraw_spectrum()

    def _on_spacing_changed(self) -> None:
        if self.mode_combo.currentData() != FitMode.DISTRIBUTION.value:
            return
        xmin = self.comb_xmin_spin.value() if not self.comb_same.isChecked() else self.xmin_spin.value()
        xmax = self.comb_xmax_spin.value() if not self.comb_same.isChecked() else self.xmax_spin.value()
        step = self.spacing_spin.value() * self.fwhm_spin.value()
        if step <= 0.0 or xmax <= xmin:
            return
        n_peaks = int(np.floor((xmax - xmin) / step)) + 1
        n_peaks = min(max(n_peaks, 2), MAX_COMB_PEAKS)
        self.n_spin.blockSignals(True)
        self.n_spin.setValue(n_peaks)
        self.n_spin.blockSignals(False)
        self.session.seed_method = SeedMethod.SYMMETRIC

    def _view_options(self) -> ViewOptions:
        return ViewOptions(
            show_fill=self.show_fill.isChecked(),
            show_outlines=True,
            show_bkg=self.show_bkg.isChecked(),
            show_residual=self.show_res.isChecked(),
            data_as_points=self.data_points.isChecked(),
            fill_alpha=self.alpha_slider.value() / 100.0,
            show_strain=self.show_strain.isChecked(),
            strain_as_stress=self.strain_as_stress.isChecked(),
        )

    def _update_mechanics(self, *_args) -> None:
        if self.result is None:
            return
        session = self._collect_session()
        if session is None:
            return
        try:
            report = analyze_components(self.result.components, session)
        except Exception as exc:
            self.mech_label.setText(str(exc))
            return
        self.result.mechanics = report
        self.result.session = session
        if session.d0_source != D0_MANUAL:
            self.d0_spin.blockSignals(True)
            self.d0_spin.setValue(report.d0)
            self.d0_spin.blockSignals(False)
        mean_stress = "—" if not report.compute_stress else f"{report.mean_stress_mpa:.0f} MPa"
        self.mech_label.setText(
            f"d0={report.d0:.5f} Å\n"
            f"ε̄={100.0 * report.mean_strain:.3f}%    σ̄_h={mean_stress}\n"
            f"拉 {100.0 * report.tension_weight:.0f}% / 压 {100.0 * report.compression_weight:.0f}%\n"
            f"{report.note}"
        )
        self._fill_table(self.result)
        self.plot.draw_result(self.result, options=self._view_options())

    def _redraw_view(self, *_args) -> None:
        self.alpha_label.setText(f"填色  {self.alpha_slider.value()}%")
        if self.result is not None:
            self.plot.draw_result(self.result, options=self._view_options())

    def _on_n_slider_moved(self, value: int) -> None:
        self.n_label.setText(f"N  {value}")
        self.n_spin.blockSignals(True)
        self.n_spin.setValue(value)
        self.n_spin.blockSignals(False)

    def _on_n_slider_released(self) -> None:
        self.session.n_peaks = self.n_slider.value()
        if self.mode_combo.currentData() == FitMode.DISTRIBUTION.value:
            self.session.seed_method = SeedMethod.SYMMETRIC
            self.seed_combo.setCurrentText("对称密梳")
        else:
            self.session.seed_method = SeedMethod.EQUAL
        if self.spectrum is not None:
            self.reseed()

    def _on_fwhm_slider_moved(self, value: int) -> None:
        fwhm = value / 10000.0
        self.fwhm_label.setText(f"FWHM  {fwhm:.3f}")
        self.fwhm_spin.blockSignals(True)
        self.fwhm_spin.setValue(fwhm)
        self.fwhm_spin.blockSignals(False)

    def _on_fwhm_slider_released(self) -> None:
        self.session.fwhm = self.fwhm_slider.value() / 10000.0
        if self.spectrum is not None:
            self.reseed()

    def _on_n_changed(self) -> None:
        if self.spectrum is None:
            return
        session = self._collect_session()
        if session is None:
            return
        if self.n_spin.value() == len(self.session.peaks) and session.mode is FitMode.DOMAIN:
            return
        self.session.peaks = []
        self.session.n_peaks = self.n_spin.value()
        if session.mode is FitMode.DISTRIBUTION:
            self.session.seed_method = SeedMethod.SYMMETRIC
            self.seed_combo.setCurrentText("对称密梳")
        self.reseed()

    def _schedule_preview(self) -> None:
        if self.spectrum is None:
            return
        self._preview_timer.start()

    def refresh_preview(self) -> None:
        if self.spectrum is None:
            return
        session = self._collect_session()
        if session is None:
            return
        if session.mode is FitMode.DISTRIBUTION:
            try:
                seed_peaks(session, self.spectrum)
            except Exception as exc:
                self.statusBar().showMessage(str(exc), 8000)
                return
            self.n_spin.blockSignals(True)
            self.n_spin.setValue(len(session.peaks))
            self.n_spin.blockSignals(False)
            self.n_slider.blockSignals(True)
            self.n_slider.setValue(min(max(len(session.peaks), 3), 16))
            self.n_slider.blockSignals(False)
            self.n_label.setText(f"N  {len(session.peaks)}")
            if len(session.peaks) >= 2:
                span = session.peaks[-1].center - session.peaks[0].center
                step = span / (len(session.peaks) - 1)
                if session.fwhm > 0.0:
                    self.spacing_spin.blockSignals(True)
                    self.spacing_spin.setValue(max(step / session.fwhm, 0.1))
                    self.spacing_spin.blockSignals(False)
        if self._preview_busy:
            self._preview_pending = True
            return
        self._start_preview_job(session)

    def _start_preview_job(self, session: Session) -> None:
        self._preview_busy = True
        self._preview_pending = False
        self._preview_generation += 1
        generation = self._preview_generation
        snapshot = Session.from_dict(session.to_dict())
        self.statusBar().showMessage("正在预览…")
        QApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
        thread = QThread(self)
        worker = _PreviewWorker(snapshot, self.spectrum, generation)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._on_preview_finished)
        worker.failed.connect(self._on_preview_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        self._thread = thread
        self._worker = worker
        thread.start()

    def _finish_preview_job(self) -> None:
        QApplication.restoreOverrideCursor()
        self._preview_busy = False
        self._worker = None
        self._thread = None
        if self._preview_pending:
            self._preview_pending = False
            QTimer.singleShot(0, self.refresh_preview)

    def _on_preview_finished(self, result: FitResult, generation: int) -> None:
        if generation != self._preview_generation:
            self._finish_preview_job()
            return
        self.result = result
        self._update_mechanics()
        extra = ""
        if result.qc.n_peaks >= MAX_COMB_PEAKS:
            extra = f"  （密梳已封顶 {MAX_COMB_PEAKS} 峰）"
        self.statusBar().showMessage(
            f"NNLS  Rwp={result.qc.rwp:.4g}  N={result.qc.n_peaks}  FWHM={result.session.fwhm:g}{extra}"
        )
        self._finish_preview_job()

    def _on_preview_failed(self, message: str, generation: int) -> None:
        if generation == self._preview_generation:
            self.statusBar().showMessage(message, 8000)
        self._finish_preview_job()

    def run_fityk(self) -> None:
        if self.spectrum is None:
            QMessageBox.information(self, "PeakComb", "先打开一条谱。")
            return
        session = self._collect_session()
        if session is None:
            return
        work = Path(self.session.input_path).expanduser().resolve().parent / "_peakcomb_work"
        work.mkdir(parents=True, exist_ok=True)
        try:
            result = official_fit(session, self.spectrum, work)
        except Exception as exc:
            QMessageBox.warning(self, "拟合失败", str(exc))
            return
        self.result = result
        self._update_mechanics()
        self.statusBar().showMessage(
            f"{result.fitter}  Rwp={result.qc.rwp:.4g}  shared_FWHM={'OK' if result.qc.shared_fwhm_ok else 'CHECK'}"
        )
        if result.qc.warnings:
            QMessageBox.information(self, "QC", "\n".join(result.qc.warnings))

    def export_current(self) -> None:
        if self.result is None:
            QMessageBox.information(self, "PeakComb", "先预览或拟合。")
            return
        source = Path(self.session.input_path)
        default = source.parent / f"{source.stem}_peakcomb"
        out_dir = QFileDialog.getExistingDirectory(self, "选择导出目录", str(default))
        if not out_dir:
            out_dir = str(default)
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        bundle = export_run(self.result, out_dir, self.session)
        pretty = Path(out_dir) / "overlay_current.png"
        plot_overlay(self.result, pretty, options=self._view_options(), dpi=360)
        QMessageBox.information(self, "已导出", f"当前图：{pretty}\n其余文件：{bundle.output_dir}")

    def _fill_table(self, result: FitResult) -> None:
        rows = result.components[:MAX_TABLE_ROWS]
        self.table.setRowCount(len(rows))
        self.table.setUpdatesEnabled(False)
        try:
            for row, item in enumerate(rows):
                values = [
                    f"{item.d_a:.5f}" if item.d_a is not None else f"{item.center:.5f}",
                    f"{item.area_frac:.3f}",
                    "" if item.strain is None else f"{100.0 * item.strain:.3f}",
                    "" if item.sigma_h_mpa is None else f"{item.sigma_h_mpa:.0f}",
                    f"{item.fwhm:.4g}",
                ]
                for column, value in enumerate(values):
                    self.table.setItem(row, column, QTableWidgetItem(value))
        finally:
            self.table.setUpdatesEnabled(True)

    def delete_selected(self) -> None:
        rows = sorted({index.row() for index in self.table.selectedIndexes()}, reverse=True)
        if not rows or not self.session.peaks:
            return
        for row in rows:
            if 0 <= row < len(self.session.peaks):
                del self.session.peaks[row]
        self.session.n_peaks = len(self.session.peaks)
        self.session.seed_method = SeedMethod.MANUAL
        self.n_spin.blockSignals(True)
        self.n_spin.setValue(max(self.session.n_peaks, 1))
        self.n_spin.blockSignals(False)
        self.refresh_preview()

    def _on_click_box(self, flag: bool) -> None:
        self.click_adds = bool(flag)
        if flag:
            self.tool_add.setChecked(True)

    def _on_plot_range(self, xmin: float, xmax: float) -> None:
        if self.spectrum is None:
            return
        self._set_range_spins(xmin, xmax)
        if self.auto_baseline.isChecked() and self.baseline_combo.currentData() != "none":
            self.apply_baseline(preview=False)
        else:
            self._redraw_spectrum()
        pad = 0.12 * max(xmax - xmin, 1e-12)
        self.plot.ax_main.set_xlim(xmin - pad, xmax + pad)
        self.session.peaks = []
        self.session.seed_method = SeedMethod.SYMMETRIC
        self.seed_combo.setCurrentText("对称密梳")
        self.step_label.setText("已按峰顶左右对称铺峰。拖 N / FWHM 可改。")
        self.reseed()

    def _on_plot_click(self, x: float, y: float) -> None:
        if self.spectrum is None or not self.tool_add.isChecked():
            return
        if self.mode_combo.currentData() == FitMode.DISTRIBUTION.value:
            return
        self.session.peaks.append(
            PeakSeed(center=x, height=max(y, 1.0), label=f"P{len(self.session.peaks) + 1}")
        )
        self.session.seed_method = SeedMethod.MANUAL
        self.session.n_peaks = len(self.session.peaks)
        self.n_spin.blockSignals(True)
        self.n_spin.setValue(self.session.n_peaks)
        self.n_spin.blockSignals(False)
        self.refresh_preview()

    def _set_range_spins(self, xmin: float, xmax: float) -> None:
        self.xmin_spin.blockSignals(True)
        self.xmax_spin.blockSignals(True)
        self.comb_xmin_spin.blockSignals(True)
        self.comb_xmax_spin.blockSignals(True)
        self.xmin_spin.setValue(xmin)
        self.xmax_spin.setValue(xmax)
        if self.comb_same.isChecked():
            self.comb_xmin_spin.setValue(xmin)
            self.comb_xmax_spin.setValue(xmax)
        self.xmin_spin.blockSignals(False)
        self.xmax_spin.blockSignals(False)
        self.comb_xmin_spin.blockSignals(False)
        self.comb_xmax_spin.blockSignals(False)
        self._sync_comb_enabled()

    def _sync_comb_enabled(self, *_args) -> None:
        same = self.comb_same.isChecked()
        self.comb_xmin_spin.setEnabled(not same)
        self.comb_xmax_spin.setEnabled(not same)
        if same:
            self.comb_xmin_spin.blockSignals(True)
            self.comb_xmax_spin.blockSignals(True)
            self.comb_xmin_spin.setValue(self.xmin_spin.value())
            self.comb_xmax_spin.setValue(self.xmax_spin.value())
            self.comb_xmin_spin.blockSignals(False)
            self.comb_xmax_spin.blockSignals(False)
        if self.spectrum is not None:
            self._redraw_spectrum()

    def _window_from_view(self) -> None:
        xmin, xmax = self.plot.ax_main.get_xlim()
        self._set_range_spins(xmin, xmax)
        self._redraw_spectrum()
        self._schedule_preview()

    def _window_full_spectrum(self) -> None:
        if self.spectrum is None:
            return
        self._set_range_spins(float(self.spectrum.x.min()), float(self.spectrum.x.max()))
        self._redraw_spectrum()
        self._schedule_preview()

    def _redraw_spectrum(self) -> None:
        if self.spectrum is None:
            return
        original = None if self.spectrum_raw is None else self.spectrum_raw.y
        comb_xmin = None if self.comb_same.isChecked() else self.comb_xmin_spin.value()
        comb_xmax = None if self.comb_same.isChecked() else self.comb_xmax_spin.value()
        self.plot.draw_spectrum(
            self.spectrum,
            self.xmin_spin.value(),
            self.xmax_spin.value(),
            original=original if self.baseline_y is not None else None,
            baseline=self.baseline_y,
            comb_xmin=comb_xmin,
            comb_xmax=comb_xmax,
        )

    def apply_baseline(self, preview: bool = True) -> None:
        if self.spectrum_raw is None:
            QMessageBox.information(self, "PeakComb", "先打开一条谱。")
            return
        method = self.baseline_combo.currentData()
        if method == "none":
            self.restore_spectrum()
            return
        xmin, xmax = self.xmin_spin.value(), self.xmax_spin.value()
        try:
            corrected, baseline_full, fitted = subtract_in_window(
                self.spectrum_raw.x,
                self.spectrum_raw.y,
                xmin,
                xmax,
                method,
            )
        except Exception as exc:
            QMessageBox.warning(self, "扣背底失败", str(exc))
            return
        self.spectrum = Spectrum(
            x=self.spectrum_raw.x.copy(),
            y=corrected,
            sigma=self.spectrum_raw.sigma,
            path=self.spectrum_raw.path,
            axis=self.spectrum_raw.axis,
            wavelength_a=self.spectrum_raw.wavelength_a,
            energy_kev=self.spectrum_raw.energy_kev,
        )
        self.baseline_y = baseline_full
        self.session.baseline_method = fitted.method
        self.baseline_status.setText(f"已扣 {fitted.method}（{fitted.backend}）")
        self._redraw_spectrum()
        self.statusBar().showMessage(f"已在 [{xmin:g}, {xmax:g}] 扣背底：{fitted.method}")
        if preview:
            self._schedule_preview()

    def restore_spectrum(self) -> None:
        if self.spectrum_raw is None:
            return
        self.spectrum = Spectrum(
            x=self.spectrum_raw.x.copy(),
            y=self.spectrum_raw.y.copy(),
            sigma=self.spectrum_raw.sigma,
            path=self.spectrum_raw.path,
            axis=self.spectrum_raw.axis,
            wavelength_a=self.spectrum_raw.wavelength_a,
            energy_kev=self.spectrum_raw.energy_kev,
        )
        self.baseline_y = None
        self.session.baseline_method = "none"
        self.baseline_status.setText("未扣背底")
        self._redraw_spectrum()
        self.statusBar().showMessage("已恢复原始谱")

    def _measure_fwhm(self) -> None:
        if self.spectrum is None:
            return
        try:
            x, y, _ = slice_window(self.spectrum, self.xmin_spin.value(), self.xmax_spin.value())
            fwhm = half_max_fwhm(x, y)
        except Exception as exc:
            QMessageBox.warning(self, "FWHM", str(exc))
            return
        self.fwhm_spin.blockSignals(True)
        self.fwhm_spin.setValue(fwhm)
        self.fwhm_spin.blockSignals(False)
        self.statusBar().showMessage(f"窗口半高宽 ≈ {fwhm:.6g}")
        self._schedule_preview()


    def closeEvent(self, event) -> None:
        if self._thread is not None and self._thread.isRunning():
            self._thread.quit()
            self._thread.wait(1000)
        super().closeEvent(event)


def main() -> None:
    app = QApplication.instance() or QApplication(sys.argv)
    app.setStyleSheet(application_stylesheet())
    window = PeakCombWindow()
    window.show()
    sys.exit(app.exec())
