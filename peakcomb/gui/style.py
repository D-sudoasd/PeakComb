"""Qt stylesheet for the PeakComb workbench."""

from __future__ import annotations


def application_stylesheet() -> str:
    return """
    QMainWindow, QWidget {
        background: #f5f7f6;
        color: #1f2937;
        font-family: "Microsoft YaHei", DengXian, Arial, "Segoe UI", sans-serif;
        font-size: 12px;
    }
    QFrame#control_panel {
        background: #ffffff;
        border-right: 1px solid #d8dfdc;
    }
    QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox, QTableWidget {
        background: #ffffff;
        border: 1px solid #d8dfdc;
        border-radius: 6px;
        padding: 3px 6px;
        min-height: 24px;
    }
    QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus {
        border: 2px solid #005fcc;
    }
    QPushButton {
        background: #ecefed;
        border: 1px solid #c9d4cf;
        border-radius: 6px;
        padding: 6px 10px;
        min-height: 28px;
    }
    QPushButton:hover { background: #e3ebe7; }
    QPushButton#primary {
        background: #0b6f71;
        color: #ffffff;
        border: 1px solid #0b6f71;
        font-weight: 600;
    }
    QPushButton#primary:hover { background: #09585a; }
    QGroupBox {
        font-weight: 600;
        border: 1px solid #d8dfdc;
        border-radius: 8px;
        margin-top: 10px;
        padding-top: 8px;
        background: #ffffff;
    }
    QGroupBox::title {
        subcontrol-origin: margin;
        left: 10px;
        padding: 0 4px;
        color: #135f5b;
    }
    QStatusBar { background: #eaf0ed; }
    QLabel#step_label {
        background: #edf7f5;
        border: 1px solid #afd5ce;
        border-radius: 8px;
        padding: 8px;
        color: #135f5b;
        font-weight: 600;
    }
    """
