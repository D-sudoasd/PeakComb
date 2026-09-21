"""Locate cfityk, run a Lua script, and parse Fityk .peaks output."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .models import Component, PeakSeed, Session
from .analysis import components_from_arrays

PEAKCOMB_CFITYK_ENV = "PEAKCOMB_CFITYK"
PEAKTRACE_CFITYK_ENV = "FITYK_FLOW_CFITYK"
RAW_ID_RE = re.compile(r"^%_(\S+)$")
PEAKTRACE_ROOT = Path(r"E:\Vibe_coding\PeakTrace")


class FitykRuntimeError(RuntimeError):
    """Raised when cfityk cannot be run or its output cannot be parsed."""


@dataclass
class ParsedPeak:
    raw_id: str
    peak_type: str
    center: float
    height: float
    area: float
    fwhm: float
    parameters: list[float]


def _project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _candidate_paths(explicit: str | None = None) -> list[Path]:
    roots = [_project_root()]
    if PEAKTRACE_ROOT.exists():
        roots.append(PEAKTRACE_ROOT)
    bundled: list[Path] = []
    for root in roots:
        bundled.extend(
            [
                root / "Fityk" / "cfityk.exe",
                root / "Fityk" / "bin" / "cfityk.exe",
            ]
        )
    env_paths = [
        explicit,
        os.environ.get(PEAKCOMB_CFITYK_ENV),
        os.environ.get(PEAKTRACE_CFITYK_ENV),
        shutil.which("cfityk"),
        shutil.which("cfityk.exe"),
    ]
    system: list[Path] = []
    for env_name in ("ProgramFiles", "ProgramFiles(x86)"):
        root = os.environ.get(env_name)
        if root:
            system.append(Path(root) / "Fityk" / "cfityk.exe")
    paths: list[Path] = []
    for item in [*env_paths, *bundled, *system]:
        if not item:
            continue
        paths.append(Path(item))
    return paths


def detect_cfityk(explicit: str | None = None) -> Path | None:
    for candidate in _candidate_paths(explicit):
        if candidate.is_file():
            return candidate
        resolved = shutil.which(str(candidate))
        if resolved:
            return Path(resolved)
    return None


def parse_peaks_text(text: str) -> list[ParsedPeak]:
    peaks: list[ParsedPeak] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 6:
            continue
        if parts[1] not in {"PseudoVoigt", "PseudoVoigtA"}:
            continue
        match = RAW_ID_RE.match(parts[0])
        raw_id = match.group(1) if match else parts[0].lstrip("%_")
        try:
            parameters = [float(value) for value in parts[6:]]
            peaks.append(
                ParsedPeak(
                    raw_id=raw_id,
                    peak_type=parts[1],
                    center=float(parts[2]),
                    height=float(parts[3]),
                    area=float(parts[4]),
                    fwhm=float(parts[5]),
                    parameters=parameters,
                )
            )
        except ValueError:
            continue
    return peaks


def components_from_peaks(
    session: Session,
    parsed: list[ParsedPeak],
    seeds: list[PeakSeed],
) -> list[Component]:
    if not parsed:
        raise FitykRuntimeError("Fityk wrote no PseudoVoigt peaks.")
    ordered = sorted(parsed, key=lambda item: item.center)
    heights = [item.height for item in ordered]
    used_seeds = [
        PeakSeed(center=item.center, height=item.height, label=f"P{index + 1}")
        for index, item in enumerate(ordered)
    ]
    # Prefer Fityk's reported FWHM; shared-width QC compares these values.
    fwhm_values = [item.fwhm for item in ordered if item.fwhm > 0.0]
    fwhm = float(sum(fwhm_values) / len(fwhm_values)) if fwhm_values else session.fwhm
    shape = session.shape
    if ordered[0].parameters and len(ordered[0].parameters) >= 4:
        shape = float(ordered[0].parameters[3])
    result = components_from_arrays(session, used_seeds, heights, fwhm=fwhm, shape=shape)
    # Keep original seed count labels when centers stayed close.
    if len(seeds) == len(result):
        seed_order = sorted(range(len(seeds)), key=lambda i: seeds[i].center)
        for component, seed_index in zip(result, seed_order, strict=True):
            component.label = seeds[seed_index].label or component.label
    return result


def run_cfityk(
    script_path: Path,
    *,
    cfityk_path: str | Path | None = None,
    timeout: float = 120.0,
) -> subprocess.CompletedProcess[str]:
    executable = detect_cfityk(str(cfityk_path) if cfityk_path else None)
    if executable is None:
        raise FitykRuntimeError(
            "cfityk.exe was not found. Set PEAKCOMB_CFITYK or install Fityk."
        )
    command = [str(executable), "-q", "-I", str(script_path)]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=str(script_path.parent),
    )
    if completed.returncode != 0:
        raise FitykRuntimeError(
            f"cfityk exited {completed.returncode}: {completed.stderr or completed.stdout}"
        )
    return completed
