"""Lattice strain and hydrostatic internal-stress conversion from comb components."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from .models import Component, Session

D0_MANUAL = "manual"
D0_PEAK_MAX = "peak_max"
D0_AREA_WEIGHTED = "area_weighted"
D0_SOURCES = (D0_MANUAL, D0_PEAK_MAX, D0_AREA_WEIGHTED)

DEFAULT_E_GPA = 80.0
DEFAULT_NU = 0.33


@dataclass
class StrainStressReport:
    d0: float
    d0_source: str
    e_gpa: float
    nu: float
    k_gpa: float
    modulus_gpa: float
    mean_strain: float
    mean_stress_mpa: float
    strain_std: float
    tension_weight: float
    compression_weight: float
    n_live: int
    tautology: bool
    balance_ok: bool
    note: str
    compute_stress: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def hydrostatic_modulus_gpa(e_gpa: float, nu: float, k_gpa: float | None = None) -> float:
    """Return E/(1-2ν), equal to 3K."""
    if k_gpa not in (None, 0.0) and math.isfinite(float(k_gpa)):
        return 3.0 * float(k_gpa)
    denom = 1.0 - 2.0 * float(nu)
    if abs(denom) < 1e-12:
        raise ValueError("ν must not be 0.5 for hydrostatic conversion.")
    return float(e_gpa) / denom


def resolve_d0(components: list[Component], session: Session) -> tuple[float, str]:
    live = [item for item in components if item.area_frac > 0.0 and item.d_a not in (None, 0.0)]
    if not live:
        raise ValueError("No components with d-spacing to define d0.")
    source = (session.d0_source or D0_MANUAL).strip() or D0_MANUAL
    if source == D0_AREA_WEIGHTED:
        d0 = sum(float(item.area_frac) * float(item.d_a) for item in live)
        return d0, D0_AREA_WEIGHTED
    if source == D0_PEAK_MAX:
        strongest = max(live, key=lambda item: item.height)
        return float(strongest.d_a), D0_PEAK_MAX
    if session.d0 in (None, 0.0) or not math.isfinite(float(session.d0)):
        strongest = max(live, key=lambda item: item.height)
        return float(strongest.d_a), D0_PEAK_MAX
    return float(session.d0), D0_MANUAL


def analyze_components(components: list[Component], session: Session) -> StrainStressReport:
    """Fill strain / hydrostatic stress on components and return the summary."""
    d0, source = resolve_d0(components, session)
    e_gpa = float(session.e_gpa)
    nu = float(session.nu)
    k_gpa = float(session.k_gpa) if session.k_gpa not in (None, 0.0) else e_gpa / (3.0 * (1.0 - 2.0 * nu))
    modulus = hydrostatic_modulus_gpa(e_gpa, nu, session.k_gpa)
    live = []
    for item in components:
        if item.d_a in (None, 0.0):
            item.strain = None
            item.sigma_h_mpa = None
            continue
        eps = (float(item.d_a) - d0) / d0
        item.strain = eps
        item.sigma_h_mpa = 1000.0 * modulus * eps if session.compute_stress else None
        live.append(item)
    weights = [max(float(item.area_frac), 0.0) for item in live]
    wsum = sum(weights)
    if wsum <= 0.0:
        raise ValueError("Component area weights are all zero.")
    weights = [w / wsum for w in weights]
    strains = [float(item.strain) for item in live]
    mean_eps = sum(w * eps for w, eps in zip(weights, strains, strict=True))
    var = sum(w * (eps - mean_eps) ** 2 for w, eps in zip(weights, strains, strict=True))
    tension = sum(w for w, eps in zip(weights, strains, strict=True) if eps > 0.0)
    compression = sum(w for w, eps in zip(weights, strains, strict=True) if eps < 0.0)
    tautology = source == D0_AREA_WEIGHTED
    tol = float(session.mean_strain_tol)
    balance_ok = tautology or abs(mean_eps) <= tol
    if tautology:
        note = "d0 取面积加权平均，ε̄=0 是定义，不是卸载平衡的独立检验。"
    elif balance_ok:
        note = f"|ε̄|={abs(mean_eps):.3e} ≤ {tol:.1e}，与卸载、无宏观残余应力的图像相容。"
    else:
        note = (
            f"|ε̄|={abs(mean_eps):.3e} 大于阈值 {tol:.1e}。"
            "可能是宏观残余应力、d0 选错，或多峰模型不够物理。"
        )
    mean_stress = 1000.0 * modulus * mean_eps if session.compute_stress else math.nan
    return StrainStressReport(
        d0=d0,
        d0_source=source,
        e_gpa=e_gpa,
        nu=nu,
        k_gpa=k_gpa,
        modulus_gpa=modulus,
        mean_strain=mean_eps,
        mean_stress_mpa=mean_stress,
        strain_std=math.sqrt(max(var, 0.0)),
        tension_weight=tension,
        compression_weight=compression,
        n_live=len(live),
        tautology=tautology,
        balance_ok=balance_ok,
        note=note,
        compute_stress=bool(session.compute_stress),
    )
