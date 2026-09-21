"""Headless PeakComb fit and export."""

from __future__ import annotations

import argparse
from pathlib import Path

from .baseline import subtract_in_window
from .io import load_spectrum
from .models import BackgroundKind, FitMode, PeakSeed, SeedMethod, Session
from .pipeline import official_fit, preview_fit
from .export import export_run
from .seed import seed_peaks


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="PeakComb fixed-FWHM XRD peak comb.")
    parser.add_argument("input", help="Spectrum file (.xy/.csv/.chi/.xye/.txt)")
    parser.add_argument("--out", required=True, help="Run output directory")
    parser.add_argument("--mode", choices=["domain", "distribution"], default="domain")
    parser.add_argument("--axis", choices=["d_a", "q_inv_a", "two_theta_deg"], default="d_a")
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--fwhm", type=float, required=True)
    parser.add_argument("--n", type=int, default=None)
    parser.add_argument("--spacing", type=float, default=0.4)
    parser.add_argument("--shape", type=float, default=0.5)
    parser.add_argument("--background", choices=["none", "constant", "linear", "quadratic"], default="none")
    parser.add_argument(
        "--baseline",
        default="none",
        help="Windowed pybaselines method: none/asls/arpls/airpls/snip/modpoly/rubberband",
    )
    parser.add_argument("--comb-xmin", type=float, default=None)
    parser.add_argument("--comb-xmax", type=float, default=None)
    parser.add_argument("--seed", choices=["equal", "maxima", "grid", "manual"], default="equal")
    parser.add_argument("--centers", default="", help="Comma-separated peak centers for manual Domain seeding")
    parser.add_argument("--wavelength", type=float, default=None)
    parser.add_argument("--energy-kev", type=float, default=83.0)
    parser.add_argument("--d0", type=float, default=None)
    parser.add_argument("--preview-only", action="store_true")
    parser.add_argument("--nnls-only", action="store_true")
    parser.add_argument("--cfityk", default="")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    spectrum = load_spectrum(
        args.input,
        axis=args.axis,
        wavelength_a=args.wavelength,
        energy_kev=None if args.wavelength else args.energy_kev,
    )
    peaks = []
    seed_method = args.seed
    n_peaks = args.n if args.n is not None else 6
    if args.mode == "distribution" and not args.centers.strip():
        seed_method = "equal" if args.n is not None else "grid"
    if args.centers.strip():
        centers = [float(item) for item in args.centers.split(",") if item.strip()]
        peaks = [PeakSeed(center=center, height=1.0, label=f"P{i + 1}") for i, center in enumerate(centers)]
        seed_method = "manual"
    if args.baseline and args.baseline != "none":
        xmin = args.xmin if args.xmin is not None else float(spectrum.x.min())
        xmax = args.xmax if args.xmax is not None else float(spectrum.x.max())
        corrected, _baseline, fitted = subtract_in_window(spectrum.x, spectrum.y, xmin, xmax, args.baseline)
        spectrum.y = corrected
        print(f"baseline={fitted.method} backend={fitted.backend}")
    session = Session(
        input_path=str(Path(args.input)),
        axis=args.axis,
        wavelength_a=spectrum.wavelength_a,
        energy_kev=spectrum.energy_kev,
        xmin=args.xmin,
        xmax=args.xmax,
        comb_xmin=args.comb_xmin,
        comb_xmax=args.comb_xmax,
        baseline_method=args.baseline,
        mode=FitMode(args.mode),
        fwhm=args.fwhm,
        shape=args.shape,
        background=BackgroundKind(args.background),
        n_peaks=len(peaks) if peaks else n_peaks,
        seed_method=SeedMethod(seed_method),
        spacing_factor=args.spacing,
        d0=args.d0,
        cfityk_path=args.cfityk,
        peaks=peaks,
    )
    seed_peaks(session, spectrum)
    out = Path(args.out)
    if args.preview_only:
        result = preview_fit(session, spectrum)
    else:
        result = official_fit(session, spectrum, out, force_nnls=args.nnls_only)
    export_run(result, out, session)
    print(f"fitter={result.fitter} rwp={result.qc.rwp:.4g} n={result.qc.n_peaks} out={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
