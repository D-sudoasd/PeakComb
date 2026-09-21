# PeakComb agent notes

- Independent package. Do not import `fityk_flow` or copy PeakTrace source.
- Shared width is a single Fityk variable `$hwhm`. All PseudoVoigt components must reference it.
- Python NNLS is the live preview; `cfityk` is the official fitter. Distribution combs with more than 30 peaks stay on NNLS.
- `area_frac` is an intensity share, not an automatic volume fraction or stress histogram.
- Keep the GUI plot-first. Do not grow this into PeakTrace-style frame tracking.
