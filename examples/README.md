# Examples

Synthetic spectra used by the tests and the README figures.

| File | What it is | Use with |
|---|---|---|
| `synthetic_domain.xy` | Five equal-FWHM PseudoVoigt subpeaks on a linear background | Domain / 少数分立峰 |
| `synthetic_distribution.xy` | A Gaussian *d* distribution convolved with a fixed-width kernel | Distribution / 密梳拆分 |

Regenerate:

```powershell
py -m peakcomb.synthesize
```

