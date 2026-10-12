# I2 HODMD analysis and structural fitting

## Requirements and execution

MATLAB R2025b was tested. `step03_fit_I2` requires Optimization Toolbox (`lsqnonlin`).
Open this folder as MATLAB's Current Folder and run:

```matlab
% step01_generate_noise % Optional: overwrite/regenerate the bundled input
step02_run_HODMD
step03_fit_I2
```

The step prefix keeps filenames valid MATLAB identifiers while indicating order.
`step03_fit_I2` recomputes its own HODMD target; it does not need an export from
`step02_run_HODMD`. Their settings are independent. Scripts clear workspace variables
and close figures. Fit outputs are written under `results/` (ignored by Git).

## Files

- `step02_run_HODMD.m`: crop, q moving average, delay embedding, global-beta HODMD
  reconstruction and plots. Default moving-average window=20, delay=20, rank=9.
- `step03_fit_I2.m`: fit distance amplitude and phase. Period and damping time are
  fixed to the selected HODMD eigenvalue. Rgs=2.7 Angstrom, Req=3.3 Angstrom,
  scale=1 and IRF FWHM=100 fs are fixed. Input and theory both use MA20.
- `hodmd_core.m`: numerical helper functions required by `step03_fit_I2` (not a step).
- `step01_generate_noise.m`: optional regeneration of the included
  90% noise input using RNG seed 3 and one realization. Overwrites the matching
  input MAT file. `seedN1` means one realization, not RNG seed 1.
- `data/mock/I2_qDeltaS_noise90_seedN1.mat`: original analysis input; A is
  qDeltaS (Nq x Nt), q is in inverse Angstrom and t is in seconds.

Only crop and q moving average are used as preprocessing. The fit uses the full
included grid; changing analysis crop does not change fit settings automatically.
Global beta is fitted over all delay vectors; diagonal averaging returns the
reconstruction to data space. Dynamic modes group real eigencomponents or
conjugate pairs and rank their summed modal amplitudes.

The mock generator uses 110 cm^-1, damping=400 fs, distance amplitude=0.1
Angstrom, phase=pi and Gaussian IRF FWHM=100 fs on a 1 fs calculation grid.
The generator reproduces the included A exactly. Cleanup preserved the 90%
input's HODMD and fit results. No experimental data are included.


Jovanovic et al. (2014), Sparsity-promoting dynamic mode decomposition.
https://doi.org/10.1063/1.4863670
Only the all-snapshot amplitude formulation is used, not its sparsity penalty.
