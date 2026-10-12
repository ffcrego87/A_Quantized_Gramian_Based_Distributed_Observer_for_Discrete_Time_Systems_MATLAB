# Quantized Gramian Observer MATLAB Simulations

The numerical validation uses sender-local ranges, an actual byte codec, independent realizations, longer horizons, checked safeguard comparisons and a matched finite-precision information-filter reference. See `REVIEW_VALIDATION.md` for commands and the supplied canonical data. The original common-range diagnostics provide a separate reproducible control.

This repository contains the MATLAB scripts used to generate the numerical results for the manuscript
`A Quantized Gramian-Based Distributed Observer for Discrete-Time Systems`.

## Main Scripts

- `run_paper_figures.m` regenerates the manuscript figures with white backgrounds and copies the outputs to the LaTeX manuscript folder.
- `run_bit_sweep.m` runs the sweep over the number of communication bits.
- `run_sensitivity_study.m` runs the parameter sensitivity studies for the quantizer and regularization parameters.
- `run_method_comparison.m` runs the centralized, distributed unquantized, Gramian-based unquantized, and quantized observer comparison.
- `runmeOBS.m` is the original entry point for a single observer simulation.
- `run_diagnostics.m` runs the common-range diagnostics for conditioning, regularization, numerical failures, saturation, and payload. See `NUMERICAL_DIAGNOSTICS.md`.

## Problem Files

- `problem_noquant.m` configures the unquantized reference case.
- `problem_25bits.m` is the legacy filename for grid exponent `p=25`, which uses `b=26` fixed-length bits per scalar.
- `problem_23bits.m` is the legacy filename for grid exponent `p=23`, which uses `b=24` fixed-length bits per scalar.
- `problem.m` is the base problem setup.

## Utilities

The `util` folder contains the observer simulation, quantizer, network generation, plotting, and error computation routines.

## Reproducing the Results

From this folder, run one of the following commands in MATLAB.

```matlab
run_paper_figures
run_bit_sweep
run_sensitivity_study
run_method_comparison
```

From a terminal with MATLAB on the path, the same scripts can be run with:

```powershell
matlab -batch "run_paper_figures"
matlab -batch "run_bit_sweep"
matlab -batch "run_sensitivity_study"
matlab -batch "run_method_comparison"
matlab -batch "run_diagnostics"
```

The common-range diagnostics use `rng(7,'twister')` and reset the saved random-number-generator state before every resolution case. This gives all cases the same plant, graph, initialization, and noise sequence. The generated `.csv`, `.mat`, `.dat`, and PDF files are included for reproducibility and direct use in the manuscript.

## Quantizer code and numerical safeguard

`util/uniform_quantizer_sat.m` implements the grid with exponent `p`. Its signed index is in `-2^(p-1):2^(p-1)`. The transmitted integer code is obtained by adding `2^(p-1)`, so the valid codewords are `0:2^p`; a fixed-length container requires `p+1` bits. `util/uniform_quantizer_decode.m` implements the inverse map. The remaining `p+1`-bit patterns are unused.

The Gramian initialization uses `alpha=1e-6`. Separately, the linear solve uses `reg_epsilon=1e-8`: if the minimum eigenvalue of the symmetrized aggregate matrix satisfies `lambda_min <= reg_epsilon`, the code adds

```matlab
reg_amount = (reg_epsilon - lambda_min) + reg_epsilon;
```

before solving. The diagnostic output records every activation and the correction magnitude.
