# Quantized Gramian Observer MATLAB Simulations

This repository contains the MATLAB scripts used to generate the numerical results for the manuscript
`A Quantized Gramian-Based Distributed Observer for Discrete-Time Systems`.

## Main Scripts

- `run_paper_figures.m` regenerates the manuscript figures with white backgrounds and copies the outputs to the LaTeX manuscript folder.
- `run_bit_sweep.m` runs the sweep over the number of communication bits.
- `run_sensitivity_study.m` runs the parameter sensitivity studies for the quantizer and regularization parameters.
- `run_method_comparison.m` runs the centralized, distributed unquantized, Gramian-based unquantized, and quantized observer comparison.
- `runmeOBS.m` is the original entry point for a single observer simulation.

## Problem Files

- `problem_noquant.m` configures the unquantized reference case.
- `problem_25bits.m` configures the stable quantized case with 25 communication bits.
- `problem_23bits.m` configures the unstable quantized case with 23 communication bits.
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
```

The scripts use fixed random seeds where comparisons require the same plant, graph, initialization, and noise sequence. The generated `.dat` files and PDF figures are included for reproducibility and for direct use in the manuscript.
