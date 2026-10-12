# Common-range numerical diagnostics

Run in MATLAB from this repository:

```matlab
run_diagnostics
```

The script reuses the original seeded plant, graph, initialization, and noise realization for quantizer parameters `p=20:26`. It records continuous finite-horizon diagnostics without assigning the labels “bounded”, “borderline”, or “divergent”.

Outputs are written to `diagnostics_results/`:

- `diagnostics_summary.csv`: estimation error, normalized error, late-window growth, conditioning, regularization, residual, solver, non-finite, and transmitted-saturation diagnostics.
- `timeseries_p23.csv`, `timeseries_p24.csv`, and `timeseries_p25.csv`: aligned error and matrix diagnostics.
- `case_pXX.mat`: complete logged output for each tested value.
- `diagnostics_error_conditioning.pdf`: comparison of error, condition number, and minimum eigenvalue.
- `diagnostics_payload.csv`: payload calculation for the original full-matrix-plus-vector message.
- `diagnostics_configuration.csv`: exact seed, distributions, dimensions, graph, and numerical parameters.

## Interpretation limits

These diagnostics characterize the common-range control on one fixed 100-sample realization. Independent realizations, extended horizons and a matched finite-bit reference are documented in `REVIEW_VALIDATION.md`. Neither campaign certifies asymptotic stability or performance on arbitrary graphs and plants.

The quantizer uses step `2X/2^p` and includes both endpoints, giving `2^p+1` levels. The payload file therefore distinguishes the grid parameter `p` from the `p+1` fixed-length bits required to encode every level index. The common-range control's payload excludes the bits needed to make the two global ranges available to receivers. The sender-local validation reports symmetric packet sizes including two binary64 ranges.

The saturation rate counts only scalar values sent over external nonzero graph arcs. Quantization of a node's own contribution is logged as a numerical operation but excluded from communication saturation and payload.

`alpha` regularizes the Gramian recursion. `reg_epsilon` is a separate numerical safeguard applied to the aggregated quantized information matrix before solving. The diagnostics report its activation and magnitude.

## Finite-horizon invertibility check

The diagnostic run also logs the minimum eigenvalue of the corresponding unquantized aggregate information matrix. For each `p`, it evaluates the finite-horizon quantities

```text
lhs = n * max_t XmaxOmega(t) / 2^p
rhs = min_{t,i} lambda_min(Omega_nominal(t,i)) / 2
```

and records whether `lhs <= rhs`. This is an empirical check over the stored 100-sample realization. It is not a certification of the uniform infinite-horizon constants assumed by the theorem.
