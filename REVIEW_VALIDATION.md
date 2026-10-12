# Reproducing the numerical validation

Run MATLAB from this directory. The published plant, measurement setup and
20-agent ring are retained: global state dimension 40, `alpha=1e-6`,
`beta=0.5`, Gramian horizon 22, unit noise scales and zero initial estimates.

```matlab
run_review_validation('all')
run_review_validation('additional')
run_review_validation('selected_precision')
run_review_validation('ablation_checked')
run_review_validation('diagnostic_curves')
finish_review_validation
finalize_validation_data
make_review_validation_figure
```

The first command runs the original seed-7 controls, the sender-range
resolution sweep, 300-step extensions and five independent realizations at
each of p=23,24,25. Seeds 1:5 independently regenerate the time-varying plant,
measurement matrices, initial state and process/measurement noise. All methods
and resolutions within a realization use the same physical data. The graph
and observer parameters are fixed. The common-range controls in `diagnostics_results` are preserved.

`additional` checks p=28,30 on seed 7 and p=44 with the solve correction
disabled, using both original global ranges and sender ranges. The p=44 cases
evaluate sufficient invertibility over the recorded finite window.
`selected_precision` checks p=30 on 300 samples and the same five independent
realizations. `ablation_checked` stops each p=23,24,25 uncorrected case at its
first nonfinite estimate, before it can enter the encoder.

The 300-step runs retain the original 100 plant samples and append independent
matrices generated with seed 17007. Restoring the original post-setup noise
generator state preserves the stochastic prefix. The six freshly executed
short/long p=23,24,25 pairs have identical first-100 error curves; the audit
is saved as `paired_prefix_checks.csv`. New cases also record SHA-256 hashes
of their first-100 A/C data and physical state/process-noise/measurement-noise
logs, allowing the physical prefix to be checked separately from estimates.

The campaign uses four MATLAB computational threads. Archived p=23 errors
differ by up to 1.40159e-4 relative from the fresh control because small
floating-point differences are amplified near quantization boundaries.
The deviation is explicitly recorded; the archived-control tolerance is
0.1%. Fresh paired prefixes are checked at a relative tolerance of 1e-12.

Sender-range cases encode one symmetric matrix triangle and the information
vector with two binary64 range values. The byte serializer/decoder checks
packet reconstruction and equality to canonical full-matrix quantization.
Startup self-checks cover random symmetric and zero matrices with n=2,40 and
p=1,23,26. The fixed code length is p+1 bits; the final byte may need at most
seven zero padding bits. The finite-bit information-filter comparison uses
the same codec, sender ranges, physical data and solve safeguard, applied to
the existing distributed information-filter branch.

In the common-range diagnostic summary, pre-correction extrema were aggregated over
finite values: the field named maximum_condition_number_pre_regularization
therefore records the largest finite condition number. The original full
logs also contain infinite pre-correction condition numbers (117 entries at
p=23 and 67 at p=25); their true maxima are infinite. The audit is included as
condition_infinite_audit.csv. All corrected post-correction condition numbers
in those original controls are finite. The current checked ablations retain
infinite post-correction values in their reported maximum and record explicit
pre/post infinite counts. Missing counts for older checkpoints are NaN.

Errors are Euclidean global-state errors averaged across the agents. Their
normalized form divides each time sample by max(1,norm(actual_state)). Across
five independent realizations, reported ranges are observed minima/maxima,
with means, standard deviations and medians; they are not confidence
intervals. All large low-resolution outcomes are retained. Finite arithmetic
does not establish bounded error or asymptotic stability.

`validation_cases.csv` in the supplied data package collects the valid
campaign records. Each case has a MAT checkpoint and CSV time series. Older
unchecked `safeguard_off` records, in which nonfinite values could enter the
encoder, are archival diagnostics and are excluded from submission data;
only `safeguard_off_checked` is used for ablation conclusions. Checkpoint names
for the current default command request the checked version, so these older
records are never reused as current ablations.

The canonical table contains 47 distinct controlled scenarios, including the
three failed, explicitly stopped ablations. Three diagnostic re-executions
verify the current byte codec and collect conditioning profiles; they are
identified separately and are not counted as independent stochastic samples.
