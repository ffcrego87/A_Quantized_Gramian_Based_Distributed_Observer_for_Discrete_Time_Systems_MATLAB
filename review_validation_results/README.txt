NUMERICAL VALIDATION: REPRODUCIBLE CHECKS

Run run_review_validation('all') from the MATLAB repository directory.
Completed individual cases are reused as checkpoints. The common-range diagnostics_results
are read as controls and never overwritten. Use 'additional' for p=28,30 and
the p=44 invertibility examples, and 'selected_precision' for p=30 on the
longer horizon and five independent realizations. Use 'ablation_checked' for
the safeguard tests with explicit termination at the first nonfinite estimate.
Run finish_review_validation after the initial campaign for paired-prefix checks.

All cases retain the published 20-agent ring, n=40, alpha=1e-6, beta=0.5,
Gramian horizon k=22, unit noise scales, initial-state scale 1000, and zero
initial estimates. MATLAB uses four computational threads for this campaign.

The primary realization is the published rng(7,'twister') 100-step example.
The independent campaign uses seeds 1:5 to regenerate the random time-varying
plant, measurement matrices, initial state and process/measurement noise.
All compared resolutions and methods within a realization reuse precisely
the same plant, initial state and noise sequence.

The 300-step example preserves the original first 100 plant samples. Further
plant and measurement matrices are appended with rng(17007,'twister'). The
noise generator is reset to the original post-setup state, so the complete
first 100 samples are paired with the short run. The strict paired-prefix
check compares freshly executed short and long cases. Archived p=23 curves
may differ in the last digits because linear algebra arithmetic and the
quantization boundary can amplify rounding differences; these deviations are
reported in reproduction_max_relative_error_difference, with a 0.1% tolerance.

Global-range controls retain the original full-matrix encoding. Local-range
cases calculate two ranges at each sender, encode one matrix triangle and the
information vector, and include two binary64 range values in each packet.
Packed reconstruction is checked against canonical full quantization on
every call. Twelve random/zero matrix tests cover n=2,40 and p=1,23,26.

The finite-bit reference applies the same range, codec and solve safeguard to
the existing distributed information filter branch. It is a matched internal
comparison; it is not claimed to reproduce an external published estimator.

Reported errors are the time average/final value of the average across agents
of the Euclidean global-state error. Normalized error divides this quantity
at each time by max(1,norm(actual_state)). Reported intervals are observed
minimum and maximum across five independent realizations, not confidence
intervals. Finite numerical outputs do not imply bounded or stable errors.
Large low-resolution errors are retained as resolution-boundary results.

cases_all.csv contains the complete original campaign; summary_all.csv contains
its independent-realization aggregates. Additional and selected-precision
campaigns have separate case and summary files. Individual MAT/CSV files
contain metrics and curves. The three original uncorrected cases are retained
for provenance; use cases_ablation_checked.csv for final ablation conclusions.
