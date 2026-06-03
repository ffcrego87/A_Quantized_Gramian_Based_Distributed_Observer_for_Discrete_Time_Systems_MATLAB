%% Sweep communication bits for the paper numerical table
clc;
close all;
addpath('./util');
set(0,'DefaultFigureVisible','off');

repoDir = fileparts(mfilename('fullpath'));
resubmissionDir = fileparts(fileparts(repoDir));
paperFigureDir = fullfile(resubmissionDir, ...
    'latex_artigo', ...
    'A_Quantized_Gramian_Based_Distributed_Observer_for_Discrete_Time_Systems');

bit_values = 20:26;
n_cases = numel(bit_values);

mean_error = zeros(n_cases,1);
final_error = zeros(n_cases,1);
max_error = zeros(n_cases,1);
mean_sat_rate = zeros(n_cases,1);
max_delta_omega = zeros(n_cases,1);
max_delta_msg = zeros(n_cases,1);
max_xmax_msg = zeros(n_cases,1);
outcome = strings(n_cases,1);

% Build the plant, graph, and observer data once. The random generator state
% immediately after setup is then reused so each bit value sees the same
% initial condition and noise realization.
rng(7,'twister');
problem;
gen_networks;
setup_script;
sample_rng_state = rng;

for idx = 1:n_cases
    current_bits = bit_values(idx);
    fprintf('Running bit sweep case with %i communication bits\n', current_bits);

    rng(sample_rng_state);
    use_quantization = true;
    nbits = current_bits;
    simulation_OBS;

    active_idx = 2:ifinal;
    total_quantized_scalars = NAg * NAg * (nx^2 + nx);
    sat_total = sim_output.satOmega_log(active_idx) + sim_output.satMsg_log(active_idx);

    mean_error(idx) = mean(nedistg);
    final_error(idx) = nedistg(end);
    max_error(idx) = max(nedistg);
    mean_sat_rate(idx) = mean(sat_total) / total_quantized_scalars;
    max_delta_omega(idx) = max(sim_output.DeltaOmega_log(active_idx));
    max_delta_msg(idx) = max(sim_output.DeltaMsg_log(active_idx));
    max_xmax_msg(idx) = max(sim_output.XmaxMsg_log(active_idx));

    if isfinite(max_error(idx)) && max_error(idx) < 1e6
        outcome(idx) = "bounded";
    elseif isfinite(final_error(idx)) && final_error(idx) < 1e6
        outcome(idx) = "borderline";
    else
        outcome(idx) = "divergent";
    end
end

metrics = table(bit_values(:), mean_error, final_error, max_error, ...
    mean_sat_rate, max_delta_omega, max_delta_msg, max_xmax_msg, outcome, ...
    'VariableNames', {'communication_bits','mean_gramian_error', ...
    'final_gramian_error','max_gramian_error','mean_saturation_rate', ...
    'max_delta_omega','max_delta_m','max_xmax_m','outcome'});

csv_file = fullfile(paperFigureDir,'bit_sweep_metrics.csv');
writetable(metrics,csv_file);

disp(metrics);
fprintf('Saved bit-sweep metrics to %s\n', csv_file);
