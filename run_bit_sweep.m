%% Sweep quantizer grid resolution for the paper numerical table
clc;
close all;
addpath('./util');
set(0,'DefaultFigureVisible','off');

repoDir = fileparts(mfilename('fullpath'));
resubmissionDir = fileparts(fileparts(repoDir));
paperFigureDir = fullfile(resubmissionDir, ...
    'latex_artigo', ...
    'A_Quantized_Gramian_Based_Distributed_Observer_for_Discrete_Time_Systems');

p_values = 20:26;
n_cases = numel(p_values);

mean_error = zeros(n_cases,1);
final_error = zeros(n_cases,1);
max_error = zeros(n_cases,1);
mean_sat_rate = zeros(n_cases,1);
max_delta_omega = zeros(n_cases,1);
max_delta_msg = zeros(n_cases,1);
max_xmax_msg = zeros(n_cases,1);

% Build the plant, graph, and observer data once. The random generator state
% immediately after setup is then reused so each bit value sees the same
% initial condition and noise realization.
rng(7,'twister');
problem;
gen_networks;
setup_script;
sample_rng_state = rng;

for idx = 1:n_cases
    current_p = p_values(idx);
    fprintf('Running grid sweep case with p=%i (b=%i bits/scalar)\n', current_p,current_p+1);

    rng(sample_rng_state);
    use_quantization = true;
    nbits = current_p;
    simulation_OBS;

    active_idx = 2:ifinal;
    total_transmitted_scalars = sum(sim_output.transmittedScalars_log(active_idx));
    sat_total = sim_output.satOmegaTransmitted_log(active_idx) + ...
        sim_output.satMsgTransmitted_log(active_idx);

    mean_error(idx) = mean(nedistg);
    final_error(idx) = nedistg(end);
    max_error(idx) = max(nedistg);
    mean_sat_rate(idx) = sum(sat_total) / total_transmitted_scalars;
    max_delta_omega(idx) = max(sim_output.DeltaOmega_log(active_idx));
    max_delta_msg(idx) = max(sim_output.DeltaMsg_log(active_idx));
    max_xmax_msg(idx) = max(sim_output.XmaxMsg_log(active_idx));

end

fixed_length_bits_per_scalar = p_values(:) + 1;
metrics = table(p_values(:), fixed_length_bits_per_scalar,mean_error, final_error, max_error, ...
    mean_sat_rate, max_delta_omega, max_delta_msg, max_xmax_msg, ...
    'VariableNames', {'quantizer_parameter_p','fixed_length_bits_per_scalar','mean_gramian_error', ...
    'final_gramian_error','max_gramian_error','mean_saturation_rate', ...
    'max_delta_omega','max_delta_m','max_xmax_m'});

csv_file = fullfile(paperFigureDir,'bit_sweep_metrics.csv');
writetable(metrics,csv_file);

disp(metrics);
fprintf('Saved bit-sweep metrics to %s\n', csv_file);
