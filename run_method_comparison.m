%% Compare centralized, distributed, unquantized Gramian, and quantized observers
clc;
close all;
addpath('./util');
set(0,'DefaultFigureVisible','off');

repoDir = fileparts(mfilename('fullpath'));
resubmissionDir = fileparts(fileparts(repoDir));
paperFigureDir = fullfile(resubmissionDir, ...
    'latex_artigo', ...
    'A_Quantized_Gramian_Based_Distributed_Observer_for_Discrete_Time_Systems');

p_values = [20:26, 28, 30, 32];
gamma_nominal = 1.0;

% Build the plant, graph, and observer data once. The random generator state
% after setup is reused so all methods see the same initial condition and
% noise sequence.
rng(7,'twister');
problem;
use_quantization = false;
gammaOmega = gamma_nominal;
gammaMsg = gamma_nominal;
gen_networks;
setup_script;
sample_rng_state = rng;

n_rows = 3 + numel(p_values);
methods = strings(n_rows,1);
reference = strings(n_rows,1);
quantizer_parameter_p = strings(n_rows,1);
fixed_length_bits_per_scalar = strings(n_rows,1);
quantized_quantity = strings(n_rows,1);
mean_error = zeros(n_rows,1);
final_error = zeros(n_rows,1);
max_error = zeros(n_rows,1);
mean_saturation_rate = nan(n_rows,1);
max_delta_msg = nan(n_rows,1);

fprintf('Running unquantized comparison references\n');
rng(sample_rng_state);
use_quantization = false;
simulation_OBS;

methods(1) = "Centralized information filter";
reference(1) = "centralized";
quantizer_parameter_p(1) = "real-valued";
fixed_length_bits_per_scalar(1) = "real-valued";
quantized_quantity(1) = "none";
[mean_error(1), final_error(1), max_error(1)] = ...
    collect_error_metrics(necent);

methods(2) = "Distributed information filter";
reference(2) = "distributed unquantized";
quantizer_parameter_p(2) = "real-valued";
fixed_length_bits_per_scalar(2) = "real-valued";
quantized_quantity(2) = "none";
[mean_error(2), final_error(2), max_error(2)] = ...
    collect_error_metrics(nedist);

methods(3) = "Gramian-based distributed observer";
reference(3) = "Gramian unquantized";
quantizer_parameter_p(3) = "real-valued";
fixed_length_bits_per_scalar(3) = "real-valued";
quantized_quantity(3) = "none";
[mean_error(3), final_error(3), max_error(3)] = ...
    collect_error_metrics(nedistg);

for idx = 1:numel(p_values)
    row = 3 + idx;
    current_p = p_values(idx);
    fprintf('Running proposed observer with p=%i (b=%i bits/scalar)\n',current_p,current_p+1);

    rng(sample_rng_state);
    use_quantization = true;
    nbits = current_p;
    gammaOmega = gamma_nominal;
    gammaMsg = gamma_nominal;
    simulation_OBS;

    methods(row) = "Proposed quantized Gramian observer";
    reference(row) = "proposed";
    quantizer_parameter_p(row) = sprintf('%i',current_p);
    fixed_length_bits_per_scalar(row) = sprintf('%i',current_p+1);
    quantized_quantity(row) = "Gramian matrix and information vector";
    [mean_error(row), final_error(row), max_error(row)] = ...
        collect_error_metrics(nedistg);

    [mean_saturation_rate(row), max_delta_msg(row)] = collect_quantizer_metrics( ...
        sim_output, ifinal);
end

metrics = table(methods,reference,quantizer_parameter_p,fixed_length_bits_per_scalar,quantized_quantity, ...
    mean_error,final_error,max_error,mean_saturation_rate,max_delta_msg, ...
    'VariableNames', {'method','reference','quantizer_parameter_p','fixed_length_bits_per_scalar', ...
    'quantized_quantity','mean_error','final_error','max_error', ...
    'mean_saturation_rate','max_delta_m'});

csv_file = fullfile(paperFigureDir,'method_comparison_metrics.csv');
writetable(metrics,csv_file);

disp(metrics);
fprintf('Saved method comparison metrics to %s\n', csv_file);

function [mean_error, final_error, max_error] = collect_error_metrics(error_curve)
    mean_error = mean(error_curve);
    final_error = error_curve(end);
    max_error = max(error_curve);
end

function [mean_saturation_rate, max_delta_msg] = collect_quantizer_metrics( ...
    sim_output, ifinal)

    active_idx = 2:ifinal;
    total_transmitted_scalars = sum(sim_output.transmittedScalars_log(active_idx));
    sat_total = sim_output.satOmegaTransmitted_log(active_idx) + ...
        sim_output.satMsgTransmitted_log(active_idx);

    mean_saturation_rate = sum(sat_total) / total_transmitted_scalars;
    max_delta_msg = max(sim_output.DeltaMsg_log(active_idx));
end
