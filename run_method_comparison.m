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

bit_values = [20:26, 28, 30, 32];
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

methods = strings(3 + numel(bit_values),1);
reference = strings(3 + numel(bit_values),1);
communication_bits = strings(3 + numel(bit_values),1);
quantized_quantity = strings(3 + numel(bit_values),1);
mean_error = zeros(3 + numel(bit_values),1);
final_error = zeros(3 + numel(bit_values),1);
max_error = zeros(3 + numel(bit_values),1);
mean_saturation_rate = nan(3 + numel(bit_values),1);
max_delta_msg = nan(3 + numel(bit_values),1);
outcome = strings(3 + numel(bit_values),1);

fprintf('Running unquantized comparison references\n');
rng(sample_rng_state);
use_quantization = false;
simulation_OBS;

methods(1) = "Centralized information filter";
reference(1) = "centralized";
communication_bits(1) = "real-valued";
quantized_quantity(1) = "none";
[mean_error(1), final_error(1), max_error(1), outcome(1)] = ...
    collect_error_metrics(necent);

methods(2) = "Distributed information filter";
reference(2) = "distributed unquantized";
communication_bits(2) = "real-valued";
quantized_quantity(2) = "none";
[mean_error(2), final_error(2), max_error(2), outcome(2)] = ...
    collect_error_metrics(nedist);

methods(3) = "Gramian-based distributed observer";
reference(3) = "Gramian unquantized";
communication_bits(3) = "real-valued";
quantized_quantity(3) = "none";
[mean_error(3), final_error(3), max_error(3), outcome(3)] = ...
    collect_error_metrics(nedistg);

for idx = 1:numel(bit_values)
    row = 3 + idx;
    current_bits = bit_values(idx);
    fprintf('Running proposed quantized observer with %i communication bits\n', current_bits);

    rng(sample_rng_state);
    use_quantization = true;
    nbits = current_bits;
    gammaOmega = gamma_nominal;
    gammaMsg = gamma_nominal;
    simulation_OBS;

    methods(row) = "Proposed quantized Gramian observer";
    reference(row) = "proposed";
    communication_bits(row) = sprintf('%i', current_bits);
    quantized_quantity(row) = "Gramian matrix and information vector";
    [mean_error(row), final_error(row), max_error(row), outcome(row)] = ...
        collect_error_metrics(nedistg);

    [mean_saturation_rate(row), max_delta_msg(row)] = collect_quantizer_metrics( ...
        sim_output, ifinal, NAg, nx);
end

metrics = table(methods, reference, communication_bits, quantized_quantity, ...
    mean_error, final_error, max_error, mean_saturation_rate, max_delta_msg, outcome, ...
    'VariableNames', {'method','reference','communication_bits', ...
    'quantized_quantity','mean_error','final_error','max_error', ...
    'mean_saturation_rate','max_delta_m','outcome'});

csv_file = fullfile(paperFigureDir,'method_comparison_metrics.csv');
writetable(metrics,csv_file);

disp(metrics);
fprintf('Saved method comparison metrics to %s\n', csv_file);

function [mean_error, final_error, max_error, outcome] = collect_error_metrics(error_curve)
    mean_error = mean(error_curve);
    final_error = error_curve(end);
    max_error = max(error_curve);

    if isfinite(max_error) && max_error < 1e6
        outcome = "bounded";
    elseif isfinite(final_error) && final_error < 1e6
        outcome = "borderline";
    else
        outcome = "divergent";
    end
end

function [mean_saturation_rate, max_delta_msg] = collect_quantizer_metrics( ...
    sim_output, ifinal, NAg, nx)

    active_idx = 2:ifinal;
    total_quantized_scalars = NAg * NAg * (nx^2 + nx);
    sat_total = sim_output.satOmega_log(active_idx) + sim_output.satMsg_log(active_idx);

    mean_saturation_rate = mean(sat_total) / total_quantized_scalars;
    max_delta_msg = max(sim_output.DeltaMsg_log(active_idx));
end
