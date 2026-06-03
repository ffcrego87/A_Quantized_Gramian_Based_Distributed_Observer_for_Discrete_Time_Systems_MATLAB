%% Sensitivity study for observer and quantizer parameters
clc;
close all;
addpath('./util');
set(0,'DefaultFigureVisible','off');

repoDir = fileparts(mfilename('fullpath'));
resubmissionDir = fileparts(fileparts(repoDir));
paperFigureDir = fullfile(resubmissionDir, ...
    'latex_artigo', ...
    'A_Quantized_Gramian_Based_Distributed_Observer_for_Discrete_Time_Systems');

nbits_nominal = 25;
alpha_nominal = 1e-6;
beta_nominal = 0.5;
gamma_nominal = 1.0;

% Build the plant, graph, and observer data once. The random generator state
% immediately after setup is reused so each parameter choice sees the same
% initial condition and noise realization.
rng(7,'twister');
problem;
use_quantization = true;
nbits = nbits_nominal;
alpha = alpha_nominal;
beta = beta_nominal;
gammaOmega = gamma_nominal;
gammaMsg = gamma_nominal;
gen_networks;
setup_script;
sample_rng_state = rng;

alpha_values = [1e-8, 1e-6, 1e-4];
beta_values = [0.3, 0.5, 0.7];
gamma_values = [1.0, 2.0, 4.0];

alpha_beta_rows = numel(alpha_values) * numel(beta_values);
alpha_beta = zeros(alpha_beta_rows,1);
beta_column = zeros(alpha_beta_rows,1);
gamma_omega_column = gamma_nominal * ones(alpha_beta_rows,1);
gamma_msg_column = gamma_nominal * ones(alpha_beta_rows,1);
mean_error = zeros(alpha_beta_rows,1);
final_error = zeros(alpha_beta_rows,1);
max_error = zeros(alpha_beta_rows,1);
mean_sat_rate = zeros(alpha_beta_rows,1);
max_delta_omega = zeros(alpha_beta_rows,1);
max_delta_msg = zeros(alpha_beta_rows,1);
outcome = strings(alpha_beta_rows,1);

row = 0;
for alpha_idx = 1:numel(alpha_values)
    for beta_idx = 1:numel(beta_values)
        row = row + 1;
        alpha = alpha_values(alpha_idx);
        beta = beta_values(beta_idx);
        gammaOmega = gamma_nominal;
        gammaMsg = gamma_nominal;
        nbits = nbits_nominal;
        Delta_cell = recompute_delta_cell(NAg, Pi, yidxs, kfix, Cglob, beta);

        fprintf('Sensitivity alpha=%g beta=%g gamma=%g\n', alpha, beta, gammaOmega);
        rng(sample_rng_state);
        simulation_OBS;

        [mean_error(row), final_error(row), max_error(row), ...
            mean_sat_rate(row), max_delta_omega(row), max_delta_msg(row), ...
            outcome(row)] = collect_metrics(nedistg, sim_output, ifinal, NAg, nx);

        alpha_beta(row) = alpha;
        beta_column(row) = beta;
    end
end

alpha_beta_metrics = table(alpha_beta, beta_column, gamma_omega_column, ...
    gamma_msg_column, mean_error, final_error, max_error, mean_sat_rate, ...
    max_delta_omega, max_delta_msg, outcome, ...
    'VariableNames', {'alpha','beta','gammaOmega','gammaMsg', ...
    'mean_gramian_error','final_gramian_error','max_gramian_error', ...
    'mean_saturation_rate','max_delta_omega','max_delta_m','outcome'});

gamma_rows = numel(gamma_values);
gamma_alpha = alpha_nominal * ones(gamma_rows,1);
gamma_beta = beta_nominal * ones(gamma_rows,1);
gamma_omega = zeros(gamma_rows,1);
gamma_msg = zeros(gamma_rows,1);
gamma_mean_error = zeros(gamma_rows,1);
gamma_final_error = zeros(gamma_rows,1);
gamma_max_error = zeros(gamma_rows,1);
gamma_mean_sat_rate = zeros(gamma_rows,1);
gamma_max_delta_omega = zeros(gamma_rows,1);
gamma_max_delta_msg = zeros(gamma_rows,1);
gamma_outcome = strings(gamma_rows,1);

for idx = 1:gamma_rows
    alpha = alpha_nominal;
    beta = beta_nominal;
    gammaOmega = gamma_values(idx);
    gammaMsg = gamma_values(idx);
    nbits = nbits_nominal;
    Delta_cell = recompute_delta_cell(NAg, Pi, yidxs, kfix, Cglob, beta);

    fprintf('Sensitivity gammaOmega=gammaMsg=%g\n', gammaOmega);
    rng(sample_rng_state);
    simulation_OBS;

    [gamma_mean_error(idx), gamma_final_error(idx), gamma_max_error(idx), ...
        gamma_mean_sat_rate(idx), gamma_max_delta_omega(idx), ...
        gamma_max_delta_msg(idx), gamma_outcome(idx)] = collect_metrics( ...
        nedistg, sim_output, ifinal, NAg, nx);

    gamma_omega(idx) = gammaOmega;
    gamma_msg(idx) = gammaMsg;
end

gamma_metrics = table(gamma_alpha, gamma_beta, gamma_omega, gamma_msg, ...
    gamma_mean_error, gamma_final_error, gamma_max_error, ...
    gamma_mean_sat_rate, gamma_max_delta_omega, gamma_max_delta_msg, ...
    gamma_outcome, ...
    'VariableNames', {'alpha','beta','gammaOmega','gammaMsg', ...
    'mean_gramian_error','final_gramian_error','max_gramian_error', ...
    'mean_saturation_rate','max_delta_omega','max_delta_m','outcome'});

alpha_beta_file = fullfile(paperFigureDir,'sensitivity_alpha_beta_metrics.csv');
gamma_file = fullfile(paperFigureDir,'sensitivity_gamma_metrics.csv');
writetable(alpha_beta_metrics,alpha_beta_file);
writetable(gamma_metrics,gamma_file);

disp(alpha_beta_metrics);
disp(gamma_metrics);
fprintf('Saved sensitivity metrics to %s\n', alpha_beta_file);
fprintf('Saved gamma sensitivity metrics to %s\n', gamma_file);

function Delta_cell = recompute_delta_cell(NAg, Pi, yidxs, kfix, Cglob, beta)
    Delta_cell = cell(NAg,1);
    block_size = size(Cglob{1},1);
    for i = 1:NAg
        auxdiag = zeros(block_size*(kfix+1),1);
        for k = 0:kfix
            Pik = Pi^k;
            for j = 1:NAg
                row_idx = block_size*k + yidxs{j};
                auxdiag(row_idx,1) = beta^(k+1) * Pik(i,j) * ones(length(yidxs{j}),1);
            end
        end
        Delta_cell{i} = diag(auxdiag);
    end
end

function [mean_error, final_error, max_error, mean_sat_rate, ...
    max_delta_omega, max_delta_msg, outcome] = collect_metrics( ...
    nedistg, sim_output, ifinal, NAg, nx)

    active_idx = 2:ifinal;
    total_quantized_scalars = NAg * NAg * (nx^2 + nx);
    sat_total = sim_output.satOmega_log(active_idx) + sim_output.satMsg_log(active_idx);

    mean_error = mean(nedistg);
    final_error = nedistg(end);
    max_error = max(nedistg);
    mean_sat_rate = mean(sat_total) / total_quantized_scalars;
    max_delta_omega = max(sim_output.DeltaOmega_log(active_idx));
    max_delta_msg = max(sim_output.DeltaMsg_log(active_idx));

    if isfinite(max_error) && max_error < 1e6
        outcome = "bounded";
    elseif isfinite(final_error) && final_error < 1e6
        outcome = "borderline";
    else
        outcome = "divergent";
    end
end
