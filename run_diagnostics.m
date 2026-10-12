%% Common-range numerical diagnostics
% Reuses the original plant, graph, initialization, and noise realization.
% It does not constitute a Monte Carlo study or an asymptotic stability test.

clc;
close all;
addpath('./util');
set(0,'DefaultFigureVisible','off');

repoDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(repoDir,'diagnostics_results');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

bit_values = 20:26;
diagnostic_bits = [23 24 25];
n_cases = numel(bit_values);

% Build the shared numerical realization once. Resetting to sample_rng_state
% before each case gives every resolution the same initial condition and noise.
rng(7,'twister');
problem;
gen_networks;
setup_script;
sample_rng_state = rng;

mean_error = nan(n_cases,1);
final_error = nan(n_cases,1);
max_error = nan(n_cases,1);
mean_normalized_error = nan(n_cases,1);
final_normalized_error = nan(n_cases,1);
log10_error_slope = nan(n_cases,1);
minimum_eigenvalue_pre_regularization = nan(n_cases,1);
minimum_nominal_eigenvalue = nan(n_cases,1);
maximum_information_matrix_range = nan(n_cases,1);
empirical_invertibility_lhs = nan(n_cases,1);
empirical_invertibility_rhs = nan(n_cases,1);
empirical_invertibility_satisfied = false(n_cases,1);
maximum_condition_number_pre_regularization = nan(n_cases,1);
maximum_condition_number_post_regularization = nan(n_cases,1);
regularization_activations = zeros(n_cases,1);
regularization_rate = zeros(n_cases,1);
maximum_regularization_amount = zeros(n_cases,1);
maximum_relative_residual = nan(n_cases,1);
solver_failures = zeros(n_cases,1);
nonfinite_solutions = zeros(n_cases,1);
transmitted_saturations = zeros(n_cases,1);
transmitted_scalar_count = zeros(n_cases,1);
transmitted_saturation_rate = nan(n_cases,1);
all_error_values_finite = false(n_cases,1);

diagnostic_curves = cell(numel(diagnostic_bits),1);

for idx = 1:n_cases
    current_bits = bit_values(idx);
    fprintf('Running common-range diagnostic case with quantizer parameter p=%i\n',current_bits);

    rng(sample_rng_state);
    use_quantization = true;
    nbits = current_bits;
    simulation_OBS;

    active_idx = 2:ifinal;
    error_curve = nedistg(:);
    state_norm = sqrt(sum(sim_output.xlog.^2,1))';
    normalized_error = error_curve ./ max(1,state_norm);

    mean_error(idx) = mean(error_curve);
    final_error(idx) = error_curve(end);
    max_error(idx) = max(error_curve);
    mean_normalized_error(idx) = mean(normalized_error);
    final_normalized_error(idx) = normalized_error(end);
    all_error_values_finite(idx) = all(isfinite(error_curve));

    window_start = max(2,ifinal-19);
    window_idx = (window_start:ifinal)';
    fit_coefficients = polyfit(window_idx,log10(max(error_curve(window_idx),realmin)),1);
    log10_error_slope(idx) = fit_coefficients(1);

    min_eig_values = sim_output.minEigPreAgent_log(active_idx,:);
    nominal_eig_values = sim_output.minEigNominalAgent_log(active_idx,:);
    cond_pre_values = sim_output.condPreAgent_log(active_idx,:);
    cond_post_values = sim_output.condPostAgent_log(active_idx,:);
    residual_values = sim_output.relativeResidualAgent_log(active_idx,:);
    reg_values = sim_output.regularizationAgent_log(active_idx,:);

    minimum_eigenvalue_pre_regularization(idx) = finite_min(min_eig_values);
    minimum_nominal_eigenvalue(idx) = finite_min(nominal_eig_values);
    maximum_information_matrix_range(idx) = max(sim_output.XmaxOmega_log(active_idx));
    empirical_invertibility_lhs(idx) = nx * maximum_information_matrix_range(idx) / 2^current_bits;
    empirical_invertibility_rhs(idx) = minimum_nominal_eigenvalue(idx) / 2;
    empirical_invertibility_satisfied(idx) = ...
        empirical_invertibility_lhs(idx) <= empirical_invertibility_rhs(idx);
    maximum_condition_number_pre_regularization(idx) = finite_max(cond_pre_values);
    maximum_condition_number_post_regularization(idx) = finite_max(cond_post_values);
    regularization_activations(idx) = nnz(reg_values > 0);
    regularization_rate(idx) = regularization_activations(idx) / numel(reg_values);
    maximum_regularization_amount(idx) = finite_max([0; reg_values(:)]);
    maximum_relative_residual(idx) = finite_max(residual_values);
    solver_failures(idx) = nnz(sim_output.solverFailureAgent_log(active_idx,:));
    nonfinite_solutions(idx) = nnz(sim_output.nonfiniteSolutionAgent_log(active_idx,:));

    transmitted_saturations(idx) = sum( ...
        sim_output.satOmegaTransmitted_log(active_idx) + ...
        sim_output.satMsgTransmitted_log(active_idx));
    transmitted_scalar_count(idx) = sum(sim_output.transmittedScalars_log(active_idx));
    if transmitted_scalar_count(idx) > 0
        transmitted_saturation_rate(idx) = ...
            transmitted_saturations(idx) / transmitted_scalar_count(idx);
    end

    case_file = fullfile(outputDir,sprintf('case_p%02i.mat',current_bits));
    save(case_file,'sim_output','error_curve','normalized_error','current_bits','Pi');

    selected_idx = find(diagnostic_bits == current_bits,1);
    if ~isempty(selected_idx)
        time = (1:ifinal)';
        min_eigenvalue_pre = row_min(min_eig_values,active_idx,ifinal);
        min_eigenvalue_nominal = row_min(nominal_eig_values,active_idx,ifinal);
        max_condition_number_pre = row_max(cond_pre_values,active_idx,ifinal);
        max_condition_number_post = row_max(cond_post_values,active_idx,ifinal);
        regularization_count = zeros(ifinal,1);
        regularization_count(active_idx) = sum(reg_values > 0,2);
        max_relative_residual = row_max(residual_values,active_idx,ifinal);

        time_metrics = table(time,error_curve,normalized_error,min_eigenvalue_nominal, ...
            min_eigenvalue_pre, ...
            max_condition_number_pre,max_condition_number_post, ...
            regularization_count,max_relative_residual);
        writetable(time_metrics,fullfile(outputDir, ...
            sprintf('timeseries_p%02i.csv',current_bits)));
        diagnostic_curves{selected_idx} = time_metrics;
    end
end

summary = table(bit_values(:),mean_error,final_error,max_error, ...
    mean_normalized_error,final_normalized_error,log10_error_slope, ...
    minimum_eigenvalue_pre_regularization, ...
    minimum_nominal_eigenvalue,maximum_information_matrix_range, ...
    empirical_invertibility_lhs,empirical_invertibility_rhs, ...
    empirical_invertibility_satisfied, ...
    maximum_condition_number_pre_regularization, ...
    maximum_condition_number_post_regularization, ...
    regularization_activations,regularization_rate, ...
    maximum_regularization_amount,maximum_relative_residual, ...
    solver_failures,nonfinite_solutions,transmitted_saturations, ...
    transmitted_scalar_count,transmitted_saturation_rate, ...
    all_error_values_finite, ...
    'VariableNames',{'quantizer_parameter_p','mean_error','final_error', ...
    'max_error','mean_normalized_error','final_normalized_error', ...
    'late_log10_error_slope','minimum_eigenvalue_pre_regularization', ...
    'minimum_nominal_eigenvalue','maximum_information_matrix_range', ...
    'empirical_invertibility_lhs','empirical_invertibility_rhs', ...
    'empirical_invertibility_satisfied', ...
    'maximum_condition_number_pre_regularization', ...
    'maximum_condition_number_post_regularization', ...
    'regularization_activations','regularization_rate', ...
    'maximum_regularization_amount','maximum_relative_residual', ...
    'solver_failures','nonfinite_solutions','transmitted_saturations', ...
    'transmitted_scalar_count','transmitted_saturation_rate', ...
    'all_error_values_finite'});
writetable(summary,fullfile(outputDir,'diagnostics_summary.csv'));

% The implemented grid has 2^p+1 levels. A fixed-length representation of
% all level indices therefore needs p+1 bits. Range metadata is excluded
% because the current implementation obtains the ranges globally.
quantizer_parameter_p = bit_values(:);
quantization_levels = 2.^quantizer_parameter_p + 1;
fixed_length_bits_per_scalar = ceil(log2(quantization_levels));
scalars_per_message = repmat(nx^2 + nx,n_cases,1);
payload_bits_per_message_excluding_range = ...
    scalars_per_message .* fixed_length_bits_per_scalar;
external_directed_arcs = nnz(Pi) - nnz(diag(Pi));
network_payload_bits_per_step_excluding_range = ...
    payload_bits_per_message_excluding_range * external_directed_arcs;
range_metadata_included = false(n_cases,1);

payload = table(quantizer_parameter_p,quantization_levels, ...
    fixed_length_bits_per_scalar,scalars_per_message, ...
    payload_bits_per_message_excluding_range, ...
    repmat(external_directed_arcs,n_cases,1), ...
    network_payload_bits_per_step_excluding_range,range_metadata_included, ...
    'VariableNames',{'quantizer_parameter_p','quantization_levels', ...
    'fixed_length_bits_per_scalar','scalars_per_message', ...
    'payload_bits_per_message_excluding_range','external_directed_arcs', ...
    'network_payload_bits_per_step_excluding_range','range_metadata_included'});
writetable(payload,fullfile(outputDir,'diagnostics_payload.csv'));

parameter = {'rng_seed';'rng_generator';'agents';'local_dimension'; ...
    'global_dimension';'time_samples';'graph';'neighbors_each_side'; ...
    'consensus_weights';'A_local';'C_local';'process_noise'; ...
    'measurement_noise';'initial_state';'initial_estimates'; ...
    'alpha';'beta';'horizon_k';'gamma_Omega';'gamma_m'; ...
    'range_floor';'regularization_floor';'resolution_sweep'};
value = {'7';'twister';'20';'2';'40';'100';'circular';'1'; ...
    'Metropolis';'I_2 + 0.1 N(0,1)';'I_2 + 1.0 N(0,1)'; ...
    'N(0,I_40)';'N(0,I)';'1000 N(0,I_40)';'zero'; ...
    '1e-6';'0.5';'22';'1';'1';'1e-12';'1e-8';'p=20:26'};
configuration = table(parameter,value);
writetable(configuration,fullfile(outputDir,'diagnostics_configuration.csv'));

make_diagnostic_figure(diagnostic_curves,diagnostic_bits,outputDir);
save(fullfile(outputDir,'diagnostics_workspace_summary.mat'),'summary','payload', ...
    'configuration','diagnostic_bits','diagnostic_curves','sample_rng_state');

disp(summary);
disp(payload);
fprintf('Common-range diagnostics saved in %s\n',outputDir);

function value = finite_min(values)
    values = values(isfinite(values));
    if isempty(values)
        value = nan;
    else
        value = min(values);
    end
end

function value = finite_max(values)
    values = values(isfinite(values));
    if isempty(values)
        value = nan;
    else
        value = max(values);
    end
end

function output = row_min(values,active_idx,ifinal)
    output = nan(ifinal,1);
    for row = 1:size(values,1)
        output(active_idx(row)) = finite_min(values(row,:));
    end
end

function output = row_max(values,active_idx,ifinal)
    output = nan(ifinal,1);
    for row = 1:size(values,1)
        output(active_idx(row)) = finite_max(values(row,:));
    end
end

function make_diagnostic_figure(curves,bits,outputDir)
    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];
    fig = figure('Visible','off','Color','w','Position',[100 100 900 900]);

    subplot(3,1,1);
    hold on;
    for idx = 1:numel(bits)
        plot_values = curves{idx}.error_curve;
        plot_values(~isfinite(plot_values) | plot_values <= 0) = nan;
        semilogy(curves{idx}.time,plot_values, ...
            'LineWidth',1.4,'Color',colors(idx,:));
    end
    set(gca,'YScale','log');
    grid on;
    ylabel('Average error norm');
    legend(compose('p=%i',bits),'Location','best');

    subplot(3,1,2);
    hold on;
    for idx = 1:numel(bits)
        plot_values = curves{idx}.max_condition_number_post;
        plot_values(~isfinite(plot_values) | plot_values <= 0) = nan;
        semilogy(curves{idx}.time,plot_values, ...
            'LineWidth',1.4,'Color',colors(idx,:));
    end
    set(gca,'YScale','log');
    grid on;
    ylabel('max_i cond_2 after correction');

    subplot(3,1,3);
    hold on;
    for idx = 1:numel(bits)
        plot(curves{idx}.time,curves{idx}.min_eigenvalue_pre, ...
            'LineWidth',1.4,'Color',colors(idx,:));
    end
    yline(0,'k--');
    grid on;
    xlabel('Time step');
    ylabel('min_i lambda_{min}');

    exportgraphics(fig,fullfile(outputDir,'diagnostics_error_conditioning.pdf'), ...
        'ContentType','vector');
    savefig(fig,fullfile(outputDir,'diagnostics_error_conditioning.fig'));
    close(fig);
end
