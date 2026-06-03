%% Regenerate paper figures with a white-background style
clc;
close all;
addpath('./util');
set(0,'DefaultFigureVisible','off');

repoDir = fileparts(mfilename('fullpath'));
resubmissionDir = fileparts(fileparts(repoDir));
paperFigureDir = fullfile(resubmissionDir, ...
    'latex_artigo', ...
    'A_Quantized_Gramian_Based_Distributed_Observer_for_Discrete_Time_Systems');

cases = { ...
    'noquant', 'problem_noquant'; ...
    '25bits',  'problem_25bits'; ...
    '23bits',  'problem_23bits'};

for case_idx = 1:size(cases,1)
    paperFigurePrefix = cases{case_idx,1};
    problem_script = cases{case_idx,2};

    fprintf('Generating paper figures for %s\n', paperFigurePrefix);
    close all;

    rng(7,'twister');
    eval(problem_script);
    gen_networks;
    setup_script;
    simulation_OBS;
    plot_OBS;
end

disp('Paper figures regenerated with white backgrounds.');
