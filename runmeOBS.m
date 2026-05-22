%% Quantized distributed estimation
clc;close all;
addpath('./util')

%% Problem definition
problem_23bits

%% Generate networks
gen_networks

%% Setup
setup_script

%% Simulation
simulation_OBS

%% Plot
plot_OBS