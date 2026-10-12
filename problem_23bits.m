%% Example problem

%% System
NAg   = 20;   % Number of agents
dim   = 2;    % Dimension of the local systems
Amag  = 0.1;  % Magnitude of the dynamics matrix perturbation
Cmag  = 1;    % Magnitude of the measurement matrix perturbation

% process noise
wmag  = 1;    % Noise magnitude
% measurement noise
vmag  = 1;    % Measurement noise magnitude
% initial uncertainty
x0mag = 1000;

% Observer parameters
% process noise
wobs  = 1;    % Noise magnitude
% measurement noise
vobs  = 1;    % Measurement noise magnitude
% initial uncertainty
x0obs = 1000;

%% Network
Circ_order = 1; % communicate with circ_order neighbours on each side

%% Final iteration
ifinal = 100;

%% Number of samples
samples = 1;

%% Design parameters
beta = 0.5;
alpha = 10^-6;
kfix = ceil((NAg/Circ_order)+dim);

%% Quantization parameters
use_quantization = true;
nbits = 23; % grid exponent p; fixed-length cost is b=24 bits/scalar
gammaOmega = 1.0;
gammaMsg = 1.0;
eps_quant = 1e-12;

%% Numerical safeguard before solving the quantized information equation
% This is distinct from alpha, which regularizes the Gramian recursion.
reg_epsilon = 1e-8;

%% Debug
check_equivalence = 0;
