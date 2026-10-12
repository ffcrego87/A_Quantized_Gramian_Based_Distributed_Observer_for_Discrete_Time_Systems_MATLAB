function sim_output = sim_sample(sim_input)

ifinal = sim_input.ifinal;
NAg = sim_input.NAg;
nx = sim_input.nx;
ny = sim_input.ny;
A = sim_input.A;
Pi = sim_input.Pi;
Cglob = sim_input.Cglob;
yidxs = sim_input.yidxs;
wmag = sim_input.wmag;
vmag = sim_input.vmag;
x0mag = sim_input.x0mag;
wobs = sim_input.wobs;
vobs = sim_input.vobs;
x0obs = sim_input.x0obs;
E_cell = sim_input.E_cell;
D_cell = sim_input.D_cell;
Delta_cell = sim_input.Delta_cell;
kfix = sim_input.kfix;
beta = sim_input.beta;
alpha = sim_input.alpha;
check_equivalence = sim_input.check_equivalence;

use_quantization = sim_input.use_quantization;
nbits = sim_input.nbits;
gammaOmega = sim_input.gammaOmega;
gammaMsg = sim_input.gammaMsg;
eps_quant = sim_input.eps_quant;
reg_epsilon = sim_input.reg_epsilon;
range_mode = 'global';
pack_symmetric = false;
quantize_reference = false;
if isfield(sim_input,'range_mode'), range_mode = sim_input.range_mode; end
if isfield(sim_input,'pack_symmetric'), pack_symmetric = sim_input.pack_symmetric; end
if isfield(sim_input,'quantize_reference'), quantize_reference = sim_input.quantize_reference; end

if check_equivalence
    Omegabar = sim_input.Omegabar;
end

% initialization
Cglobt = Cglob{1};
xdhat = zeros(NAg*nx,1);
xdghat = zeros(NAg*nx,1);
xhat = zeros(nx,1);
xdqhat = zeros(NAg*nx,1);
Omega = (1/x0obs)^2 * eye(nx);
W = (1/wobs)^2 * eye(nx);
V = (1/vobs)^2 * eye(size(Cglob{1},1));

x = x0mag * randn(nx,1);
w = wmag * randn(nx,1);
v = vmag * randn((2*NAg-1)*ny,1);
y = Cglobt*x + v;

T_cell = cell(NAg,1);
for ii = 1:NAg
    T_cell{ii} = zeros(size(Cglobt,1),nx);
end

O_cell = cell(NAg,1);
for ii = 1:NAg
    O_cell{ii} = zeros(size(Cglob{1},1)*(kfix+1),nx);
end

Omegacell = cell(NAg,1);
for ii = 1:NAg
    Omegacell{ii} = (1/x0obs)^2 * eye(nx);
end

OmegacellGramian = cell(NAg,1);
OmegacellReference = Omegacell;
for ii = 1:NAg
    OmegacellGramian{ii} = eye(nx);
end

Chi = alpha * eye(nx);

% log variables
xdhatlog = zeros(NAg*nx,ifinal);
xdghatlog = zeros(NAg*nx,ifinal);
xhatlog = zeros(nx,ifinal);
xdqhatlog = zeros(NAg*nx,ifinal);
xlog = zeros(nx,ifinal);
wlog = zeros(nx,ifinal);
vlog = zeros((2*NAg-1)*ny,ifinal);
ylog = zeros((2*NAg-1)*ny,ifinal);

% quantization logs
XmaxOmega_log = zeros(ifinal,1);
XmaxMsg_log = zeros(ifinal,1);
DeltaOmega_log = zeros(ifinal,1);
DeltaMsg_log = zeros(ifinal,1);
satOmega_log = zeros(ifinal,1);
satMsg_log = zeros(ifinal,1);
satOmegaTransmitted_log = zeros(ifinal,1);
satMsgTransmitted_log = zeros(ifinal,1);
quantizedScalars_log = zeros(ifinal,1);
transmittedScalars_log = zeros(ifinal,1);
minEigPreAgent_log = nan(ifinal,NAg);
minEigNominalAgent_log = nan(ifinal,NAg);
condPreAgent_log = nan(ifinal,NAg);
condPostAgent_log = nan(ifinal,NAg);
regularizationAgent_log = zeros(ifinal,NAg);
relativeResidualAgent_log = nan(ifinal,NAg);
solverFailureAgent_log = false(ifinal,NAg);
nonfiniteSolutionAgent_log = false(ifinal,NAg);
rangeMetadataBits_log = zeros(ifinal,1);
codecRoundtripVerified_log = true(ifinal,1);
referenceRegularizationAgent_log = zeros(ifinal,NAg);
referenceNonfiniteSolutionAgent_log = false(ifinal,NAg);
termination_step = ifinal;
termination_reason = 'completed';
XmaxOmegaAgent_log = zeros(ifinal,NAg);
XmaxMsgAgent_log = zeros(ifinal,NAg);

xdhatlog(:,1) = xdhat;
xdghatlog(:,1) = xdghat;
xhatlog(:,1) = xhat;
xlog(:,1) = x;
wlog(:,1) = w;
vlog(:,1) = v;
ylog(:,1) = y;

% cycle
for ii = 2:ifinal
    At = A{ii};

    % Distributed observer
    [xdhat,Omegacell] = DKF(xdhat,Pi,y,yidxs,Omegacell,Cglobt,V);
    xdhat = kron(eye(NAg),At) * xdhat;

    for j = 1:NAg
        Omegacell{j} = W - W*At*((Omegacell{j} + At'*W*At)^(-1))*At'*W;
    end

    if quantize_reference
        [refDO,refXO,refDM,refXM] = information_pair_ranges( ...
            OmegacellReference,xdqhat,nbits,gammaOmega,gammaMsg,eps_quant,'local');
        [xdqhat,OmegacellReference,referenceLog] = DKF_quantized( ...
            xdqhat,Pi,y,yidxs,OmegacellReference,Cglobt,V, ...
            refDO,refXO,refDM,refXM,reg_epsilon,true);
        xdqhat = kron(eye(NAg),At) * xdqhat;
        for j = 1:NAg
            information = OmegacellReference{j};
            OmegacellReference{j} = W - W*At*((information + At'*W*At)\(At'*W));
        end
        referenceRegularizationAgent_log(ii,:) = referenceLog.regularization_log(:)';
        referenceNonfiniteSolutionAgent_log(ii,:) = referenceLog.nonfinite_solution_log(:)';
    end

    % Distributed observer gramian
    Chi = beta * (At^(-1))' * Chi * At^(-1);
    T_prev_cell = T_cell;
    O_prev_cell = O_cell;

    for j = 1:NAg
        % T update
        T_cell{j} = 0 * T_cell{j};
        for k = 1:NAg
            if j == k
                T_cell{j} = T_cell{j} + D_cell{j,k} * Cglobt * At^(-1);
            else
                T_cell{j} = T_cell{j} + D_cell{j,k} * T_prev_cell{k} * At^(-1);
            end
        end

        % O update
        O_update = [zeros(size(Cglobt)); O_prev_cell{j}] * At^(-1);
        O_cell{j} = O_update(1:(size(Cglob{1},1)*(kfix+1)),:) + E_cell{j} * T_cell{j};

        OmegacellGramian{j} = O_prev_cell{j}' * Delta_cell{j} * O_prev_cell{j} + Chi;

        if check_equivalence
            if max(max(abs(OmegacellGramian{j} - Omegabar{j,ii}))) > 1e-5
                fprintf('Disagreement at iteration %i, agent %i\n', ii, j);
            end
        end
    end

    % Adaptive quantization parameters from current Gramian/messages
    if use_quantization
        [DeltaOmega,XmaxOmega,DeltaMsg,XmaxMsg] = information_pair_ranges( ...
            OmegacellGramian,xdghat,nbits,gammaOmega,gammaMsg,eps_quant,range_mode);

        [xdghat,~,qlog] = DKF_quantized( ...
            xdghat, Pi, y, yidxs, OmegacellGramian, Cglobt, eye(size(Cglobt,1)), ...
            DeltaOmega, XmaxOmega, DeltaMsg, XmaxMsg, reg_epsilon,pack_symmetric);

        XmaxOmega_log(ii) = max(XmaxOmega);
        XmaxMsg_log(ii) = max(XmaxMsg);
        DeltaOmega_log(ii) = max(DeltaOmega);
        DeltaMsg_log(ii) = max(DeltaMsg);
        XmaxOmegaAgent_log(ii,:) = XmaxOmega(:)';
        XmaxMsgAgent_log(ii,:) = XmaxMsg(:)';
        rangeMetadataBits_log(ii) = qlog.range_metadata_bits_transmitted;
        codecRoundtripVerified_log(ii) = qlog.codec_roundtrip_verified;
        satOmega_log(ii) = qlog.sat_count_Omega_total;
        satMsg_log(ii) = qlog.sat_count_msg_total;
        satOmegaTransmitted_log(ii) = qlog.sat_count_Omega_transmitted;
        satMsgTransmitted_log(ii) = qlog.sat_count_msg_transmitted;
        quantizedScalars_log(ii) = qlog.quantized_Omega_scalars_total + ...
            qlog.quantized_msg_scalars_total;
        transmittedScalars_log(ii) = qlog.transmitted_Omega_scalars_total + ...
            qlog.transmitted_msg_scalars_total;
        minEigNominalAgent_log(ii,:) = qlog.min_eig_nominal_log(:)';
        minEigPreAgent_log(ii,:) = qlog.min_eig_pre_log(:)';
        condPreAgent_log(ii,:) = qlog.cond_pre_log(:)';
        condPostAgent_log(ii,:) = qlog.cond_post_log(:)';
        regularizationAgent_log(ii,:) = qlog.regularization_log(:)';
        relativeResidualAgent_log(ii,:) = qlog.relative_residual_log(:)';
        solverFailureAgent_log(ii,:) = qlog.solver_failed_log(:)';
        nonfiniteSolutionAgent_log(ii,:) = qlog.nonfinite_solution_log(:)';
    else
        [xdghat,~] = DKF(xdghat,Pi,y,yidxs,OmegacellGramian,Cglobt,eye(size(Cglobt,1)));
    end

    xdghat = kron(eye(NAg),At) * xdghat;

    % Centralized observer
    xhat = (Omega + Cglobt' * V * Cglobt)^(-1) * (Omega * xhat + Cglobt' * V * y);
    Omega = Omega + Cglobt' * V * Cglobt;
    xhat = At * xhat;
    Omega = W - W * At * (Omega + At' * W * At)^(-1) * At' * W;

    % Real system
    x = At * x + w;
    w = wmag * randn(nx,1);
    v = vmag * randn((2*NAg-1)*ny,1);
    Cglobt = Cglob{ii};
    y = Cglobt * x + v;

    % log variables
    xdhatlog(:,ii) = xdhat;
    xdghatlog(:,ii) = xdghat;
    xhatlog(:,ii) = xhat;
    xdqhatlog(:,ii) = xdqhat;
    xlog(:,ii) = x;
    wlog(:,ii) = w;
    vlog(:,ii) = v;
    ylog(:,ii) = y;
    if any(~isfinite(xdghat))
        termination_step = ii;
        termination_reason = 'nonfinite Gramian-observer estimate';
        % Stop before a nonfinite estimate can be silently converted to an
        % integer code. Unobserved future samples remain explicitly missing.
        xdghatlog(:,ii+1:end) = nan;
        xdqhatlog(:,ii+1:end) = nan;
        xdhatlog(:,ii+1:end) = nan;
        xhatlog(:,ii+1:end) = nan;
        xlog(:,ii+1:end) = nan;
        break
    end
end

sim_output = struct;
sim_output.xdhatlog = xdhatlog;
sim_output.xdghatlog = xdghatlog;
sim_output.xhatlog = xhatlog;
sim_output.xdqhatlog = xdqhatlog;
sim_output.xlog = xlog;
sim_output.wlog = wlog;
sim_output.vlog = vlog;
sim_output.ylog = ylog;

sim_output.XmaxOmega_log = XmaxOmega_log;
sim_output.XmaxMsg_log = XmaxMsg_log;
sim_output.DeltaOmega_log = DeltaOmega_log;
sim_output.DeltaMsg_log = DeltaMsg_log;
sim_output.satOmega_log = satOmega_log;
sim_output.satMsg_log = satMsg_log;
sim_output.satOmegaTransmitted_log = satOmegaTransmitted_log;
sim_output.satMsgTransmitted_log = satMsgTransmitted_log;
sim_output.quantizedScalars_log = quantizedScalars_log;
sim_output.transmittedScalars_log = transmittedScalars_log;
sim_output.minEigNominalAgent_log = minEigNominalAgent_log;
sim_output.minEigPreAgent_log = minEigPreAgent_log;
sim_output.condPreAgent_log = condPreAgent_log;
sim_output.condPostAgent_log = condPostAgent_log;
sim_output.regularizationAgent_log = regularizationAgent_log;
sim_output.relativeResidualAgent_log = relativeResidualAgent_log;
sim_output.solverFailureAgent_log = solverFailureAgent_log;
sim_output.nonfiniteSolutionAgent_log = nonfiniteSolutionAgent_log;
sim_output.reg_epsilon = reg_epsilon;
sim_output.range_mode = range_mode;
sim_output.pack_symmetric = pack_symmetric;
sim_output.XmaxOmegaAgent_log = XmaxOmegaAgent_log;
sim_output.XmaxMsgAgent_log = XmaxMsgAgent_log;
sim_output.rangeMetadataBits_log = rangeMetadataBits_log;
sim_output.codecRoundtripVerified_log = codecRoundtripVerified_log;
sim_output.referenceRegularizationAgent_log = referenceRegularizationAgent_log;
sim_output.referenceNonfiniteSolutionAgent_log = referenceNonfiniteSolutionAgent_log;
sim_output.termination_step = termination_step;
sim_output.termination_reason = termination_reason;

end
