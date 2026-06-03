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

if check_equivalence
    Omegabar = sim_input.Omegabar;
end

% initialization
Cglobt = Cglob{1};
xdhat = zeros(NAg*nx,1);
xdghat = zeros(NAg*nx,1);
xhat = zeros(nx,1);
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
for ii = 1:NAg
    OmegacellGramian{ii} = eye(nx);
end

Chi = alpha * eye(nx);

% log variables
xdhatlog = zeros(NAg*nx,ifinal);
xdghatlog = zeros(NAg*nx,ifinal);
xhatlog = zeros(nx,ifinal);
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
        maxOmegaVal = 0;
        maxMsgVal = 0;

        for j = 1:NAg
            idxsxn = (j-1)*nx + (1:nx);
            msg_current = OmegacellGramian{j} * xdghat(idxsxn);

            maxOmegaVal = max(maxOmegaVal, norm(OmegacellGramian{j}(:), inf));
            maxMsgVal = max(maxMsgVal, norm(msg_current(:), inf));
        end

        XmaxOmega = gammaOmega * max(maxOmegaVal, eps_quant);
        XmaxMsg = gammaMsg * max(maxMsgVal, eps_quant);
        DeltaOmega = 2 * XmaxOmega / (2^nbits);
        DeltaMsg = 2 * XmaxMsg / (2^nbits);

        [xdghat,~,qlog] = DKF_quantized( ...
            xdghat, Pi, y, yidxs, OmegacellGramian, Cglobt, eye(size(Cglobt,1)), ...
            DeltaOmega, XmaxOmega, DeltaMsg, XmaxMsg);

        XmaxOmega_log(ii) = XmaxOmega;
        XmaxMsg_log(ii) = XmaxMsg;
        DeltaOmega_log(ii) = DeltaOmega;
        DeltaMsg_log(ii) = DeltaMsg;
        satOmega_log(ii) = qlog.sat_count_Omega_total;
        satMsg_log(ii) = qlog.sat_count_msg_total;
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
    xlog(:,ii) = x;
    wlog(:,ii) = w;
    vlog(:,ii) = v;
    ylog(:,ii) = y;
end

sim_output = struct;
sim_output.xdhatlog = xdhatlog;
sim_output.xdghatlog = xdghatlog;
sim_output.xhatlog = xhatlog;
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

end
