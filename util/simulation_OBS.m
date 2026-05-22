%% Simulation

sim_input=struct;

sim_input.ifinal = ifinal;
sim_input.NAg = NAg;
sim_input.nx = nx;
sim_input.ny = ny;
sim_input.A = A;
sim_input.Pi = Pi;
sim_input.Cglob = Cglob;
sim_input.yidxs = yidxs;
sim_input.wmag = wmag;
sim_input.vmag = vmag;
sim_input.x0mag = x0mag;
sim_input.vobs = vobs;
sim_input.wobs = wobs;
sim_input.x0obs = x0obs;
sim_input.E_cell = E_cell;
sim_input.D_cell = D_cell;
sim_input.Delta_cell = Delta_cell;
sim_input.kfix = kfix;
sim_input.alpha = alpha;
sim_input.beta = beta;
sim_input.check_equivalence = check_equivalence;

% novos parâmetros da quantização
sim_input.use_quantization = use_quantization;
sim_input.nbits = nbits;
sim_input.gammaOmega = gammaOmega;
sim_input.gammaMsg = gammaMsg;
sim_input.eps_quant = eps_quant;

if check_equivalence
    sim_input.Omegabar = Omegabar;
end 

mnecentlog = zeros(samples,1);
mnedistlog = zeros(samples,1);
mnedistglog = zeros(samples,1);

% opcionais: logs médios da quantização
mean_satOmega_log = zeros(samples,1);
mean_satMsg_log = zeros(samples,1);
mean_DeltaOmega_log = zeros(samples,1);
mean_DeltaMsg_log = zeros(samples,1);
mean_XmaxOmega_log = zeros(samples,1);
mean_XmaxMsg_log = zeros(samples,1);

for i=1:samples
    fprintf('run number: %i\n',i)
    
    sim_output = sim_sample(sim_input);
    
    xdhatlog = sim_output.xdhatlog;
    xdghatlog = sim_output.xdghatlog;
    xhatlog = sim_output.xhatlog;
    xlog = sim_output.xlog;
    
    ecent = xhatlog - xlog;
    edist = xdhatlog - repmat(xlog,NAg,1);
    edistg = xdghatlog - repmat(xlog,NAg,1);
    
    necent = sqrt(sum(ecent.^2,1));
    nedistag = zeros(NAg,ifinal);
    nedistgag = zeros(NAg,ifinal);
    for j=1:NAg
        nedistag(j,:) = sqrt(sum(edist(nx*(j-1)+(1:nx),:).^2,1));
        nedistgag(j,:) = sqrt(sum(edistg(nx*(j-1)+(1:nx),:).^2,1));
    end
    nedist = mean(nedistag,1);
    nedistg = mean(nedistgag,1);
    
    mnecent = mean(necent);
    mnedist = mean(nedist);
    mnedistg = mean(nedistg);
    
    mnecentlog(i) = mnecent;
    mnedistlog(i) = mnedist;
    mnedistglog(i) = mnedistg;

    % logs da quantização
    mean_satOmega_log(i) = mean(sim_output.satOmega_log);
    mean_satMsg_log(i) = mean(sim_output.satMsg_log);
    mean_DeltaOmega_log(i) = mean(sim_output.DeltaOmega_log);
    mean_DeltaMsg_log(i) = mean(sim_output.DeltaMsg_log);
    mean_XmaxOmega_log(i) = mean(sim_output.XmaxOmega_log);
    mean_XmaxMsg_log(i) = mean(sim_output.XmaxMsg_log);
end