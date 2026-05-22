function [xhat,Omeganextcell,qlog] = DKF_quantized( ...
    xpred,Pi,y,yidxs,Omegacell,Cglob,V,DeltaOmega,XmaxOmega,DeltaMsg,XmaxMsg)

NAg = length(Omegacell);
nx = length(xpred)/NAg;

xhat = xpred;
Omeganextcell = Omegacell;

qlog = struct;
qlog.sat_count_Omega_total = 0;
qlog.sat_count_msg_total = 0;
qlog.regularization_count_total = 0;
qlog.min_eig_log = zeros(NAg,1);
qlog.regularization_log = zeros(NAg,1);

reg_eps = 1e-8;

for i = 1:NAg
    Cloc = Cglob(yidxs{i},:);
    Vloc = V(yidxs{i},yidxs{i});
    Omeganext = Cloc' * Vloc * Cloc;

    idxsx = (i-1)*nx + (1:nx);
    qloc = Cloc' * Vloc * y(yidxs{i});

    for j = 1:NAg
        idxsxn = (j-1)*nx + (1:nx);

        Omega_raw = Omegacell{j};
        msg_raw = Omegacell{j} * xpred(idxsxn);

        [Omegaq,satOmega] = uniform_quantizer_sat(Omega_raw,DeltaOmega,XmaxOmega);
        [msgq,satMsg] = uniform_quantizer_sat(msg_raw,DeltaMsg,XmaxMsg);

        qlog.sat_count_Omega_total = qlog.sat_count_Omega_total + satOmega;
        qlog.sat_count_msg_total = qlog.sat_count_msg_total + satMsg;

        qloc = qloc + Pi(i,j) * msgq;
        Omeganext = Omeganext + Pi(i,j) * Omegaq;
    end

    % Force symmetry after quantization
    Omeganext = 0.5 * (Omeganext + Omeganext');

    % Check conditioning / positive definiteness
    mineig = min(eig(Omeganext));
    qlog.min_eig_log(i) = mineig;

    if mineig <= reg_eps
        reg_amount = (reg_eps - mineig) + reg_eps;
        Omeganext = Omeganext + reg_amount * eye(nx);
        qlog.regularization_count_total = qlog.regularization_count_total + 1;
        qlog.regularization_log(i) = reg_amount;
    end

    xhat(idxsx) = Omeganext \ qloc;
    Omeganextcell{i} = Omeganext;
end

end