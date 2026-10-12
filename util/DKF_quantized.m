function [xhat,Omeganextcell,qlog] = DKF_quantized( ...
    xpred,Pi,y,yidxs,Omegacell,Cglob,V,DeltaOmega,XmaxOmega,DeltaMsg,XmaxMsg,reg_eps,pack_symmetric)

if nargin < 12
    reg_eps = 1e-8;
end
if nargin < 13
    pack_symmetric = false;
end
if isscalar(DeltaOmega)
    DeltaOmega = repmat(DeltaOmega,length(Omegacell),1);
    XmaxOmega = repmat(XmaxOmega,length(Omegacell),1);
    DeltaMsg = repmat(DeltaMsg,length(Omegacell),1);
    XmaxMsg = repmat(XmaxMsg,length(Omegacell),1);
end

NAg = length(Omegacell);
nx = length(xpred)/NAg;

xhat = xpred;
Omeganextcell = Omegacell;

qlog = struct;
qlog.sat_count_Omega_total = 0;
qlog.sat_count_msg_total = 0;
qlog.sat_count_Omega_transmitted = 0;
qlog.sat_count_msg_transmitted = 0;
qlog.quantized_Omega_scalars_total = 0;
qlog.quantized_msg_scalars_total = 0;
qlog.transmitted_Omega_scalars_total = 0;
qlog.transmitted_msg_scalars_total = 0;
qlog.regularization_count_total = 0;
qlog.range_metadata_bits_transmitted = 0;
qlog.codec_roundtrip_verified = true;
qlog.min_eig_nominal_log = nan(NAg,1);
qlog.min_eig_pre_log = nan(NAg,1);
qlog.cond_pre_log = nan(NAg,1);
qlog.cond_post_log = nan(NAg,1);
qlog.regularization_log = zeros(NAg,1);
qlog.relative_residual_log = nan(NAg,1);
qlog.solver_failed_log = false(NAg,1);
qlog.nonfinite_solution_log = false(NAg,1);

for i = 1:NAg
    Cloc = Cglob(yidxs{i},:);
    Vloc = V(yidxs{i},yidxs{i});
    Omeganext = Cloc' * Vloc * Cloc;
    Omega_nominal = Omeganext;

    idxsx = (i-1)*nx + (1:nx);
    qloc = Cloc' * Vloc * y(yidxs{i});

    for j = 1:NAg
        weight = Pi(i,j);
        if weight == 0
            continue
        end

        idxsxn = (j-1)*nx + (1:nx);

        Omega_raw = Omegacell{j};
        msg_raw = Omegacell{j} * xpred(idxsxn);
        Omega_nominal = Omega_nominal + weight * Omega_raw;

        if pack_symmetric
            [Omegaq,msgq,packet,codec] = quantized_information_pair( ...
                Omega_raw,msg_raw,DeltaOmega(j),XmaxOmega(j), ...
                DeltaMsg(j),XmaxMsg(j));
            satOmega = codec.full_matrix_saturations;
            satMsg = codec.vector_saturations;
            matrix_symbols = numel(packet.matrix_codes);
            transmitted_matrix_saturations = codec.matrix_code_saturations;
            qlog.codec_roundtrip_verified = qlog.codec_roundtrip_verified && codec.roundtrip_verified;
        else
            [Omegaq,satOmega] = uniform_quantizer_sat(Omega_raw,DeltaOmega(j),XmaxOmega(j));
            [msgq,satMsg] = uniform_quantizer_sat(msg_raw,DeltaMsg(j),XmaxMsg(j));
            matrix_symbols = numel(Omega_raw);
            transmitted_matrix_saturations = satOmega;
        end

        qlog.sat_count_Omega_total = qlog.sat_count_Omega_total + satOmega;
        qlog.sat_count_msg_total = qlog.sat_count_msg_total + satMsg;
        qlog.quantized_Omega_scalars_total = ...
            qlog.quantized_Omega_scalars_total + numel(Omega_raw);
        qlog.quantized_msg_scalars_total = ...
            qlog.quantized_msg_scalars_total + numel(msg_raw);

        if i ~= j
            qlog.sat_count_Omega_transmitted = ...
                qlog.sat_count_Omega_transmitted + transmitted_matrix_saturations;
            qlog.sat_count_msg_transmitted = ...
                qlog.sat_count_msg_transmitted + satMsg;
            qlog.transmitted_Omega_scalars_total = ...
                qlog.transmitted_Omega_scalars_total + matrix_symbols;
            qlog.transmitted_msg_scalars_total = ...
                qlog.transmitted_msg_scalars_total + numel(msg_raw);
            if pack_symmetric
                qlog.range_metadata_bits_transmitted = qlog.range_metadata_bits_transmitted + 128;
            end
        end

        qloc = qloc + weight * msgq;
        Omeganext = Omeganext + weight * Omegaq;
    end

    % Force symmetry after quantization
    Omeganext = 0.5 * (Omeganext + Omeganext');
    Omega_nominal = 0.5 * (Omega_nominal + Omega_nominal');

    % Check conditioning / positive definiteness
    qlog.min_eig_nominal_log(i) = min(real(eig(Omega_nominal)));
    mineig = min(real(eig(Omeganext)));
    qlog.min_eig_pre_log(i) = mineig;
    qlog.cond_pre_log(i) = cond(Omeganext,2);

    if reg_eps > 0 && mineig <= reg_eps
        reg_amount = (reg_eps - mineig) + reg_eps;
        Omeganext = Omeganext + reg_amount * eye(nx);
        qlog.regularization_count_total = qlog.regularization_count_total + 1;
        qlog.regularization_log(i) = reg_amount;
    end

    qlog.cond_post_log(i) = cond(Omeganext,2);

    try
        local_solution = Omeganext \ qloc;
    catch
        local_solution = nan(nx,1);
        qlog.solver_failed_log(i) = true;
    end

    qlog.nonfinite_solution_log(i) = any(~isfinite(local_solution));
    if ~qlog.nonfinite_solution_log(i)
        residual_denominator = norm(Omeganext,2) * norm(local_solution,2) + norm(qloc,2);
        qlog.relative_residual_log(i) = ...
            norm(Omeganext * local_solution - qloc,2) / max(residual_denominator,realmin);
    end

    xhat(idxsx) = local_solution;
    Omeganextcell{i} = Omeganext;
end

end
