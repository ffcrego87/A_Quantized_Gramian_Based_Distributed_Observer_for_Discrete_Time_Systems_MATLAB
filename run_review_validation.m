function run_review_validation(mode)
%RUN_REVIEW_VALIDATION Paired stochastic, longer-horizon and safeguard checks.
% Common-range diagnostics_results are read as reproduction controls.
% Independent seeds change the random time-varying plant, initial state and
% noise while retaining the graph/design. Longer runs append data and preserve the
% original first 100 samples and original stochastic realization exactly.

if nargin == 0
    mode = 'all';
end
repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir,'util'));
outputDir = fullfile(repoDir,'review_validation_results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
set(0,'DefaultFigureVisible','off');
maxNumCompThreads(4);

% Generate the published seed-7 realization before extending its plant data.
rng(7,'twister');
problem;
setup_script;
originalNoiseState = rng;
baseInput = struct;
inputNames = {'ifinal','NAg','nx','ny','A','Pi','Cglob','yidxs','wmag', ...
    'vmag','x0mag','vobs','wobs','x0obs','E_cell','D_cell','Delta_cell', ...
    'kfix','alpha','beta','check_equivalence','use_quantization','nbits', ...
    'gammaOmega','gammaMsg','eps_quant','reg_epsilon'};
for i = 1:numel(inputNames)
    baseInput.(inputNames{i}) = eval(inputNames{i});
end
baseInput.quantization_range_mode = 'global';
baseInput.range_mode = 'global';
baseInput.safeguard_enabled = true;
baseInput.pack_symmetric = false;
baseInput.quantize_reference = false;

longInput = baseInput;
longInput.ifinal = 300;
% Use a separate generator seed solely for the appended time-varying matrices.
% The prefix A{1:100}, Cglob{1:100} is retained bit for bit.
rng(17007,'twister');
for t = 101:longInput.ifinal
    Aloc = eye(dim) + Amag * randn(dim,dim);
    longInput.A{t} = kron(eye(NAg),Aloc);
    Cloc = eye(dim) + Cmag * randn(dim,dim);
    blocks = cell(NAg,1);
    G = zeros(2,NAg); G(1,1)=1; G(2,1)=1; G(1,2)=-1; G(2,NAg)=-1;
    blocks{1} = kron(G,Cloc);
    for j = 2:NAg-1
        G=zeros(2,NAg); G(1,j)=1; G(2,j)=1; G(1,j+1)=-1; G(2,j-1)=-1;
        blocks{j}=kron(G,Cloc);
    end
    G=zeros(1,NAg); G(NAg)=1;
    blocks{NAg}=kron(G,Cloc);
    longInput.Cglob{t}=cell2mat(blocks);
end

codecTestsPassed=codec_selfcheck();
switch mode
    case 'benchmark'
        cases = table("reproduction",25,7,100,"global",true, ...
            'VariableNames',{'experiment','p','noise_seed','samples','range_mode','safeguard'});
    case 'all'
        cases = empty_case_table();
        cases(end+1,:) = {"local_range",25,7,100,"local",true};
        for p = [23 24 25]
            cases(end+1,:) = {"reproduction",p,7,100,"global",true};
            cases(end+1,:) = {"long_horizon",p,7,300,"global",true};
            if p~=25,cases(end+1,:) = {"local_range",p,7,100,"local",true};end
            cases(end+1,:) = {"long_horizon",p,7,300,"local",true};
            cases(end+1,:) = {"safeguard_off_checked",p,7,100,"local",false};
        end
        for p=[20 21 22 26]
            cases(end+1,:) = {"local_range",p,7,100,"local",true};
        end
        for seed = 1:5
            for p = [23 24 25]
                cases(end+1,:) = {"independent_realization",p,seed,100,"local",true};
            end
        end
    case 'local'
        cases = empty_case_table();
        for p = [23 24 25]
            cases(end+1,:) = {"local_range",p,7,100,"local",true};
        end
    case 'additional'
        cases = empty_case_table();
        for p=[28 30]
            cases(end+1,:)={"higher_precision",p,7,100,"local",true};
        end
        cases(end+1,:)={"invertibility_margin",44,7,100,"global",false};
        cases(end+1,:)={"invertibility_margin",44,7,100,"local",false};
    case 'ablation_checked'
        cases = empty_case_table();
        for p=[23 24 25]
            cases(end+1,:)={"safeguard_off_checked",p,7,100,"local",false};
        end
    case 'selected_precision'
        cases = empty_case_table();
        cases(end+1,:)={"long_horizon",30,7,300,"local",true};
        for seed=1:5
            cases(end+1,:)={"independent_realization",30,seed,100,"local",true};
        end
    case 'diagnostic_curves'
        cases = empty_case_table();
        for p=[23 24 25]
            cases(end+1,:)={"local_diagnostics",p,7,100,"local",true};
        end
    otherwise
        error('Unknown validation mode: %s',mode);
end

records = struct([]);
cachedPlantSeed = nan;
for row = 1:height(cases)
    c = cases(row,:);
    tag=sprintf('%s_%s_p%i_seed%i_T%i',c.experiment,c.range_mode,c.p,c.noise_seed,c.samples);
    caseFile=fullfile(outputDir,[tag '.mat']);
    fprintf('REVIEW CASE %i/%i: %s\n',row,height(cases),tag);
    reuseCheckpoint=exist(caseFile,'file')~=0;
    if reuseCheckpoint && c.experiment=="safeguard_off_checked"
        old=load(caseFile,'metrics');
        reuseCheckpoint=isfield(old.metrics,'condition_pre_infinite_count');
    end
    if reuseCheckpoint
        stored=load(caseFile,'metrics');
        metrics=stored.metrics;
        fprintf('Using completed checkpoint (%g seconds).\n',metrics.runtime_seconds);
    else
        if c.noise_seed ~= 7
            if cachedPlantSeed ~= c.noise_seed
                [independentInput,independentNoiseState]=build_realization(c.noise_seed);
                cachedPlantSeed=c.noise_seed;
            end
            input=independentInput;
            noiseState=independentNoiseState;
        else
            if c.samples == 300, input=longInput; else, input=baseInput; end
            noiseState=originalNoiseState;
        end
        input.nbits=c.p;
        input.quantization_range_mode=char(c.range_mode);
        input.range_mode=char(c.range_mode);
        input.safeguard_enabled=c.safeguard;
        input.pack_symmetric=c.range_mode == "local";
        input.quantize_reference=c.range_mode == "local" && c.p<44;
        if ~c.safeguard
            input.reg_epsilon=0;
        end
        rng(noiseState);
        started=tic;
        sim_output=sim_sample(input);
        [metrics,curves]=collect_metrics(sim_output,input,c);
        metrics.runtime_seconds=toc(started);
        if c.experiment == "reproduction"
            original=load(fullfile(repoDir,'diagnostics_results',sprintf('case_p%02i.mat',c.p)),'error_curve');
            difference=max(abs(curves.error-original.error_curve(:)) ./ max(1,abs(original.error_curve(:))));
            metrics.reproduction_max_relative_error_difference=difference;
            assert(difference < 1e-3,'Original case reproduction failed (relative tolerance 0.1%%): %g',difference);
        elseif c.experiment == "long_horizon" && c.range_mode == "global"
            original=load(fullfile(repoDir,'diagnostics_results',sprintf('case_p%02i.mat',c.p)),'error_curve');
            difference=max(abs(curves.error(1:100)-original.error_curve(:)) ./ max(1,abs(original.error_curve(:))));
            metrics.reproduction_max_relative_error_difference=difference;
            assert(difference < 1e-3,'Long-horizon original prefix failed (relative tolerance 0.1%%): %g',difference);
        end
        save(caseFile,'metrics','curves','c');
        writetable(struct2table(curves),fullfile(outputDir,[tag '.csv']));
        fprintf('Completed %s in %.1f seconds; final error %.5g, normalized %.5g.\n', ...
            tag,metrics.runtime_seconds,metrics.final_error,metrics.final_normalized_error);
    end
    metrics=normalize_metrics(metrics);
    if isempty(records),records=metrics;else,records(end+1)=metrics;end
    writetable(struct2table(records),fullfile(outputDir,['cases_' mode '.csv']));
end
caseMetrics=struct2table(records);
summaries=aggregate_independent(caseMetrics);
if ~isempty(summaries)
    writetable(summaries,fullfile(outputDir,['summary_' mode '.csv']));
end
configuration=struct('plant_seed',7,'append_plant_seed',17007, ...
    'independent_plant_noise_seeds',1:5,'original_noise_rng_state',originalNoiseState, ...
    'agents',baseInput.NAg,'global_dimension',baseInput.nx,'initial_state_scale',baseInput.x0mag, ...
    'alpha',baseInput.alpha,'beta',baseInput.beta,'horizon_k',baseInput.kfix,'gammaOmega',baseInput.gammaOmega, ...
    'gammaMsg',baseInput.gammaMsg,'regularization_floor',baseInput.reg_epsilon, ...
    'fixed_plant_across_noise_seeds',false,'long_plant_prefix_preserved',true, ...
    'codec_selfcheck_cases_passed',codecTestsPassed,'maximum_computation_threads',4);
save(fullfile(outputDir,['validation_' mode '.mat']),'caseMetrics','summaries','configuration');
disp(caseMetrics);
disp(summaries);
end

function metrics=normalize_metrics(metrics)
% Keep pre-extension corrected checkpoints compatible with the current schema.
% Unchecked safeguard-off checkpoint names are never requested by 'all'.
defaults=struct('minimum_nominal_eigenvalue',nan,'maximum_information_matrix_range',nan, ...
    'empirical_invertibility_lhs',nan,'empirical_invertibility_rhs',nan, ...
    'empirical_invertibility_satisfied',false,'pointwise_invertibility_condition_satisfied',false, ...
    'termination_step',nan,'termination_reason',"", ...
    'physical_prefix_sha256',"",'plant_prefix_sha256',"", ...
    'condition_pre_infinite_count',nan,'condition_post_infinite_count',nan);
names=fieldnames(defaults);
for j=1:numel(names)
    if ~isfield(metrics,names{j}),metrics.(names{j})=defaults.(names{j});end
end
metrics=orderfields(metrics);
end

function [input,noiseState]=build_realization(plantSeed)
rng(plantSeed,'twister');
problem;
setup_script;
noiseState=rng;
input=struct;
names={'ifinal','NAg','nx','ny','A','Pi','Cglob','yidxs','wmag','vmag', ...
    'x0mag','vobs','wobs','x0obs','E_cell','D_cell','Delta_cell','kfix', ...
    'alpha','beta','check_equivalence','use_quantization','nbits', ...
    'gammaOmega','gammaMsg','eps_quant','reg_epsilon'};
for j=1:numel(names),input.(names{j})=eval(names{j});end
end

function cases=empty_case_table()
cases=table(strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1),false(0,1), ...
    'VariableNames',{'experiment','p','noise_seed','samples','range_mode','safeguard'});
end

function [m,curves]=collect_metrics(s,input,c)
error=s.xdghatlog-repmat(s.xlog,input.NAg,1);
agentError=zeros(input.NAg,input.ifinal);
for j=1:input.NAg
    agentError(j,:)=sqrt(sum(error((j-1)*input.nx+(1:input.nx),:).^2,1));
end
curve=mean(agentError,1)';
normalized=curve ./ max(1,sqrt(sum(s.xlog.^2,1))');
active=2:input.ifinal;
reg=s.regularizationAgent_log(active,:);
sat=s.satOmegaTransmitted_log(active)+s.satMsgTransmitted_log(active);
res=s.relativeResidualAgent_log(active,:);
minEig=s.minEigPreAgent_log(active,:);
postCond=s.condPostAgent_log(active,:);
m=struct('experiment',c.experiment,'p',c.p,'noise_seed',c.noise_seed, ...
    'samples',c.samples,'range_mode',c.range_mode,'safeguard',c.safeguard, ...
    'mean_error',mean(curve),'final_error',curve(end),'max_error',max(curve), ...
    'mean_normalized_error',mean(normalized),'final_normalized_error',normalized(end), ...
    'last20_mean_normalized_error',mean(normalized(max(2,end-19):end)), ...
    'all_error_values_finite',all(isfinite(curve)), ...
    'all_states_and_estimates_finite',all(isfinite(s.xlog(:))) && all(isfinite(s.xdghatlog(:))), ...
    'solver_failures',nnz(s.solverFailureAgent_log(active,:)), ...
    'nonfinite_solutions',nnz(s.nonfiniteSolutionAgent_log(active,:)), ...
    'regularization_rate',nnz(reg>0)/numel(reg),'maximum_regularization_amount',max(reg(:)), ...
    'maximum_relative_residual',finite_max(res),'minimum_eigenvalue_pre',finite_min(minEig), ...
    'maximum_condition_number_post',maximum_with_infinite(postCond), ...
    'transmitted_saturations',sum(sat),'transmitted_scalars',sum(s.transmittedScalars_log(active)), ...
    'runtime_seconds',nan,'reproduction_max_relative_error_difference',nan);
m.transmitted_saturation_rate=m.transmitted_saturations/m.transmitted_scalars;
m.codec_roundtrip_verified=all(s.codecRoundtripVerified_log(active));
m.range_metadata_bits_per_network_step=max(s.rangeMetadataBits_log(active));
m.reference_mean_error=nan;
m.reference_final_error=nan;
m.reference_mean_normalized_error=nan;
m.reference_final_normalized_error=nan;
m.reference_all_error_values_finite=false;
m.reference_regularization_rate=nan;
m.reference_nonfinite_solutions=nan;
m.condition_pre_infinite_count=nnz(isinf(s.condPreAgent_log(active,:)));
m.condition_post_infinite_count=nnz(isinf(postCond));
m.minimum_nominal_eigenvalue=finite_min(s.minEigNominalAgent_log(active,:));
m.maximum_information_matrix_range=max(s.XmaxOmega_log(active));
m.empirical_invertibility_lhs=input.nx*m.maximum_information_matrix_range/2^input.nbits;
m.empirical_invertibility_rhs=m.minimum_nominal_eigenvalue/2;
m.empirical_invertibility_satisfied=m.empirical_invertibility_lhs<=m.empirical_invertibility_rhs;
m.pointwise_invertibility_condition_satisfied=all( ...
    input.nx*s.XmaxOmega_log(active)/2^input.nbits <= ...
    min(s.minEigNominalAgent_log(active,:),[],2)/2);
m.termination_step=nan;
m.termination_reason="";
if isfield(s,'termination_step'),m.termination_step=s.termination_step;end
if isfield(s,'termination_reason'),m.termination_reason=string(s.termination_reason);end
prefix=1:min(100,input.ifinal);
m.physical_prefix_sha256=sha256_values([s.xlog(:,prefix);s.wlog(:,prefix);s.vlog(:,prefix)]);
plantPrefix=[];
for t=prefix
    plantPrefix=[plantPrefix;input.A{t}(:);input.Cglob{t}(:)]; %#ok<AGROW>
end
m.plant_prefix_sha256=sha256_values(plantPrefix);
referenceCurve=nan(input.ifinal,1);
referenceNormalized=nan(input.ifinal,1);
if input.quantize_reference
    referenceError=s.xdqhatlog-repmat(s.xlog,input.NAg,1);
    referenceAgentError=zeros(input.NAg,input.ifinal);
    for j=1:input.NAg
        referenceAgentError(j,:)=sqrt(sum(referenceError((j-1)*input.nx+(1:input.nx),:).^2,1));
    end
    referenceCurve=mean(referenceAgentError,1)';
    referenceNormalized=referenceCurve ./ max(1,sqrt(sum(s.xlog.^2,1))');
    m.reference_mean_error=mean(referenceCurve);
    m.reference_final_error=referenceCurve(end);
    m.reference_mean_normalized_error=mean(referenceNormalized);
    m.reference_final_normalized_error=referenceNormalized(end);
    m.reference_all_error_values_finite=all(isfinite(referenceCurve));
    m.reference_regularization_rate=nnz(s.referenceRegularizationAgent_log(active,:)>0)/numel(reg);
    m.reference_nonfinite_solutions=nnz(s.referenceNonfiniteSolutionAgent_log(active,:));
end
curves=struct('time',(1:input.ifinal)','error',curve,'normalized_error',normalized, ...
    'regularization_count',sum(s.regularizationAgent_log>0,2), ...
    'reference_error',referenceCurve,'reference_normalized_error',referenceNormalized, ...
    'minimum_eigenvalue_pre',min(s.minEigPreAgent_log,[],2), ...
    'maximum_condition_number_post',max(s.condPostAgent_log,[],2));
end

function out=aggregate_independent(cases)
rows=struct([]);
for p=unique(cases.p(cases.experiment=="independent_realization"))'
    c=cases(cases.experiment=="independent_realization" & cases.p==p,:);
    if isempty(c),continue;end
    summary=struct('p',p,'independent_runs',height(c), ...
        'finite_runs',nnz(c.all_states_and_estimates_finite), ...
        'mean_error_median',median(c.mean_error),'mean_error_min',min(c.mean_error), ...
        'mean_error_max',max(c.mean_error),'mean_error_mean',mean(c.mean_error), ...
        'mean_error_std',std(c.mean_error), ...
        'final_error_median',median(c.final_error),'final_error_min',min(c.final_error), ...
        'final_error_max',max(c.final_error), ...
        'mean_normalized_error_median',median(c.mean_normalized_error), ...
        'mean_normalized_error_mean',mean(c.mean_normalized_error), ...
        'mean_normalized_error_std',std(c.mean_normalized_error), ...
        'mean_normalized_error_min',min(c.mean_normalized_error), ...
        'mean_normalized_error_max',max(c.mean_normalized_error), ...
        'final_normalized_error_median',median(c.final_normalized_error), ...
        'final_normalized_error_min',min(c.final_normalized_error), ...
        'final_normalized_error_max',max(c.final_normalized_error), ...
        'regularization_rate_mean',mean(c.regularization_rate), ...
        'solver_failures',sum(c.solver_failures),'nonfinite_solutions',sum(c.nonfinite_solutions), ...
        'transmitted_saturations',sum(c.transmitted_saturations), ...
        'reference_mean_normalized_error_median',median(c.reference_mean_normalized_error), ...
        'reference_mean_normalized_error_min',min(c.reference_mean_normalized_error), ...
        'reference_mean_normalized_error_max',max(c.reference_mean_normalized_error));
    if isempty(rows),rows=summary;else,rows(end+1)=summary;end
end
if isempty(rows),out=table;else,out=struct2table(rows);end
end

function count=codec_selfcheck()
count=0;
rng(91001,'twister');
for p=[1 23 26]
    for n=[2 40]
        for zeroCase=[false true]
            if zeroCase
                matrix=zeros(n);vector=zeros(n,1);
            else
                raw=randn(n);matrix=0.5*(raw+raw');vector=randn(n,1);
            end
            matrixRange=max(max(abs(matrix(:))),1e-12);
            vectorRange=max(max(abs(vector(:))),1e-12);
            [~,~,packet,d]=quantized_information_pair(matrix,vector, ...
                2*matrixRange/2^p,matrixRange,2*vectorRange/2^p,vectorRange);
            assert(d.roundtrip_verified);
            assert(numel(packet.matrix_codes)==n*(n+1)/2);
            assert(numel(packet.vector_codes)==n);
            count=count+1;
        end
    end
end
fprintf('Codec selfcheck passed %i random/zero size/grid cases.\n',count);
end

function value=finite_max(x)
x=x(isfinite(x));if isempty(x),value=nan;else,value=max(x);end
end
function value=finite_min(x)
x=x(isfinite(x));if isempty(x),value=nan;else,value=min(x);end
end

function value=maximum_with_infinite(x)
x=x(~isnan(x));if isempty(x),value=nan;else,value=max(x);end
end

function hash=sha256_values(values)
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(double(values(:)),'uint8'));
bytes=typecast(digest.digest(),'uint8');
hash=string(lower(reshape(dec2hex(bytes,2)',1,[])));
end
