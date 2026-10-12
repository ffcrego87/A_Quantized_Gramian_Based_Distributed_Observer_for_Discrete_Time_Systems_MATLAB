function finish_review_validation()
%FINISH_REVIEW_VALIDATION Audit paired prefixes and collect submission tables.
repoDir=fileparts(mfilename('fullpath'));
outputDir=fullfile(repoDir,'review_validation_results');
cases=readtable(fullfile(outputDir,'cases_all.csv'),'TextType','string');
checks=struct([]);
for mode=["global" "local"]
    for p=[23 24 25]
        if mode=="global",shortName='reproduction';else,shortName='local_range';end
        short=load(fullfile(outputDir,sprintf('%s_%s_p%i_seed7_T100.mat',shortName,mode,p)),'curves');
        long=load(fullfile(outputDir,sprintf('long_horizon_%s_p%i_seed7_T300.mat',mode,p)),'curves');
        relative=max(abs(short.curves.error-long.curves.error(1:100)) ./ max(1,abs(short.curves.error)));
        assert(relative<1e-12,'Paired short/long prefix differs: p%i,%s,%g',p,mode,relative);
        check=struct('range_mode',mode,'p',p,'maximum_relative_prefix_error_difference',relative,'passed',true);
        if isempty(checks),checks=check;else,checks(end+1)=check;end
    end
end
primary=cases(cases.experiment=="local_range",:);
if exist(fullfile(outputDir,'cases_additional.csv'),'file')
    extra=readtable(fullfile(outputDir,'cases_additional.csv'),'TextType','string');
    extra=extra(extra.experiment=="higher_precision",:);
    common=intersect(primary.Properties.VariableNames,extra.Properties.VariableNames,'stable');
    primary=[primary(:,common);extra(:,common)];
end
primary=sortrows(primary,'p');
writetable(primary,fullfile(outputDir,'local_resolution_sweep.csv'));
long=cases(cases.experiment=="long_horizon",:);
if exist(fullfile(outputDir,'cases_selected_precision.csv'),'file')
    extra=readtable(fullfile(outputDir,'cases_selected_precision.csv'),'TextType','string');
    extra=extra(extra.experiment=="long_horizon",:);
    common=intersect(long.Properties.VariableNames,extra.Properties.VariableNames,'stable');
    long=[long(:,common);extra(:,common)];
end
writetable(long,fullfile(outputDir,'long_horizon_checks.csv'));
ablation=cases(cases.experiment=="safeguard_off_checked",:);
if exist(fullfile(outputDir,'cases_ablation_checked.csv'),'file')
    ablation=readtable(fullfile(outputDir,'cases_ablation_checked.csv'),'TextType','string');
end
writetable(ablation,fullfile(outputDir,'safeguard_ablation.csv'));
if exist(fullfile(outputDir,'cases_additional.csv'),'file') && ...
        exist(fullfile(outputDir,'cases_selected_precision.csv'),'file')
    short=load(fullfile(outputDir,'higher_precision_local_p30_seed7_T100.mat'),'metrics','curves');
    extended=load(fullfile(outputDir,'long_horizon_local_p30_seed7_T300.mat'),'metrics','curves');
    assert(short.metrics.physical_prefix_sha256==extended.metrics.physical_prefix_sha256, ...
        'Actual state/process-noise/measurement-noise prefix differs.');
    assert(short.metrics.plant_prefix_sha256==extended.metrics.plant_prefix_sha256, ...
        'Actual A/C plant prefix differs.');
    physical=table(true,true,short.metrics.physical_prefix_sha256,short.metrics.plant_prefix_sha256, ...
        'VariableNames',{'physical_logs_prefix_identical','plant_matrices_prefix_identical', ...
        'physical_prefix_sha256','plant_prefix_sha256'});
    writetable(physical,fullfile(outputDir,'physical_prefix_checks.csv'));
    relative=max(abs(short.curves.error-extended.curves.error(1:100)) ./ max(1,abs(short.curves.error)));
    assert(relative<1e-12,'Fresh p30 short/long error prefix differs: %g',relative);
    checks(end+1)=struct('range_mode',"local",'p',30, ...
        'maximum_relative_prefix_error_difference',relative,'passed',true);
    fprintf('Independent hash audit: actual A/C and x/w/v prefixes match exactly.\n');
end
writetable(struct2table(checks),fullfile(outputDir,'paired_prefix_checks.csv'));
disp(struct2table(checks));
fprintf('All %i paired 100/300-step error prefixes agree within 1e-12.\n',numel(checks));
end
