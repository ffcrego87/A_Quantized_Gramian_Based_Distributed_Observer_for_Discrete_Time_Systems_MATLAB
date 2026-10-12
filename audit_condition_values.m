function audit_condition_values()
%AUDIT_CONDITION_VALUES Distinguish singular (Inf) from largest finite values.
repoDir=fileparts(mfilename('fullpath'));
records=struct([]);
for p=20:26
    stored=load(fullfile(repoDir,'diagnostics_results',sprintf('case_p%02i.mat',p)),'sim_output');
    pre=stored.sim_output.condPreAgent_log(2:end,:);
    post=stored.sim_output.condPostAgent_log(2:end,:);
    finitePre=pre(isfinite(pre));finitePost=post(isfinite(post));
    r=struct('p',p,'pre_infinite_count',nnz(isinf(pre)), ...
        'post_infinite_count',nnz(isinf(post)), ...
        'largest_finite_pre',max(finitePre),'largest_finite_post',max(finitePost), ...
        'true_maximum_pre',max(pre(:)),'true_maximum_post',max(post(:)));
    if isempty(records),records=r;else,records(end+1)=r;end
end
writetable(struct2table(records),fullfile(repoDir,'review_validation_results','condition_infinite_audit.csv'));
disp(struct2table(records));
end
