function finalize_validation_data()
%FINALIZE_VALIDATION_DATA Collect every valid scenario once for submission.
repoDir=fileparts(mfilename('fullpath'));
outputDir=fullfile(repoDir,'review_validation_results');
files={'cases_ablation_checked.csv','cases_all.csv','cases_additional.csv', ...
    'cases_selected_precision.csv','cases_diagnostic_curves.csv'};
combined=table;
for j=1:numel(files)
    rows=readtable(fullfile(outputDir,files{j}),'TextType','string');
    rows=rows(rows.experiment~="safeguard_off",:);
    for field=["condition_pre_infinite_count" "condition_post_infinite_count"]
        if ~ismember(field,string(rows.Properties.VariableNames))
            rows.(field)=nan(height(rows),1);
        end
    end
    rows=rows(:,sort(rows.Properties.VariableNames));
    if isempty(combined),combined=rows;else,combined=[combined;rows];end
end
key=combined.experiment+"|"+string(combined.p)+"|"+string(combined.noise_seed)+ ...
    "|"+string(combined.samples)+"|"+combined.range_mode+"|"+string(combined.safeguard);
[~,first]=unique(key,'stable');
combined=combined(first,:);
assert(height(combined)==47,'Unexpected current scenario count: %i',height(combined));
assert(nnz(combined.experiment=="safeguard_off_checked")==3);
assert(all(combined.termination_step(combined.experiment=="safeguard_off_checked")==3));
assert(all(combined.all_error_values_finite(combined.safeguard==1)));
assert(all(combined.codec_roundtrip_verified));
writetable(combined,fullfile(outputDir,'validation_cases.csv'));
fprintf('Canonical validation_cases.csv: %i distinct scenarios, including three checked failed ablations.\n',height(combined));
end
