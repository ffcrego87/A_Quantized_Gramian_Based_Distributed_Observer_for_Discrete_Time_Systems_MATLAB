function make_review_validation_figure()
% Recreate the paired horizon/resolution figure from supplied case-level CSV.
root=fileparts(mfilename('fullpath'));
results=fullfile(root,'review_validation_results');
fig=figure('Visible','off','Color','w','Position',[100 100 900 730]);
colors=[0.16 0.47 0.71;0.88 0.53 0.15;0.21 0.57 0.42;0.55 0.35 0.66];
subplot(2,1,1);hold on;
pvals=[25 30];
colorIndex=[3 4];
for j=1:numel(pvals)
    c=readtable(fullfile(results,sprintf('long_horizon_local_p%i_seed7_T300.csv',pvals(j))));
    semilogy(c.time,c.normalized_error,'Color',colors(colorIndex(j),:),'LineWidth',1.4);
end
semilogy(c.time,c.reference_normalized_error,'--','Color',[.35 .35 .35],'LineWidth',1.3);
set(gca,'YScale','log');grid on;
xlabel('Sample');ylabel('Average error / max(1, state norm)');
legend('Gramian, p=25','Gramian, p=30', ...
    'Finite-precision DIF, p=30','NumColumns',2,'Location','northwest','FontSize',9);
subplot(2,1,2);hold on;
pvals=[25 28 30 44];
colorIndex=[3 2 4 1];
for j=1:numel(pvals)
    if pvals(j)==25,tag='local_range';elseif pvals(j)==44,tag='invertibility_margin';else,tag='higher_precision';end
    c=readtable(fullfile(results,sprintf('%s_local_p%i_seed7_T100.csv',tag,pvals(j))));
    semilogy(c.time,c.error,'Color',colors(colorIndex(j),:),'LineWidth',1.4);
end
c=readtable(fullfile(results,'higher_precision_local_p30_seed7_T100.csv'));
semilogy(c.time,c.reference_error,'--','Color',[.35 .35 .35],'LineWidth',1.3);
set(gca,'YScale','log');grid on;
xlabel('Sample');ylabel('Average error norm');
legend('Gramian, p=25','Gramian, p=28','Gramian, p=30','Gramian, p=44','Finite-precision DIF, p=30', ...
    'NumColumns',2,'Location','northeast','FontSize',9);
exportgraphics(fig,fullfile(results,'review_validation_figure.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(results,'review_validation_figure.png'),'Resolution',150);
close(fig);
end
