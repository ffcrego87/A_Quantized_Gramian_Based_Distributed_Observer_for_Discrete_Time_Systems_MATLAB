%% Post-processing and plots

fig_error = figure('Color','w','InvertHardcopy','off');
set(fig_error,'Units','centimeters','Position',[2 2 13 8.2]);
blue = [0.0000 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
green = [0.4660 0.6740 0.1880];
semilogy(necent,'LineWidth',1.5,'Color',blue);
hold on;
semilogy(nedist,'LineWidth',1.5,'Color',orange);
semilogy(nedistg,'LineWidth',1.5,'Color',green);

grid on;
box on;

xlabel('t','FontSize',14);
ylabel('Estimation error norm','FontSize',14);
show_error_legend = ~exist('paperFigurePrefix','var') || strcmp(paperFigurePrefix,'noquant');
if show_error_legend
    legend_handle = legend({'Centralized','Distributed','Gramian-based distributed'}, ...
        'Location','northeast','FontSize',9);
    style_legend(legend_handle);
end
set(gca,'FontSize',14);
style_paper_axes(gca);

% Save data for paper
save('necent.dat','necent','-ascii');
save('nedist.dat','nedist','-ascii');
save('nedistg.dat','nedistg','-ascii');

export_if_requested(fig_error,'fig_baseline');

%% Optional quantization plots and data export
do_quantization_plots = exist('use_quantization','var') && use_quantization;
do_saturation_plots = do_quantization_plots && exist('plot_saturation_diagnostics','var') && plot_saturation_diagnostics;

if do_quantization_plots && isfield(sim_output,'DeltaOmega_log') && isfield(sim_output,'DeltaMsg_log')
    
    DeltaOmega_log = sim_output.DeltaOmega_log;
    DeltaMsg_log = sim_output.DeltaMsg_log;

    fig_delta = figure('Color','w','InvertHardcopy','off');
    set(fig_delta,'Units','centimeters','Position',[2 2 13 8.2]);
    semilogy(DeltaOmega_log,'LineWidth',1.5,'Color',blue);
    hold on;
    semilogy(DeltaMsg_log,'LineWidth',1.5,'Color',orange);

    grid on;
    box on;

    xlabel('t','FontSize',14);
    ylabel('Quantization step','FontSize',14);
    legend_handle = legend({'\Delta_\Omega','\Delta_m'}, ...
        'Location','southwest','FontSize',9);
    style_legend(legend_handle);
    set(gca,'FontSize',14);
    style_paper_axes(gca);

    save('DeltaOmega.dat','DeltaOmega_log','-ascii');
    save('DeltaMsg.dat','DeltaMsg_log','-ascii');

    export_if_requested(fig_delta,'fig_delta');
end

if do_quantization_plots && isfield(sim_output,'XmaxOmega_log') && isfield(sim_output,'XmaxMsg_log')
    
    XmaxOmega_log = sim_output.XmaxOmega_log;
    XmaxMsg_log = sim_output.XmaxMsg_log;
    XmaxOmega_plot = positive_for_log_plot(XmaxOmega_log);
    XmaxMsg_plot = positive_for_log_plot(XmaxMsg_log);

    fig_xmax = figure('Color','w','InvertHardcopy','off');
    set(fig_xmax,'Units','centimeters','Position',[2 2 13 8.2]);
    yyaxis left;
    semilogy(XmaxOmega_plot,'LineWidth',1.5,'Color',blue);
    hold on;
    ylabel('X_{max}^{\Omega}','FontSize',14);
    yyaxis right;
    semilogy(XmaxMsg_plot,'LineWidth',1.5,'Color',orange);
    ylabel('X_{max}^{m}','FontSize',14);

    grid on;
    grid minor;
    box on;

    xlabel('t','FontSize',14);
    set(gca,'FontSize',14);
    style_paper_axes(gca);
    style_dual_axes(gca,blue,orange);

    save('XmaxOmega.dat','XmaxOmega_log','-ascii');
    save('XmaxMsg.dat','XmaxMsg_log','-ascii');

    export_if_requested(fig_xmax,'fig_xmax');
end

if do_saturation_plots && isfield(sim_output,'satOmega_log') && isfield(sim_output,'satMsg_log')
    
    satOmega_log = sim_output.satOmega_log;
    satMsg_log = sim_output.satMsg_log;

    fig_sat = figure('Color','w','InvertHardcopy','off');
    set(fig_sat,'Units','centimeters','Position',[2 2 13 8.2]);
    plot(satOmega_log,'LineWidth',1.5,'Color',blue);
    hold on;
    plot(satMsg_log,'LineWidth',1.5,'Color',orange);

    grid on;
    box on;

    xlabel('t','FontSize',14);
    ylabel('Saturation count','FontSize',14);
    legend_handle = legend({'\Omega saturation','Message saturation'}, ...
        'Location','best','FontSize',9);
    style_legend(legend_handle);
    set(gca,'FontSize',14);
    style_paper_axes(gca);

    save('satOmega.dat','satOmega_log','-ascii');
    save('satMsg.dat','satMsg_log','-ascii');

    export_if_requested(fig_sat,'fig_sat');
end

function values = positive_for_log_plot(values)
    safeguard_floor = 1e-9;
    values(values <= safeguard_floor) = NaN;
end

function style_paper_axes(ax)
    set(ax,'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'GridColor',[0.75 0.75 0.75], ...
        'MinorGridColor',[0.88 0.88 0.88], ...
        'GridAlpha',0.45, ...
        'MinorGridAlpha',0.35, ...
        'LineWidth',0.8);
end

function style_dual_axes(ax,left_color,right_color)
    ax.YAxis(1).Color = left_color;
    ax.YAxis(2).Color = right_color;
    ax.YAxis(1).Scale = 'log';
    ax.YAxis(2).Scale = 'log';
end

function style_legend(legend_handle)
    set(legend_handle,'Color','w', ...
        'TextColor','k', ...
        'EdgeColor',[0.35 0.35 0.35]);
end

function export_if_requested(fig_handle,base_name)
    if evalin('base','exist(''paperFigureDir'',''var'')') && evalin('base','exist(''paperFigurePrefix'',''var'')')
        paperFigureDir = evalin('base','paperFigureDir');
        paperFigurePrefix = evalin('base','paperFigurePrefix');
        if ~exist(paperFigureDir,'dir')
            mkdir(paperFigureDir);
        end
        filename = fullfile(paperFigureDir,[base_name '_' paperFigurePrefix '.pdf']);
        set(fig_handle,'PaperUnits','centimeters');
        fig_pos = get(fig_handle,'Position');
        set(fig_handle,'PaperSize',fig_pos(3:4));
        set(fig_handle,'PaperPosition',[0 0 fig_pos(3) fig_pos(4)]);
        set(fig_handle,'InvertHardcopy','off');
        print(fig_handle,filename,'-dpdf','-painters');
    end
end
