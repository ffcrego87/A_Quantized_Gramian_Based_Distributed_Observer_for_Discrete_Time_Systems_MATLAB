%% Post-processing and plots

figure;
semilogy(necent,'LineWidth',1.5);
hold on;
semilogy(nedist,'LineWidth',1.5);
semilogy(nedistg,'LineWidth',1.5);

grid on;
box on;

xlabel('t','FontSize',14);
ylabel('Estimation error norm','FontSize',14);
legend({'Centralized','Distributed','Gramian-based distributed'}, ...
    'Location','northeast','FontSize',12);
set(gca,'FontSize',14);

% Save data for paper
save('necent.dat','necent','-ascii');
save('nedist.dat','nedist','-ascii');
save('nedistg.dat','nedistg','-ascii');

%% Optional quantization plots and data export
if isfield(sim_output,'DeltaOmega_log') && isfield(sim_output,'DeltaMsg_log')
    
    DeltaOmega_log = sim_output.DeltaOmega_log;
    DeltaMsg_log = sim_output.DeltaMsg_log;

    figure;
    semilogy(DeltaOmega_log,'LineWidth',1.5);
    hold on;
    semilogy(DeltaMsg_log,'LineWidth',1.5);

    grid on;
    box on;

    xlabel('t','FontSize',14);
    ylabel('Quantization step','FontSize',14);
    legend({'\Delta_\Omega','\Delta_m'}, ...
        'Location','northeast','FontSize',12);
    set(gca,'FontSize',14);

    save('DeltaOmega.dat','DeltaOmega_log','-ascii');
    save('DeltaMsg.dat','DeltaMsg_log','-ascii');
end

if isfield(sim_output,'XmaxOmega_log') && isfield(sim_output,'XmaxMsg_log')
    
    XmaxOmega_log = sim_output.XmaxOmega_log;
    XmaxMsg_log = sim_output.XmaxMsg_log;

    figure;
    plot(XmaxOmega_log,'LineWidth',1.5);
    hold on;
    plot(XmaxMsg_log,'LineWidth',1.5);

    grid on;
    box on;

    xlabel('t','FontSize',14);
    ylabel('Adaptive quantization range','FontSize',14);
    legend({'Xmax_\Omega','Xmax_m'}, ...
        'Location','northeast','FontSize',12);
    set(gca,'FontSize',14);

    save('XmaxOmega.dat','XmaxOmega_log','-ascii');
    save('XmaxMsg.dat','XmaxMsg_log','-ascii');
end

if isfield(sim_output,'satOmega_log') && isfield(sim_output,'satMsg_log')
    
    satOmega_log = sim_output.satOmega_log;
    satMsg_log = sim_output.satMsg_log;

    figure;
    plot(satOmega_log,'LineWidth',1.5);
    hold on;
    plot(satMsg_log,'LineWidth',1.5);

    grid on;
    box on;

    xlabel('t','FontSize',14);
    ylabel('Saturation count','FontSize',14);
    legend({'\Omega saturation','Message saturation'}, ...
        'Location','northeast','FontSize',12);
    set(gca,'FontSize',14);

    save('satOmega.dat','satOmega_log','-ascii');
    save('satMsg.dat','satMsg_log','-ascii');
end