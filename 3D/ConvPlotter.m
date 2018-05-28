close all
PlottetLabels = {};
normalize = 0;
AddFile = '3mmWave';
Amount = size(EzyCrossZ);
start=1;
Jcount = 2*start;
Minus=0;

% fonts properties
    iFontSize       = 17;
    strFontName     = 'Helvetica';      % [Times | Courier | ] TODO complete the list
    strInterpreter  = 'latex';          % [{tex} | latex]
    fXLabelRotation = 0.0;
    fYLabelRotation = 90;

figure(1)
hold on
    xlabel( '\bf{y [cm]}', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    ylabel( '$\bf{|J_y|}$ (arb. units)', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    set(get(gca, 'XLabel'), 'Rotation', fXLabelRotation);
    set(get(gca, 'YLabel'), 'Rotation', fYLabelRotation);
    % in order to make matlab to do not "cut" latex-interpreted axes labels
    set(gca, 'Position', [0.15 0.15 0.75 0.75]);
    % general properties
    set(gca, 'FontName', strFontName, 'FontSize', 12);
%%
figure(2)
hold on
    xlabel( '\bf{x [m]}', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    ylabel( '$\bf{|E_{xy}|}$ (arb. units)',  'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    set(get(gca, 'XLabel'), 'Rotation', fXLabelRotation);
    set(get(gca, 'YLabel'), 'Rotation', fYLabelRotation);
    % in order to make matlab to do not "cut" latex-interpreted axes labels
    set(gca, 'Units', 'normalized', 'Position', [0.15 0.15 0.75 0.75]);
    set(gca, 'FontName', strFontName, 'FontSize', 12);
%%
figure(3)
hold on
    xlabel( '\bf{Theta}', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    ylabel( '$\bf{|E_{sc}|}$ (arb. units)', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    set(get(gca, 'XLabel'), 'Rotation', fXLabelRotation);
    set(get(gca, 'YLabel'), 'Rotation', fYLabelRotation);
    % in order to make matlab to do not "cut" latex-interpreted axes labels
    set(gca, 'Position', [0.15 0.15 0.75 0.75]);
    % general properties
    set(gca, 'FontName', strFontName, 'FontSize', 12);
    %%
figure(4)
hold on
    xlabel( '\bf{Theta}', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    ylabel( '$\bf{|E_{sc}|}$ (arb. units)', 'FontName', strFontName, 'FontSize', iFontSize, 'Interpreter', strInterpreter);
    set(get(gca, 'XLabel'), 'Rotation', fXLabelRotation);
    set(get(gca, 'YLabel'), 'Rotation', fYLabelRotation);
    % in order to make matlab to do not "cut" latex-interpreted axes labels
    set(gca, 'Position', [0.15 0.15 0.75 0.75]);
    % general properties
    set(gca, 'FontName', strFontName, 'FontSize', 12);
%%
if normalize
for i=start:Amount(2)-Minus
    TriAmount = J(:,Jcount);
    TriAmount(TriAmount==0)= [];
    TriCent = center(:,Jcount);
    TriCent(TriCent==0)= [];
    
    if isempty(TriAmount)
        Minus=Minus+1;
    else
    figure(1)
    plot(TriCent*100,abs(TriAmount)/max(abs(TriAmount)),'.', 'markersize', 11)
    figure(2)
    plot(linspace(-2,2,200),abs(EzyCrossZ(:,i))/max(abs(EzyCrossZ(:,i))), 'linewidth', 1.5)
    figure(3)
    plot(linspace(-pi,2*pi,200),1/2*(abs(ESCAng(i,:))/max(abs(ESCAng(i,:)))).^2*10^2, 'linewidth', 1.5)
    polarplot(theta, 1/2*abs(ESCAng(i,:))*10^2, 'linewidth', 1.5);
            
    figure(4)
    plot(linspace(-pi,2*pi,200),1/2*abs(EscRef(i,:))/max(abs(EscRef(i,:)))*10^2, 'linewidth', 1.5)
    
            
    TriAmount = int2str(length(TriAmount));
    PlottetLabels = [PlottetLabels, strcat(TriAmount, 'T')];
    end
    
    Jcount = Jcount+3;
end
else
for i=start:Amount(2)-Minus
        TriAmount = J(:,Jcount);
        TriAmount(TriAmount==0) = [];
        TriCent = center(:,Jcount);
        TriCent(TriCent==0)= [];
    
    if isempty(TriAmount)
        Minus=Minus+1;
    else
        figure(1)
        plot(TriCent,abs(TriAmount),'.', 'markersize', 8)
        figure(2)
        plot(linspace(-2,2,200),abs(EzyCrossZ(:,i)), 'linewidth', 1.5) 
        figure(3)
        plot(linspace(-pi,2*pi,200),1/2*10^2*abs(ESCAng(i,:)).^2, 'linewidth', 1.5)
        figure(4)
        plot(linspace(-pi,2*pi,200),1/2*10^2*abs(EscRef(i,:)), 'linewidth', 1.5)
    
        TriAmount = int2str(length(TriAmount));
        PlottetLabels = [PlottetLabels, strcat(TriAmount, 'T')];
    end
    Jcount = Jcount+3;
end
end
%%
h=figure(1);
legend(PlottetLabels, 'FontSize', 11, 'box', 'off')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('Jy',AddFile),'-dpdf','-r0')
h=figure(2);
legend(PlottetLabels, 'FontSize', 11, 'box', 'off')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('Eyx',AddFile),'-dpdf','-r0')
h=figure(3);
legend(PlottetLabels, 'FontSize', 11, 'box', 'off')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('Angular',AddFile),'-dpdf','-r0')
h=figure(4);
legend(PlottetLabels, 'FontSize', 11, 'box', 'off')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat('AngularRef', AddFile),'-dpdf','-r0')




