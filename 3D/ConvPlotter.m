close all
PlottetLabels = {};
normalize = 0;

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
Amount = size(ExyCrossX);
start=1;
Minus=5;
Jcount = 2*start;
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
    plot(linspace(-2,2,200),abs(ExyCrossX(:,i))/max(abs(ExyCrossX(:,i))), 'linewidth', 1.5)
   
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
        plot(TriCent,abs(TriAmount),'.', 'markersize', 15)
        figure(2)
        plot(linspace(-2,2,200),abs(ExyCrossX(:,i)), 'linewidth', 1.5)
    
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
print(h,'JyPoint1mmPoint','-dpdf','-r0')
h=figure(2);
legend(PlottetLabels, 'FontSize', 11, 'box', 'off')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'EyxPoint1mmPoint','-dpdf','-r0')




