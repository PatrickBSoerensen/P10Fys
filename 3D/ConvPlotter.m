% close all
normalize = 0;

figure(1)
title('Convergence of J_y')
xlabel('y')
ylabel('abs(J_y)')
PlottetLabels = {};
hold on

figure(2)
title('Convergence of E_{xy}')
xlabel('x')
ylabel('Size of E_{xy}')
hold on

figure(3)
title('Convergence of E_{xz}')
xlabel('z')
ylabel('Size of E_{xz}')
hold on

Amount = size(ExyCrossX);
Jcount = 2;
if normalize
for i=1:Amount(2)
    figure(1)
    plot(center(:,Jcount),abs(J(:,Jcount))/max(abs(J(:,Jcount))),'*')
    figure(2)
    plot(abs(ExyCrossX(:,i))/max(abs(ExyCrossX(:,i))))
    figure(3)
    plot(abs(ExzCrossZ(:,i))/max(abs(ExzCrossZ(:,i))))
    TriAmount = J(:,Jcount);
    TriAmount(TriAmount==0) = [];
    TriAmount = int2str(length(TriAmount));
    PlottetLabels = [PlottetLabels, strcat(TriAmount, 'T')];
    Jcount = Jcount+3;
end
else
for i=1:Amount(2)
    figure(1)
    plot(center(:,Jcount),abs(J(:,Jcount)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,i)))
    figure(3)
    plot(abs(ExzCrossZ(:,i)))
    TriAmount = J(:,Jcount);
    TriAmount(TriAmount==0) = [];
    TriAmount = int2str(length(TriAmount));
    PlottetLabels = [PlottetLabels, strcat(TriAmount, 'T')];
    Jcount = Jcount+3;
end
end
%%
figure(1)
legend(PlottetLabels)
figure(2)
legend(PlottetLabels)
figure(3)
legend(PlottetLabels)

% saveas(figure(1), 'JyWave.jpg')
% saveas(figure(2), 'EyxWave.jpg')


