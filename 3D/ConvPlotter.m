close all
%% First batch
% load('conv/ConvFirstBatch/ConvSlow.mat')
% load('conv/ConvFirstBatch/ConvSlowSub.mat')
% load('conv/ConvFirstBatch/ConvFast.mat')
% load('conv/ConvFirstBatch/ConvFastSub.mat')
%% Wave/Dipole initiatl batch
% load('conv/DipoleandWave/ConvSlowSubWave.mat')
%% Test batch
% load('ConvFastHalfTestWave.mat')
% load('ConvFastHalfTestDipole.mat')
% load('ConvFastSubHalfTestWave.mat')
% load('ConvFastSubHalfTestDipole.mat')
% load('ConvSlowHalfTestWave.mat')
% load('ConvSlowHalfTestDipole.mat')
% load('ConvSlowSubHalfTestWave.mat')
% load('ConvSlowSubHalfTestDi.mat')
% load('ConvSlowAspecWave.mat')
%% test low aspect ratio i.e. more uniform across antenna
% load('conv/test/ConvSlowtestWave.mat')
% load('conv/test/ConvSlowtestWaveFirstFive.mat')
% load('conv/test/ConvSlowtestDipole.mat')
% load('conv/test/ConvSlowSubtestWave.mat')
% load('conv/test/ConvSlowSubtestDipole.mat')
%% Tomatchtest Longer center triangles
% load('conv/Tomatchtest/ConvSlowWave.mat')
% load('conv/Tomatchtest/ConvSlowDipole.mat')
% load('conv/Tomatchtest/ConvSlowSubWave.mat')
% load('conv/Tomatchtest/ConvSlowSubDipole.mat')
%%
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
%%
figure(1)
legend(PlottetLabels)
figure(2)
legend(PlottetLabels)
figure(3)
legend(PlottetLabels)

saveas(figure(1), 'JyWave.jpg')
saveas(figure(2), 'EyxWave.jpg')


