P264=1;
P580=1;
P722=1;
P744=1;
P904=1;
P924=0;
P1060=0;
P1104=0;
P1458 =0;
P1680 =0;
P1922 =0;
P2312 =0;
P2888 =0;
P3528 =0;

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
load('conv/test/ConvSlowtestWave.mat')
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

figure(4)
title('Convergence of max(abs(J_y))')
xlabel('y')
ylabel('max(abs(J_{y}))')
hold on

if P264
    figure(1)
    plot(center(:,2),abs(J(:,2)),'*')
    figure(4)
    Jm = max(abs(J(:,2)));
    plot(264,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,1)))
    figure(3)
    plot(abs(ExzCrossZ(:,1)))
    PlottetLabels = [PlottetLabels, '1444 T'];
end
if P580
    figure(1)
    plot(center(:,5),abs(J(:,5)),'*')
    figure(4)
    [Jm, Ji] = max(abs(J(:,5)));
    plot(580,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,2)))
    figure(3)
    plot(abs(ExzCrossZ(:,2)))
    PlottetLabels = [PlottetLabels, '1634 T'];
end
if P722
    figure(1)
    plot(center(:,8),abs(J(:,8)),'*')
    figure(4)
    [Jm, Ji] = max(abs(J(:,8)));
    plot(722,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,3)))
    figure(3)
    plot(abs(ExzCrossZ(:,3)))
    PlottetLabels = [PlottetLabels, '1900 T' ];
end
if P744
    figure(1)
    plot(center(:,11),abs(J(:,11)),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,11)));
    plot(744,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,4)))
    figure(3)
    plot(abs(ExzCrossZ(:,4)))
    PlottetLabels = [PlottetLabels, '2280 T' ];
end
if P904
    figure(1)
    plot(center(:,14),abs(J(:,14)),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,14)));
    plot(904,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,5)))
    figure(3)
    plot(abs(ExzCrossZ(:,5)))
    PlottetLabels = [PlottetLabels, '2546 T'];
end
if P924
    figure(1)
    plot(center(:,17),abs(J(:,17))/max(abs(J(:,17))),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,17)));
    plot(924,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,5)))
    figure(3)
    plot(abs(ExzCrossZ(:,5)))
    PlottetLabels = [PlottetLabels, '924 T'];
end
if P1060
    figure(1)
    plot(center(:,20),abs(J(:,20)),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,20)));
    plot(1060,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1060 T'];
end
if P1104
    figure(1)
    plot(center(:,23),abs(J(:,23)),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,23)));
    plot(1104,abs(Jm),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1104 T'];
end
if P1458
    figure(1)
    plot(center(:,26),abs(J(:,26)),'*')
    figure(4)
    [Jm, Ji] = max(abs(J(:,26)));
    plot(1458,abs(Jm),'*')
    
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1458 T'];
end
if P1680
    figure(1)
    plot(center(:,29),abs(J(:,29)),'*')
    figure(4)
    [Jm, Ji] = max(abs(J(:,29)));
    plot(1680,abs(Jm),'*')
    
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1680 T'];
end
if P1922
    figure(1)
    plot(center(:,32),abs(J(:,32)),'*')
    figure(4)
    [Jm, Ji] = max(abs(J(:,32)));
    plot(1922,abs(Jm),'*')
    
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1922 T'];
end
if P2312
    figure(1)
    plot(center(:,35),abs(J(:,35)),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,35)));
    plot(2312,abs(Jm),'*')
    
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '2312 T'];
end
if P2888
    figure(1)
    plot(center(:,38),abs(J(:,38)),'*')
    figure(4)
    [Jm, Ji] = max(abs(J(:,38)));
    plot(2888,abs(Jm),'*')
    
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '2888 T'];
end
if P3528
    figure(1)
    plot(center(:,41),abs(J(:,41)),'*')
    
    figure(4)
    [Jm, Ji] = max(abs(J(:,41)));
    plot(3528,abs(Jm),'*')
    
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '3528 T'];
end

figure(1)
legend(PlottetLabels)
figure(2)
legend(PlottetLabels)
figure(3)
legend(PlottetLabels)
figure(4)
legend(PlottetLabels)