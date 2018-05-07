P264=0;
P580=0;
P722=1;
P744=0;
P904=0;
P924=1;
P1060=0;
P1104=1;
P1458 =0;
P1680 =0;
P1922 =0;
P2312 =0;
P2888 =0;
P3528 =0;

close all
%% First batch
% load('ConvSlow.mat')
% load('ConvSlowSub.mat')
% load('ConvFast.mat')
% load('ConvFastSub.mat')
%% Wave/Dipole initiatl batch
% load('ConvSlowSubWave.mat')
%% Test batch
% load('ConvFastHalfTestWave.mat')
% load('ConvFastHalfTestDipole.mat')
% load('ConvFastSubHalfTestWave.mat')
% load('ConvFastSubHalfTestDipole.mat')
load('ConvSlowHalfTestWave.mat')
% load('ConvSlowHalfTestDipole.mat')
% load('ConvSlowSubHalfTestWave.mat')
% load('ConvSlowSubHalfTestDi.mat')


figure(1)
title('Convergence of J_y')
xlabel('y')
ylabel('Size of J_y')
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

if P264
    figure(1)
    plot(center(:,2),abs(J(:,2)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,1)))
    figure(3)
    plot(abs(ExzCrossZ(:,1)))
    PlottetLabels = [PlottetLabels, '264 T'];
end
if P580
    figure(1)
    plot(center(:,5),abs(J(:,5)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,2)))
    figure(3)
    plot(abs(ExzCrossZ(:,2)))
    PlottetLabels = [PlottetLabels, '580 T'];
end
if P722
    figure(1)
    plot(center(:,8),abs(J(:,8)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,3)))
    figure(3)
    plot(abs(ExzCrossZ(:,3)))
    PlottetLabels = [PlottetLabels, '722 T' ];
end
if P744
    figure(1)
    plot(center(:,11),abs(J(:,11)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,4)))
    figure(3)
    plot(abs(ExzCrossZ(:,4)))
    PlottetLabels = [PlottetLabels, '744 T' ];
end
if P904
    figure(1)
    plot(center(:,14),abs(J(:,14)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,5)))
    figure(3)
    plot(abs(ExzCrossZ(:,5)))
    PlottetLabels = [PlottetLabels, '904 T'];
end
if P924
    figure(1)
    plot(center(:,17),abs(J(:,17)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,5)))
    figure(3)
    plot(abs(ExzCrossZ(:,5)))
    PlottetLabels = [PlottetLabels, '924 T'];
end
if P1060
    figure(1)
    plot(center(:,20),abs(J(:,20)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1060 T'];
end
if P1104
    figure(1)
    plot(center(:,23),abs(J(:,23)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1104 T'];
end
if P1458
    figure(1)
    plot(center(:,26),abs(J(:,26)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1458 T'];
end
if P1680
    figure(1)
    plot(center(:,26),abs(J(:,26)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1680 T'];
end
if P1922
    figure(1)
    plot(center(:,29),abs(J(:,29)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1922 T'];
end
if P2312
    figure(1)
    plot(center(:,32),abs(J(:,32)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '2312 T'];
end
if P2888
    figure(1)
    plot(center(:,35),abs(J(:,35)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '2888 T'];
end
if P3528
    figure(1)
    plot(center(:,38),abs(J(:,38)),'*')
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