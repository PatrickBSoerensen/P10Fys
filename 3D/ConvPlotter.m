P264=1;
P580=1;
P722=1;
P744=1;
P904=1;
P1104=1;
close all

load('ConvSlow.mat')
% load('ConvSlowSub.mat')
% load('ConvFast.mat')
% load('ConvFastSub.mat')

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
if P1104
    figure(1)
    plot(center(:,17),abs(J(:,17)),'*')
    figure(2)
    plot(abs(ExyCrossX(:,6)))
    figure(3)
    plot(abs(ExzCrossZ(:,6)))
    PlottetLabels = [PlottetLabels, '1104 T'];
end
figure(1)
legend(PlottetLabels)
figure(2)
legend(PlottetLabels)
figure(3)
legend(PlottetLabels)