%% load STL file into matlab
% stl = stlread('AntBinMesh.stl');
stl = stlread('Dipole10cm.stl');
% stl = stlread('Dipole10cmT1104.stl');
% stl = stlread('AntBinMesh2556.stl');
% stl = stlread('BinMeshHigh.stl');
% stl = stlread('HalfAntMany.stl');
%% faces and unique vertices
tic;
fprintf('\n')
disp('Removing duplicate points')
[p, t] = ArbitraryAntenna.RemoveEqualPoints(stl);
%% Visual check
figure(1)
plot3(p(:,1),p(:,2),p(:,3),'*')
axis image
toc;
%% Calculating dimensions of dipole
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(max(p));
radius = abs(min(min(p))+min(min(p)));
Length = (maxmaxp-minp(maxaxis));
%% constants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s
%2/5, 2/3, 2 
lambda=2/3*Length;
f=c/lambda;
w=2*pi*f;
k=w/c;
%% Connectivity list
tic;
fprintf('\n')
disp('Connectivity Cell')
ConnectCell = ArbitraryAntenna.Connectivity(p, t);
toc;
%% Calculating areas
tic;
fprintf('\n')
disp('Calculating areals for triangles')
[Area, Center] = ArbitraryAntenna.TriangleAreas(p, t);
toc;
tic;
fprintf('\n')
disp('Calculating areals for subtriangles')
[SubTri, SubTriArea] = ArbitraryAntenna.SubTriangles(p, t, Center);
toc;
%% Basis Function setup
tic;
fprintf('\n')
disp('Defining basis functions')
[EdgeList, Basis, BasisLA, BasisDeriv, BasisNumber, BasisArea] = ArbitraryAntenna.BasisFunc(p, t, ConnectCell);
toc;
tic;
fprintf('\n')
disp('Evaluating basis functions in center points')
[RhoP, RhoM, RhoP_, RhoM_] = ArbitraryAntenna.BasisEvalCenter(t, EdgeList, Basis, Center, SubTri);
toc;
%% pre analytic calculations
%Self Terms
tic;
fprintf('\n')
disp('Pre-Calculating self-coupling terms')
I2 = ArbitraryAntenna.SelfTerm(p, t);
toc;
%% MoM
pJ = (p(EdgeList(:,1),:)-p(EdgeList(:,2),:))/2;
tic;
fprintf('\n')
disp('MoM')
[ZN,aN, bN ] = ArbitraryAntenna.MoMLoopCut(p, t, EdgeList, BasisNumber, BasisLA, BasisArea, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, w, mu0, eps0, SubTri, 0, 1, 0);
toc;
[Z, b, J, a] = ArbitraryAntenna.MoM(p, t, EdgeList, BasisNumber, BasisLA, BasisArea, RhoP, RhoM, RhoP_, RhoM_, I2, Center, pJ, k,  SubTri, 0, 1, 0);
toc;
sum(sum(ZN==Z))
%% Current calc in center Triangle
Jt = zeros(length(t),3);
for i=1:length(t)
    Edge(1,:) = t(i,1:2);
    Edge(2,:) = t(i,2:3);
    Edge(3,1) = t(i,1);
    Edge(3,2) = t(i,3);
    [BasisNumberJ] = ArbitraryAntenna.EdgeNumbering(EdgeList, Edge);
    
    Jt(i,:) = sum(J(BasisNumberJ,:),1)/3;
end
%% Plot current Edges
figure(2)
plot(pJ(:,2),abs(J(:,2)),'*');
title('Current in y, plottet at midpoint of edge')
%% Plot current Triangles
figure(3)
plot(Center(:,2),abs(Jt(:,2)),'*');
title('Current in y')
%% Surf plot Current
[TrianglePlus, TriangleMinus] = ArbitraryAntenna.PMTri(t, EdgeList);
clear xthree; clear ythree; clear zthree;
JT = sqrt(Jt(:,1).^2+Jt(:,2).^2+Jt(:,3).^2);
J1 = sqrt(J(:,1).^2+J(:,2).^2+J(:,3).^2);

for n=1:length(t)
    i=[0 0 0];
    for m=1:length(EdgeList)
        IE=a(m)*BasisLA(m,2);
        if(TrianglePlus(m)==n)
            i=i+IE*RhoP(m,:)/(2*Area(TrianglePlus(m)));
        end
        if(TriangleMinus(m)==n)
            i=i+IE*RhoM(m,:)/(2*Area(TriangleMinus(m)));
        end
    end
    CurrentNorm(n)=abs(norm(i));
end

Jmax=max(CurrentNorm);
MaxCurrent=strcat(num2str(Jmax),'[A/m]')
CurrentNorm1=CurrentNorm/max(CurrentNorm);
for m=1:length(t)
    N=t(m,1:3);
    xthree(m,1:3) = p(N,1);
    ythree(m,1:3) = p(N,2);
    zthree(m,1:3) = p(N,3);
end
ourC=repmat(CurrentNorm1,3,1);

figure(5)
h=fill3(xthree', ythree', zthree', ourC); %linear scale
colormap gray;
colorbar;
axis('equal');
rotate3d

%% Calculating E
%x-min/max, z-min/max, y-min/max
normalize = 1;
PlotComp = 0;
tic;
fprintf('\n')
disp('Calculating E-field')
[Eyx, Ezx, Eyz, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ArbitraryAntenna.EField(Center, w, k, mu0, J, -5, 5, -5, 5, -5, 5, 100, BasisArea, BasisLA, RhoP_, RhoM_, a, SubTri, t, EdgeList);
toc;
if normalize
Exyx = Exyx/max(max(Exyx));
Exyy = Exyy/max(max(Exyy));
Exyz = Exyz/max(max(Exyz));

Exzx = Exzx/max(max(Exzx));
Exzy = Exzy/max(max(Exzy));
Exzz = Exzz/max(max(Exzz));

Eyzx = Eyzx/max(max(Eyzx));
Eyzy = Eyzy/max(max(Eyzy));
Eyzz = Eyzz/max(max(Eyzz));
end
%% Plotting E
if PlotComp
figure(4)
pcolor(x, y, abs(Exyx).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('y');
title('xy plane- x comp');

figure(5)
pcolor(x, y, abs(Exyy).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('y');
title('xy plane - y comp');

figure(6)
pcolor(x, y, abs(Exyz).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('y');
title('xy plane - z comp');
end

figure(7)
pcolor(x, y, sqrt(abs(Exyx).^2+abs(Exyy).^2+abs(Exyz).^2).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('y');
title('xy plane E-size');

if PlotComp
figure(8)
pcolor(x, z, abs(Exzx).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('z');
title('xz plane - x comp');

figure(9)
pcolor(x, z, abs(Exzy).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('z');
title('xz plane - y comp');

figure(10)
pcolor(x, z, abs(Exzz).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('z');
title('xz plane - z comp');
end

figure(11)
pcolor(x, z, sqrt(abs(Exzx).^2+abs(Exzy).^2+abs(Exzz).^2).')
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('z');
title('xz  plane E-size');

if PlotComp
figure(12)
pcolor(z, y, abs(Eyzx).')
shading interp
colorbar
caxis([0 0.1])
xlabel('z');
ylabel('y');
title('xzx plane');

figure(13)
pcolor(z, y, abs(Eyzy).')
shading interp
colorbar
caxis([0 0.1])
xlabel('z');
ylabel('y');
title('xzy plane');

figure(14)
pcolor(z, y, abs(Eyzz).')
shading interp
colorbar
caxis([0 0.1])
xlabel('z');
ylabel('y');
title('xzz plane');
end

figure(15)
pcolor(z, y, sqrt(abs(Eyzx).^2+abs(Eyzy).^2+abs(Eyzz).^2).')
shading interp
colorbar
caxis([0 0.1])
xlabel('z');
ylabel('y');
title('xz plane E-size');