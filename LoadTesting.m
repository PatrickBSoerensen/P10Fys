%% load STL file into matlab
% stl = stlread('AntBinMesh.stl');
% stl = stlread('Dipole10cm.stl');
% stl = stlread('Dipole10cmT580.stl');
stl = stlread('Dipole10cmT1104.stl');
% stl = stlread('AntBinMesh2556.stl');
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
lambda=2*Length;
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
[ZN,aN, bN ] = ArbitraryAntenna.MoMLoopCut(t, EdgeList, BasisNumber, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, 0, 1, 0);
toc;
%%
% [Z, b, a] = ArbitraryAntenna.MoM(p, t, EdgeList, BasisNumber, BasisLA, Area, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k,  SubTri, 0, 1, 0);
% toc;
%% Current calc in center Triangle

sub = 0;
Dipole = 0;
tic;
fprintf('\n')
disp('Calculating Current')
% [Jface] = ArbitraryAntenna.CurrentCalc(t, p, EdgeList, w, mu0, a, BasisLA, RhoP, RhoM, RhoP_, RhoM_, sub, Dipole);

[Jface] = ArbitraryAntenna.CurrentCalc(t, p, EdgeList, w, mu0, aN, BasisLA, RhoP, RhoM, RhoP_, RhoM_, sub, Dipole);
toc;
%% Surf plot Current
JfaceSize = sqrt(sum(Jface.^2,2));
xthree = zeros(size(t)); ythree = zeros(size(t)); zthree = zeros(size(t));
Jmax=max(JfaceSize);
CurrentNorm1=JfaceSize/max(JfaceSize);
for m=1:length(t)
    N=t(m,1:3);
    xthree(m,1:3) = p(N,1);
    ythree(m,1:3) = p(N,2);
    zthree(m,1:3) = p(N,3);
end
C=repmat(CurrentNorm1,1,3);
figure(3)
h=fill3(xthree', ythree', zthree', Jface(:,3)'); %linear scale
colormap gray;
colorbar;
axis('equal');
rotate3d

%% Calculating E
%x-min/max, z-min/max, y-min/max
normalize = 1;
PlotComp = 1;
tic;
fprintf('\n')
disp('Calculating E-field')
[Exy, Exz, Ezy, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ...
    ArbitraryAntenna.EField(Center, w, k, mu0, Jface, -2, 2, -2, 2, -2, 2, 100, Area);
toc;

if normalize
Exyx = Exyx/max(max(Exy));
Exyy = Exyy/max(max(Exy));
Exyz = Exyz/max(max(Exy));

Exzx = Exzx/max(max(Exz));
Exzy = Exzy/max(max(Exz));
Exzz = Exzz/max(max(Exz));

Eyzx = Eyzx/max(max(Ezy));
Eyzy = Eyzy/max(max(Ezy));
Eyzz = Eyzz/max(max(Ezy));

Exy = Exy/max(max(Exy));
Exz = Exz/max(max(Exz));
Ezy = Ezy/max(max(Ezy));
end
%% Plotting E
if PlotComp
figure(4)
pcolor(x, y, real(Exyx).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('y');
title('xy plane- x comp');

figure(5)
pcolor(x, y, real(Exyy).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('y');
title('xy plane - y comp');

figure(6)
pcolor(x, y, real(Exyz).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('y');
title('xy plane - z comp');
end

figure(7)
pcolor(x, y, abs(Exy)')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('y');
title('xy plane E-size');
%%
if PlotComp
figure(8)
pcolor(x, z, real(Exzx).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('z');
title('xz plane - x comp');

figure(9)
pcolor(x, z, real(Exzy).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('z');
title('xz plane - y comp');

figure(10)
pcolor(x, z, real(Exzz).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('z');
title('xz plane - z comp');
end

figure(11)
pcolor(x, z, abs(Exz).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('x');
ylabel('z');
title('xz  plane E-size');
%% yz
if PlotComp
figure(12)
pcolor(z, y, real(Eyzx).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('z');
ylabel('y');
title('yz plane x-comp');

figure(13)
pcolor(z, y, real(Eyzy).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('z');
ylabel('y');
title('yz plane y-comp');

figure(14)
pcolor(z, y, real(Eyzz).')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('z');
ylabel('y');
title('yz plane z-comp');
end

figure(15)
pcolor(z, y, abs(Ezy)')
shading interp
colorbar
if normalize
caxis([0 0.1])
end
xlabel('z');
ylabel('y');
title('yz plane E-size');