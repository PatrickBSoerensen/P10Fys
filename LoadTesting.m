%% load STL file into matlab
stl = stlread('AntBinMesh.stl');
% stl = stlread('AntBinMesh2.stl');
% stl = stlread('AntBinMesh2556.stl');
% stl = stlread('BinMeshHigh.stl');
% stl = stlread('HalfAntMany.stl');
%% faces and unique vertices
tic;
fprintf('\n')
disp('Removing duplicate points')
[p, t] = ArbitraryAntenna.RemoveEqualPoints(stl);
toc;
%% Calculating dimensions of dipole
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(max(p));
radius = abs(min(min(p))+min(min(p)));
length = (maxmaxp+minp(maxaxis));
%% constants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

lambda=2*length;
f=c/lambda;
w=2*pi*f;
k=w/c;
%% Connectivity list
tic;
fprintf('\n')
disp('Connectivity Cell')
ConnectCell = ArbitraryAntenna.Connectivity(p, t);
toc;
%% Calculating areas for Simplex Coords
tic;
fprintf('\n')
disp('Calculating areals for triangles')
[A, Atot, Center] = ArbitraryAntenna.TriangleAreas(p, t);
toc;
tic;
fprintf('\n')
disp('Calculating areals for subtriangles')
[SubTri] = ArbitraryAntenna.SubTriangles(p, t, Center);
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
tic;
fprintf('\n')
disp('MoM')
[Z, b, J, a] = ArbitraryAntenna.MoM(p, t, EdgeList, BasisNumber, BasisLA, BasisArea, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, 0, 1, 0);
toc;
%% Calculating E
%x-min/max, z-min/max, y-min/max
tic;
fprintf('\n')
disp('Calculating E-field')
% [Eyx, Ezx, Eyz, x, y, z] = ArbitraryAntenna.EField(Center, w, k, mu0, J, -25, 325, -10, 10, -10, 10, 200, BasisArea, BasisLA, RhoP_, RhoM_, a, SubTri, t, EdgeList);

[Eyx, Ezx, Eyz, x, y, z, Exy2, Exz2, Eyz2] = ArbitraryAntenna.EField(Center, w, k, mu0, J, -5, 5, -5, 5, -5, 5, 500, BasisArea, BasisLA, RhoP_, RhoM_, a, SubTri, t, EdgeList);

toc;

Eyx = Eyx/max(max(Eyx));
Ezx = Ezx/max(max(Ezx));
Eyz = Eyz/max(max(Eyz));
%% Plotting E
figure(1)
pcolor(x, y, abs((Eyx)))
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('y');
title('yx plane');

figure(2)
pcolor(x, z, abs((Ezx)))
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('z');
title('zx plane');

figure(3)
pcolor(z, y, abs((Eyz)))
shading interp
colorbar
caxis([0 0.1])
xlabel('z');
ylabel('y');
title('yz plane');

figure(4)
pcolor(x, y, abs((Exy2)))
shading interp
colorbar
% caxis([0 0.1])
xlabel('x');
ylabel('y');
title('yx plane');

figure(5)
pcolor(x, z, abs((Exz2)))
shading interp
colorbar
% caxis([0 0.1])
xlabel('x');
ylabel('z');
title('zx plane');

figure(6)
pcolor(z, y, abs((Eyz2)))
shading interp
colorbar
% caxis([0 0.1])
xlabel('z');
ylabel('y');
title('yz plane');