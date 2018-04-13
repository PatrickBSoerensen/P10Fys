%% load STL file into matlab
stl = stlread('AntBinMesh.stl');
%% faces and unique vertices
tic;
[p, t] = ArbitraryAntenna.RemoveEqualPoints(stl);
t = sort(t,2);
toc;
%% Calculating dimensions of dipole
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(max(p));
radius = abs(min(min(p))+min(min(p)));
length = (maxmaxp+minp(maxaxis));
%% konstants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;

lambda=2/3*length;
f=c/lambda;
w=2*pi*f;
k=w/c;
%% Gibson connectivity list
tic;
ConnectCell = ArbitraryAntenna.GibsonConnect(p, t);
toc;
%% Calculating areas for Simplex Coords
tic;
[A, Atot, Center] = ArbitraryAntenna.TriangleAreas(p, t);
toc;
tic;
[SubTri, Integral] = ArbitraryAntenna.SubTriangles(p, t, Center);
toc;
%% Basis Function setup
tic;
[EdgeList, Basis, BasisLA, BasisDeriv, BasisNumber] = ArbitraryAntenna.BasisFunc(p, t, ConnectCell);
toc;
tic;
[RhoP, RhoM, RhoP_, RhoM_] = ArbitraryAntenna.BasisEvalCenter(t, EdgeList, Basis, Center, SubTri);
toc;
%% pre analytic calculations
%Self Terms
tic;
I2 = ArbitraryAntenna.SelfTerm(p, t);
toc;
%MoM
tic;
[Z, b, J] = ArbitraryAntenna.MoM(p, t, EdgeList, BasisNumber, BasisLA, A, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri);
toc;
%% Calculating E
tic;
[Eyx, Ezx, Eyz, size] = ArbitraryAntenna.EField(Center, w, k, mu0, J, -1, 1, -1, 1, -1, 1, 100, A);
toc;

Eyx = Eyx/max(max(Eyx));
Ezx = Ezx/max(max(Ezx));
Eyz = Eyz/max(max(Eyz));
%% Plotting E
figure(1)
pcolor(size, size, abs(real(Eyx.')))
shading interp
colorbar
caxis([0 0.3])
title('yx plane');

figure(2)
pcolor(size, size, abs(real(Ezx.')))
shading interp
colorbar
% caxis([0 0.3])
title('zx plane');

figure(3)
pcolor(size, size, abs(real(Eyz.')))
shading interp
colorbar
% caxis([0 0.3])
title('yz plane');
