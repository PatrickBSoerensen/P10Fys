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
tic;
fprintf('\n')
disp('MoM')
[Z, b, J, a] = ArbitraryAntenna.MoM(p, t, EdgeList, BasisNumber, BasisLA, BasisArea, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, 0, 1, 0);
toc;
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
pJ = (p(EdgeList(:,1),:)-p(EdgeList(:,2),:))/2;
figure(2)
plot(pJ(:,2),abs(J(:,2)),'*');
title('Current in y, plottet at midpoint of edge')
%% Plot current Triangles
figure(3)
plot(Center(:,2),abs(Jt(:,2)),'*');
title('Current in y')
%% Calculating E
%x-min/max, z-min/max, y-min/max
tic;
fprintf('\n')
disp('Calculating E-field')
% [Eyx, Ezx, Eyz, x, y, z] = ArbitraryAntenna.EField(Center, w, k, mu0, J, -25, 325, -10, 10, -10, 10, 200, BasisArea, BasisLA, RhoP_, RhoM_, a, SubTri, t, EdgeList);

% [Eyx, Ezx, Eyz, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ArbitraryAntenna.EField(Center, w, k, mu0, J, -5, 5, -5, 5, -5, 5, 500, BasisArea, BasisLA, RhoP_, RhoM_, a, SubTri, t, EdgeList);

% [Exy, Exz, Eyz, xrange, yrange, zrange] = ArbitraryAntenna.EFieldAlt(pJ, w, k, mu0,...
%                 -5, 5, -5, 5, -5, 5, 500, BasisArea, RhoP_, RhoM_, a, SubTri);

[Exyx, Exyy, Exyz, xrange, yrange] = ArbitraryAntenna.EFieldXY(Center, pJ, w, k, mu0,...
                -5, 5, -5, 5, 50, BasisLA, RhoP_, RhoM_, a, SubTri, t, EdgeList);
toc;

Exyx = Exyx/max(max(Exyx));
Exyy = Exyy/max(max(Exyy));
Exyz = Exyz/max(max(Exyz));

% Exy = Exy/max(max(Exy));
% Ezx = Ezx/max(max(Ezx));
% Eyz = Eyz/max(max(Eyz));
%% Plotting E
figure(5)
pcolor(xrange, yrange, abs((Exyx)))
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('y');
title('yx plane');

figure(6)
pcolor(xrange, yrange, abs((Exyy)))
shading interp
colorbar
caxis([0 0.1])
xlabel('x');
ylabel('z');
title('zx plane');

figure(7)
pcolor(xrange, yrange, abs((Exyz)))
shading interp
colorbar
caxis([0 0.1])
xlabel('z');
ylabel('y');
title('yz plane');

% figure(4)
% pcolor(x, y, abs((Exy2)))
% shading interp
% colorbar
% % caxis([0 0.1])
% xlabel('x');
% ylabel('y');
% title('yx plane');
% 
% figure(5)
% pcolor(x, z, abs((Exz2)))
% shading interp
% colorbar
% % caxis([0 0.1])
% xlabel('x');
% ylabel('z');
% title('zx plane');
% 
% figure(6)
% pcolor(z, y, abs((Eyz2)))
% shading interp
% colorbar
% % caxis([0 0.1])
% xlabel('z');
% ylabel('y');
% title('yz plane');