%% load STL file into matlab
% stl = stlread('antennas/Dipole10cmT264.stl');
% stl = stlread('antennas/Dipole10cmT580.stl'); %ok
% stl = stlread('antennas/Dipole10cmT722.stl'); %god
% stl = stlread('antennas/Dipole10cmT744.stl'); %
% stl = stlread('antennas/Dipole10cmT904.stl'); %
% stl = stlread('antennas/Dipole10cmT924.stl'); %god
% stl = stlread('antennas/Dipole10cmT1060.stl'); %god
% stl = stlread('antennas/Dipole10cmT1104.stl'); %god
% stl = stlread('antennas/Dipole10cmT1458.stl'); %god 
% stl = stlread('antennas/Dipole10cmT1680.stl'); 
% stl = stlread('antennas/Dipole10cmT1922.stl'); %god
% stl = stlread('antennas/Dipole10cmT2312.stl'); %god
stl = stlread('antennas/Dipole10cmT2888.stl'); %god
% stl = stlread('antennas/Dipole10cmT3528.stl');
% stl = stlread('antennas/AntBinMesh2556.stl');
% stl = stlread('antennas/HalfAntT212.stl');

% stl =  stlread('antennas/AspecPrio/Dipole10cmT1680.stl');
%% faces and unique vertices
tic;
fprintf('\n')
disp('Removing duplicate points')
[p, t] = ArbitraryAntenna.RemoveEqualPoints(stl);
%% Parameters
% Controls amount of antenna
p1 = p;
p2 = p;
p3 = p;
p4 = p;
% p1(:,2) = p(:,2)-0.0515;
% p2(:,2) = p(:,2)+0.0515;
% p3(:,1) = p(:,1)-0.08;
% p4(:,1) = p(:,1)-0.05;
% p = [p1; p2];%; p3; p4];
% t = [t; t+length(p1)];% t+length(p1)+length(p2); t+length(p1)+length(p2)+length(p3)];
% % p(:,1) = p(:,1)+0.03;
% Should source be dipole, if 0 a plane wave propagating in +x direction used
UseDipole = 0;
DipolePoint = [0,0,0];
% If set to one use 81 sub triangles pr element, if 0 use 9
SubSubTri = 0;
% if 1 use fast (but more inacurate) MoM
vectorized = 1;
InTest = 0;
% Use subtriangles to calculate current
sub = 1;
% Emmision parameters and size of plottet area
normalize = 0;
PlotComp = 0;
xmin = -2; xmax = 2;
ymin = -2; ymax = 2;
zmin = -2; zmax = 2;
steps = 200;
PointArea = xmax^2/steps;
%% Visual check
figure(1)
hold on
plot3(p(:,1),p(:,2),p(:,3),'*')
axis image
toc;
%% Calculating dimensions of dipole
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(max(p));
radius = 0.003;
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
%% SubTri
tic;
fprintf('\n')
disp('Calculating areals for subtriangles')
[SubTri] = ArbitraryAntenna.SubTriangles(p, t, Center, SubSubTri);
toc;
%% Basis Function setup
tic;
fprintf('\n')
disp('Defining basis functions')
[EdgeList, Basis, BasisLA] = ArbitraryAntenna.BasisFunc(p, t, ConnectCell);
toc;
%% Evaluating Basis Functions
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
%% Calculating Dipole strength on antenna points
[Ei] = ArbitraryAntenna.PointSource(0, 1, 0, w, mu0, k, Center, DipolePoint, [0,1,0], PointArea);
%% Calculating Interface params
[Interfacesurf, iG] = ArbitraryAntenna.InterfaceCalc(0.002, ymin, ymax, zmin, zmax, steps);
%% MoM
tic;
fprintf('\n')
disp('MoM')
if InTest
    [Z, b, a] = ArbitraryAntenna.MoMIG(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k,  SubTri, 0, 1, 0, UseDipole, Ei, InterfaceSurf, 3.4);
end
if vectorized
    [Z, a, b ] = ArbitraryAntenna.MoMVectorized(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, 0, 1, 0, UseDipole, Ei);
else
    [Z, b, a] = ArbitraryAntenna.MoM(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k,  SubTri, 0, 1, 0, UseDipole, Ei);
end
toc;
%% Current calc in Triangle
tic;
fprintf('\n')
disp('Calculating Current')
[Jface] = ArbitraryAntenna.CurrentCalc(t, EdgeList, w, mu0, a, BasisLA, RhoP, RhoM, RhoP_, RhoM_, sub, Area);
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
h=fill3(xthree', ythree', zthree', C', 'EdgeColor', 'none'); %linear scale
colormap gray;
colorbar;
axis('equal');
rotate3d
%%
ArbitraryAntenna.AngularFarField(w, mu0, k, 12, Center, Jface, 200)
%% Calculating E
tic;
fprintf('\n')
disp('Calculating E-field')
[Exy, Exz, Ezy, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ...
    ArbitraryAntenna.EField(Center, k, Jface, xmin, xmax, ymin, ymax, zmin, zmax, steps);
toc;
if UseDipole
fprintf('\n')
disp('Setting up Dipole')
[ExyD, ExzD, EzyD] = ...
    ArbitraryAntenna.Dipole(DipolePoint, k, [0,1,0] , xmin, xmax, ymin, ymax, zmin, zmax, steps, PointArea);
toc;
Exy=Exy+ExyD; Exz=Exz+ExzD; Ezy=Ezy+EzyD;
end
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