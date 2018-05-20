%% load STL file into matlab
stl = stlread('antennas/Dipole1mm/Dipole10cm702T1mm.stl'); %ok
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
% stl = stlread('antennas/Dipole10cmT2888.stl'); %god
% stl = stlread('antennas/Dipole10cmT3528.stl');
% stl = stlread('antennas/AntBinMesh2556.stl');
% stl = stlread('antennas/HalfAntT212.stl');
% stl =  stlread('antennas/AspecPrio/Dipole10cmT910.stl');
%% faces and unique vertices
tic;
fprintf('\n')
disp('Removing duplicate points')
[p, t] = ArbitraryAntenna.RemoveDuplicatePoints(stl);
%% Calculating dimensions of dipole
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(max(p));
radius = 0.003;
Length = (maxmaxp-minp(maxaxis));
%% Parameters
% Controls amount of antenna
p1 = p;
p1(:,2) = p(:,2)*1.038;
p2 = p;
p2(:,2) = p(:,2).*0.955;
p3 = p;
p3(:,2)  = p(:,2).*0.956;
p4 = p;
p4(:,2)  = p(:,2).*0.932;
p5 = p;
p5(:,2)  = p(:,2).*0.916;
p6 = p;
p6(:,2)  = p(:,2).*0.906;
p7 = p;
p7(:,2)  = p(:,2).*0.897;
p8 = p;
p8(:,2)  = p(:,2).*0.891;
p9 = p;
p9(:,2)  = p(:,2).*0.887;

p1(:,1) = p(:,1)-0.321;
p2(:,1) = p(:,1);
p3(:,1) = p(:,1)+0.135;
p4(:,1) = p(:,1)+0.378;
p5(:,1) = p(:,1)+0.729;
p6(:,1) = p(:,1)+1.170;
p7(:,1) = p(:,1)+1.674;
p8(:,1) = p(:,1)+2.241;
p9(:,1) = p(:,1)+2.874;
p = [p1; p2; p3; p4; p5; p6; p7; p8; p9];
t = [t; t+length(p1); t+length(p1)+length(p2);...
    t+length(p1)+length(p2)+length(p3);...
    t+length(p1)+length(p2)+length(p3)+length(p4);...
    t+length(p1)+length(p2)+length(p3)+length(p4)+length(p5);...
    t+length(p1)+length(p2)+length(p3)+length(p4)+length(p5)+length(p6);...
    t+length(p1)+length(p2)+length(p3)+length(p4)+length(p5)+length(p6)+length(p7);...
    t+length(p1)+length(p2)+length(p3)+length(p4)+length(p5)+length(p6)+length(p7)+length(p8)];
% Should source be dipole, if 0 a plane wave propagating in +x direction used
UseDipole = 0;
DipolePoint = [0,0,0];
UseFeed = 1;
FeedPos = [0,0,0];
Yagi=0;
OGSize = size(t);
OGSize = OGSize(1)*1.5;
% If set to one use 81 sub triangles pr element, if 0 use 9
SubSubTri = 0;
% if 1 use fast (but more inacurate) MoM
vectorized = 1;
InTest = 0;
% Emmision parameters and size of plottet area
normalize = 1;
PlotComp = 0;
xmin = -12; xmax = 12;
ymin = -5; ymax = 5;
zmin = -5; zmax = 5;
steps = 600;
PointArea = xmax^2/steps;
% Reflector surface params
n = 3.9;
epsR = 11.68;
Reflector = 0;
FromAnt=0;
xdist = radius/2+FromAnt;

%% Visual check
figure(1)
plot3(p(:,1),p(:,2),p(:,3),'*')
axis image
toc;
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
%% Lift
Lift = 0;
tic;
fprintf('\n')
disp('Lifting subtriangles and center points')
[Center, SubTri] = ArbitraryAntenna.CenterLift(Center, SubTri, radius, Lift);
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
%% Calculating Dipole strength on antenna points
clear Ei
if UseDipole
    [Ei] = ArbitraryAntenna.PointSource(w, mu0, k, Center, SubTri, sub, DipolePoint, [0,1,0]);
end
if UseFeed
    [Ei, v] = ArbitraryAntenna.VoltageFeed(t, p, Center, DipolePoint, 1, EdgeList, BasisLA, Yagi, OGSize);
end
if ~UseFeed && ~UseDipole
    Ei(:,1) = 0.*exp(1i*k.*(Center(:,2)));
    Ei(:,2) = 1.*exp(1i*k.*(Center(:,1)));
    Ei(:,3) = 0.*exp(1i*k.*(Center(:,1)));
end
if Reflector   
    [GIx, GIy, GIz] = ArbitraryAntenna.IDGreens(k, distx, Length, 2*radius, 0.005, 50, lambda, n, epsR, eps0, Center, SubTri);
 
    GI = GIx+GIy+GIz;
else
    GI = [];
end
%% MoM
tic;
fprintf('\n')
disp('MoM')
if vectorized
    [Z, a, b] = ArbitraryAntenna.MoMVectorized(w, mu0, t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Center, k, SubTri, Ei, Reflector, GI, n, eps0);
else
    [Z, a, b] = ArbitraryAntenna.MoM(w, mu0, t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Center, k,  SubTri, Ei, eps0);
end
toc;

if UseFeed 
    a=Z\(v)';
end
%% Current calc in Triangle
tic;
fprintf('\n')
disp('Calculating Current')
[Jface] = ArbitraryAntenna.CurrentCalc(t, EdgeList, a, BasisLA, RhoP, RhoM);
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
tic;
fprintf('\n')
disp('Angular far field calc')
ArbitraryAntenna.AngularFarField(w, mu0, k, 50, Center, Jface, 500);
toc
%% Calculating E
tic;
fprintf('\n')
disp('Calculating E-field')
[Exy, Exz, Ezy, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ...
    ArbitraryAntenna.EField(Center, w, mu0, k, Jface, xmin, xmax, ymin, ymax, zmin, zmax, steps, Area, Reflector, xdist, n, lambda);
toc;
%%
if UseDipole
fprintf('\n')
disp('Setting up Dipole')
[ExyD, ExzD, EzyD] = ...
    ArbitraryAntenna.PointSourceEmmision(DipolePoint, k, [0,1,0] , xmin, xmax, ymin, ymax, zmin, zmax, steps);
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
caxis([0 0.01])
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