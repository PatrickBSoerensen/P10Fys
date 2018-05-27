%% load STL file into matlab
% stl = stlread('antennas/Dipole1mm/Dipole10cm552T1mm.stl');
% stl = stlread('antennas/Dipole1mm/Dipole10cm702T1mm.stl'); %ok
% stl = stlread('antennas/Dipole1mm/Dipole10cm900T1mm.stl'); %god

% stl = stlread('antennas/AspecPrio/Dipole10cmT648.stl'); 

% stl = stlread('antennas/AspecPrio/Dipole10cmT910.stl'); 

% stl = stlread('antennas/test/720.stl');
% stl = stlread('antennas/Dipole10cmT180.stl');
% stl = stlread('antennas/Dipole10cmT264.stl');
stl = stlread('antennas/Dipole10cmT580.stl'); %ok
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
% stl = stlread('antennas/Dipole10cmT3528.stl');4
% stl = stlread('antennas/Dipole10cmT2888.stl'); %god
% stl = stlread('antennas/Dipole10cmT3528.stl');
% stl = stlread('antennas/AntBinMesh2556.stl');
% stl =  stlread('antennas/AspecPrio/Dipole10cmT910.stl');
%% faces and unique vertices
tic;
fprintf('\n')
disp('Removing duplicate points')
[p, t] = ArbitraryAntenna.RemoveDuplicatePoints(stl);
% [p1, t1] = ArbitraryAntenna.RemoveDuplicatePoints(stl1);
% prem = p(:,2) <= 0;
% p1rem = p1(:,2) >= 0;
% p(prem,:) = [];
% p1(p1rem,:) = [];
% trem = t>= length(p);
% t1rem = t1 >= length(p1);
% t(trem) = [];
% t1(t1rem) = [];
% t1 = t1+length(t);
% p = [p; p1];
% t = [t; t1];

OGSize = max(max(t));
p(:,1) = p(:,1);
p(:,2) = p(:,2);
p(:,3) = p(:,3);
%% Calculating dimensions of dipole
radiusdet = [1 1 1];
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(maxp);
radiusdet(maxaxis) = 0;
radiusdet = logical(radiusdet);
radius = sum(abs(maxp(radiusdet))+abs(minp(radiusdet)))/4;
Length = maxmaxp+abs(minp(maxaxis));
%% Parameters

AntFromReflector = 1*radius; %Should basically always be one radius
p(:,3) = p(:,3)+AntFromReflector;
InterAntDist = 4*radius;
% Controls amount of antenna
% p(:,1) = p(:,1)+0.03;
p1 = p;
p1(:,1) = p1(:,1)-InterAntDist;
p2 = p;
p2(:,1) = p2(:,1)+InterAntDist;
p(:,2) = p(:,2)*1;
% p = [ p1; p2];p;
% t = [t; t+length(p1) ];%t+length(p1)+length(p1);
% Should source be dipole, if 0 a plane wave propagating in +x direction used
UseDipole = 0;
DipolePoint = [0,0,0];
UseFeed = 0;
FeedPos = [0,0,0];

Yagi=0;
% If set to one use 81 sub triangles pr element, if 0 use 9
SubSubTri = 0;
sub = 0;
% if 1 use fast MoM
vectorized = 1;
% Emmision parameters and size of plottet area
normalize = 1;
PlotComp = 0;
xmin = -2; xmax = 2;
ymin = -2; ymax = 2;
zmin = -2; zmax = 2;
steps = 200;
PointArea = xmax^2/steps;
% Reflector surface params
n = 3.9;
epsR = 11.68;
Reflector = 0;
ReflectorZ = -0;%radius+FromAnt;
% Determines if points should be lifted to surf of antenna, this is semi
% hardcoded to a predetermined structure, if in doubt set to 0
Lift = 0;
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
tic;
fprintf('\n')
disp('Lifting subtriangles and center points')
[Center, SubTri] = ArbitraryAntenna.CenterLift(Center, SubTri, radius, Lift);
toc;
%% Basis Function setup
tic;
fprintf('\n')
disp('Defining basis functions')
[EdgeList, Basis, BasisLA, BasisCoord] = ArbitraryAntenna.BasisFunc(p, t, ConnectCell);
toc;
%% Evaluating Basis Functions
tic;
fprintf('\n')
disp('Evaluating basis functions in center points')
[RhoP, RhoM, RhoP_, RhoM_] = ArbitraryAntenna.BasisEvalCenter(t, EdgeList, Basis, Center, SubTri);
toc;
%% Calculating Dipole strength on antenna points
tic;
fprintf('\n')
disp('Calculating IncidentField and ID Greens if Reflector')
clear Ei
RefCoef = (1-n)/(1+n);
if Reflector
    [GIxx, GIxy, GIxz, GIyx, GIyy, GIyz, GIzx, GIzy, GIzz] = ArbitraryAntenna.IDGreens(k, ReflectorZ, Length, 2*radius, 0.003, 15, lambda, n, epsR, eps0, Center, SubTri);
%  0.003, 15 OR 0.002, 20
%     GI = GIx+GIy+GIz;
else
    GIxx=0; GIxy=0; GIxz=0; GIyx=0; GIyy=0; GIyz=0; GIzx=0; GIzy=0; GIzz=0;
    GI = [];
    GIx=0; GIy=0; GIz=0;
end
if UseDipole
    [Ei] = ArbitraryAntenna.PointSource(w, mu0, k, Center, SubTri, sub, DipolePoint, [0,1,0]);
end
if UseFeed
    [Ei, v] = ArbitraryAntenna.VoltageFeed(t, p, Center, DipolePoint, 1, EdgeList, BasisLA, Yagi, OGSize);
end
if ~UseFeed && ~UseDipole
    Ei(:,1) = 0.*exp(1i*k.*(Center(:,3)));
    Ei(:,2) = 1.*exp(1i*k.*(Center(:,3)));
    Ei(:,3) = 0.*exp(1i*k.*(Center(:,2)));
end
toc;
% if Reflector
%     Ei(:,1) = Ei(:,1) + Ei(:,1).*RefCoef;        
%     Ei(:,2) = Ei(:,2) + Ei(:,2).*RefCoef;
%     Ei(:,3) = Ei(:,3) + Ei(:,3).*RefCoef;
% end
%% MoM
% tic;
fprintf('\n')
disp('MoM')
if vectorized
    [Z, a, b] = ArbitraryAntenna.MoMVectorized(w, mu0, t, p, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, BasisCoord, Center, k, SubTri, Ei, Reflector, GIxx, GIxy, GIxz, GIyx, GIyy, GIyz, GIzx, GIzy, GIzz, eps0);
    
%     [Z, a, b] = ArbitraryAntenna.MoMOtherVectorized(w, mu0, t, p, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, BasisCoord, Center, k, SubTri, Ei, Reflector, GIxx, GIxy, GIxz, GIyx, GIyy, GIyz, GIzx, GIzy, GIzz, eps0);

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
figure(2)
h=fill3(xthree', ythree', zthree', C.', 'EdgeColor', 'none'); %linear scale
colormap gray;
colorbar;
axis('equal');
rotate3d

%%
[EscRef] = ArbitraryAntenna.AngularFarFieldSurf(w, mu0, k, 10, Center, Jface, 4*steps, Area,  ReflectorZ, epsR, eps0, n);

[Esc, EscPhi, EscTheta] = ArbitraryAntenna.AngularFarField(w, mu0, k, 10, Center, Jface, 4*steps, Area);

ProduceAnErrorHere
%% Calculating E   
tic;
fprintf('\n')
disp('Calculating E-field')
[Exy, Exz, Ezy, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ...
    ArbitraryAntenna.EField(Center, w, mu0, k, Jface, xmin, xmax, ymin, ymax, zmin, zmax, steps, Area, Reflector, ReflectorZ, n, lambda);
toc;
if UseDipole
fprintf('\n')
disp('Setting up Dipole')
[ExyD, ExzD, EzyD] = ...
    ArbitraryAntenna.PointSourceEmmision(DipolePoint, k, [0,1,0] , xmin, xmax, ymin, ymax, zmin, zmax, steps);
toc;
Exy=Exy+ExyD; Exz=Exz+ExzD; Ezy=Ezy+EzyD;
end
%%
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