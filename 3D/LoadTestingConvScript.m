%% load STL file into matlab
J =[];
center = [];
ExyCrossX = [];
ExzCrossZ =[];
ESC = [];
%3mm
% stl1 = stlread('antennas/Dipole10cmT264.stl');
% stl2 = stlread('antennas/Dipole10cmT580.stl'); %ok
% stl3 = stlread('antennas/Dipole10cmT722.stl'); %god
% stl4 = stlread('antennas/Dipole10cmT924.stl'); %god
% stl5 = stlread('antennas/Dipole10cmT1060.stl'); %god
% stl6 = stlread('antennas/Dipole10cmT1104.stl');
% stl7 = stlread('antennas/Dipole10cmT1922.stl'); %god
% Amount = 7;
%1mm
% stl1 = stlread('antennas/Dipole1mm/Dipole10cm552T1mm.stl');
% stl2 = stlread('antennas/Dipole1mm/Dipole10cm702T1mm.stl'); %ok
% stl3 = stlread('antennas/Dipole1mm/Dipole10cm900T1mm.stl'); %god
% stl4 = stlread('antennas/Dipole1mm/Dipole10cm1190T1mm.stl'); %god
% stl5 = stlread('antennas/Dipole1mm/Dipole10cm1296T1mm.stl'); %god
% stl6 = stlread('antennas/Dipole1mm/Dipole10cm1444T1mm.stl');
% stl7 = stlread('antennas/Dipole1mm/Dipole10cm1560T1mm.stl'); %god
% stl8 = stlread('antennas/Dipole1mm/Dipole10cm1936T1mm.stl'); %god
% stl9 = stlread('antennas/Dipole1mm/Dipole10cm2070T1mm.stl');
% stl10 = stlread('antennas/Dipole1mm/Dipole10cm2704T1mm.stl'); %god
% Amount = 10;
% 0.5mm
% stl1 = stlread('antennas/DipoleHalfmm/Dipole10CM576T.stl');
% stl2 = stlread('antennas/DipoleHalfmm/Dipole10CM700T.stl'); %ok
% stl3 = stlread('antennas/DipoleHalfmm/Dipole10CM836T.stl'); %god
% stl4 = stlread('antennas/DipoleHalfmm/Dipole10CM984T.stl'); %god
% stl5 = stlread('antennas/DipoleHalfmm/Dipole10CM1196T.stl'); %god
% stl6 = stlread('antennas/DipoleHalfmm/Dipole10CM1372T.stl');
% stl7 = stlread('antennas/DipoleHalfmm/Dipole10CM1560T.stl'); %god
% stl8 = stlread('antennas/DipoleHalfmm/Dipole10CM1760T.stl'); %god
% stl9 = stlread('antennas/DipoleHalfmm/Dipole10CM2040T.stl');
% stl10 = stlread('antennas/DipoleHalfmm/Dipole10CM4416T.stl'); %god
% Amount = 10;
%0.1mm
stl1 = stlread('antennas/Dipole0.1mm/Dipole10cm2502T01mm.stl');
stl2 = stlread('antennas/Dipole0.1mm/Dipole10cm3060T01mm.stl'); %ok
stl3 = stlread('antennas/Dipole0.1mm/Dipole10cm3696T01mm.stl'); %god
stl4 = stlread('antennas/Dipole0.1mm/Dipole10cm4392T01mm.stl'); %god
stl5 = stlread('antennas/Dipole0.1mm/Dipole10cm5200T01mm.stl'); %god
stl6 = stlread('antennas/Dipole0.1mm/Dipole10cm6020T01mm.stl');
stl7 = stlread('antennas/Dipole0.1mm/Dipole10cm6870T01mm.stl'); %god
stl8 = stlread('antennas/Dipole0.1mm/Dipole10cm7808T01mm.stl'); %god
stl9 = stlread('antennas/Dipole0.1mm/Dipole10cm8874T01mm.stl');
stl10 = stlread('antennas/Dipole0.1mm/Dipole10cm9936T01mm.stl'); %god
Amount = 10;
%% Parameters
% Controls amount of antenna
% p1 = p;
% p2 = p;
% p3 = p;
% p4 = p;
% p1(:,1) = p(:,1)-0.01;
% p2(:,1) = p(:,1)+0.01;
% p3(:,1) = p(:,1)-0.08;
% p4(:,1) = p(:,1)-0.05;
% p = [p1; p2];%; p3; p4];
% t = [t; t+length(p1)];% t+length(p1)+length(p2); t+length(p1)+length(p2)+length(p3)];
% p(:,1) = p(:,1)+0.03;
% Should source be dipole, if 0 a plane wave propagating in +x direction used
UseDipole = 0;

UseFeed = 1;
FeedPos = [0,0,0];

Yagi=0;
%If set to one use 81 sub triangles pr element, if 0 use 9
SubSubTri = 0;
sub =0;
% if 1 use fast (but more inacurate) MoM
vectorized = 1;
% Area of radiation
xmin = -2; xmax = 2;
ymin = -2; ymax = 2;
zmin = -2; zmax = 2;
steps = 200;
PointArea = xmax^2/steps;
% Reflector surface params
n = 3.9;
epsR = 11.68;
Reflector = 0;
Lift=0;

FileName= 'SergeyFeedPoint1mm';

%% Loop
for convloop=1:Amount
convloop
if convloop ==1
stl = stl1;
elseif convloop ==2
stl = stl2;
elseif convloop ==3
stl = stl3;
elseif convloop ==4
stl = stl4;
elseif convloop ==5
stl = stl5;
elseif convloop ==6
stl = stl6;
elseif convloop ==7
stl = stl7;
elseif convloop ==8
stl = stl8;
elseif convloop ==9
stl = stl9;
elseif convloop ==10
stl = stl10;
elseif convloop ==11
stl = stl11;
elseif convloop ==12
stl = stl12;
elseif convloop ==13
stl = stl13;
elseif convloop ==14
stl = stl14;
end
%% faces and unique vertices
tic;
fprintf('\n')
disp('Removing duplicate points')
[p, t] = ArbitraryAntenna.RemoveDuplicatePoints(stl);
%%
radiusdet = [1 1 1];
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(maxp);
radiusdet(maxaxis) = 0;
radiusdet = logical(radiusdet);
radius = sum(abs(maxp(radiusdet))+abs(minp(radiusdet)))/4;
Length = maxmaxp+abs(minp(maxaxis));

DipolePoint = [-Length,0,0];


OGSize = size(t);
OGSize = OGSize(1)*1.5;

xdist = radius+0;
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
%%
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

%% Basis Function setup
tic;
fprintf('\n')
disp('Defining basis functions')
[EdgeList, Basis, BasisLA] = ArbitraryAntenna.BasisFunc(p, t, ConnectCell);
toc;
%%
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

    Ei(:,1) = 0;
    Ei(:,2) = 1;
    Ei(:,3) = 0;

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
%%
[Esc, EscPhi, EscTheta] = ArbitraryAntenna.AngularFarField(w, mu0, k, 30, Center, Jface, 300);
close all
%% Calculating E
tic;
fprintf('\n')
disp('Calculating E-field')
[Exy, Exz, Ezy, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ...
    ArbitraryAntenna.EField(Center, w, mu0, k, Jface, xmin, xmax, ymin, ymax, zmin, zmax, steps, Area, Reflector, xdist, n, lambda);
toc;
if UseDipole
fprintf('\n')
disp('Setting up Dipole')
[ExyD, ExzD, EzyD] = ...
    ArbitraryAntenna.PointSourceEmmision(DipolePoint, k, [0,1,0] , xmin, xmax, ymin, ymax, zmin, zmax, steps);
toc;
Exy=Exy+ExyD; Exz=Exz+ExzD; Ezy=Ezy+EzyD;
end

sx = size(J);
sy = size(Jface);
a = max(sx(1),sy(1));
J = [[J;zeros(abs([a 0]-sx))],[Jface;zeros(abs([a,0]-sy))]];

sx = size(center);
sy = size(Center);
a = max(sx(1),sy(1));
center = [[center ;zeros(abs([a 0]-sx))],[Center;zeros(abs([a,0]-sy))]];
ExyCrossX= [ExyCrossX Exy(:,steps/2)];
ExzCrossZ = [ExzCrossZ  Exz(:,steps/2)];
ESC = [ESC Esc];

end
save(FileName, 'J', 'center', 'ExyCrossX', 'ExzCrossZ', 'ESC', 'UseDipole', 'SubSubTri', 'vectorized');
