%% load STL file into matlab
J =[];
center = [];
ExyCrossX = [];
ExzCrossZ =[];

stl1 = stlread('antennas/Dipole10cmT264.stl');
stl2 = stlread('antennas/Dipole10cmT580.stl');
stl3 = stlread('antennas/Dipole10cmT722.stl');
stl4 = stlread('antennas/Dipole10cmT744.stl');
stl5 = stlread('antennas/Dipole10cmT904.stl');
stl6 = stlread('antennas/Dipole10cmT924.stl'); %god
stl7 = stlread('antennas/Dipole10cmT1060.stl'); %god
stl8 = stlread('antennas/Dipole10cmT1104.stl');
stl9 = stlread('antennas/Dipole10cmT1458.stl'); %god 
stl10 = stlread('antennas/Dipole10cmT1680.stl'); 
stl11 = stlread('antennas/Dipole10cmT1922.stl'); %god
stl12 = stlread('antennas/Dipole10cmT2312.stl'); %god
stl13 = stlread('antennas/Dipole10cmT2888.stl'); %god
stl14 = stlread('antennas/Dipole10cmT3528.stl');
% stl1 = stlread('antennas/Dipole10cmT264.stl');
% stl2 = stlread('antennas/Dipole10cmT580.stl');
stl1 = stlread('antennas/Dipole10cmT722.stl');
% stl4 = stlread('antennas/Dipole10cmT744.stl');
% stl5 = stlread('antennas/Dipole10cmT904.stl');
stl2 = stlread('antennas/Dipole10cmT924.stl'); %god
stl3 = stlread('antennas/Dipole10cmT1060.stl'); %god
stl4 = stlread('antennas/Dipole10cmT1104.stl');
% stl9 = stlread('antennas/Dipole10cmT1458.stl'); %god 
% stl10 = stlread('antennas/Dipole10cmT1680.stl'); 
stl5 = stlread('antennas/Dipole10cmT1922.stl'); %god
% stl12 = stlread('antennas/Dipole10cmT2312.stl'); %god
% stl13 = stlread('antennas/Dipole10cmT2888.stl'); %god
% stl14 = stlread('antennas/Dipole10cmT3528.stl');
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
UseDipole = 1;
DipolePoint = [0,0,0];
% If set to one use 81 sub triangles pr element, if 0 use 9
SubSubTri = 1;
% if 1 use fast (but more inacurate) MoM
vectorized = 0;
sub=1;
% Area of radiation
xmin = -2; xmax = 2;
ymin = -2; ymax = 2;
zmin = -2; zmax = 2;
steps = 200;
PointArea = xmax^2/steps;

%% Loop
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
[p, t] = ArbitraryAntenna.RemoveEqualPoints(stl);
%% Calculating dimensions of dipole
minp = min(p);
maxp = max(p);
[maxmaxp, maxaxis] = max(max(p));
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
%%
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
%%
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
[Ei] = ArbitraryAntenna.PointSource(1, 0, 0, w, mu0, k, Center, DipolePoint, [0,1,0], PointArea);
%% MoM
tic;
fprintf('\n')
disp('MoM')
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
[Jface] = ArbitraryAntenna.CurrentCalc(t, EdgeList, w, mu0, a, BasisLA, RhoP, RhoM, RhoP_, RhoM_, sub);
[Jface] = ArbitraryAntenna.CurrentCalc(t, EdgeList, w, mu0, a, BasisLA, RhoP, RhoM);
toc;
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

end
save(FileName, 'J', 'center', 'ExyCrossX', 'ExzCrossZ');
