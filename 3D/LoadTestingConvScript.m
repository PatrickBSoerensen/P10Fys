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
stl6 = stlread('antennas/Dipole10cmT1104.stl');

for convloop=1:6
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
radius = abs(min(min(p))+min(min(p)));
Length = (maxmaxp-minp(maxaxis));%-2*radius;
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
SubSubTri = 0;
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
%% MoM
vectorized = 1;
tic;
fprintf('\n')
disp('MoM')
if vectorized
    [Z, a, b ] = ArbitraryAntenna.MoMVectorized(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, 0, 1, 0);
else
    [Z, b, a] = ArbitraryAntenna.MoM(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k,  SubTri, 0, 1, 0);
end
toc;
%% Current calc in Triangle
sub = 1;
tic;
fprintf('\n')
disp('Calculating Current')
[Jface] = ArbitraryAntenna.CurrentCalc(t, EdgeList, w, mu0, a, BasisLA, RhoP, RhoM, RhoP_, RhoM_, sub);
toc;

%% Calculating E
%x-min/max, z-min/max, y-min/max
tic;
fprintf('\n')
disp('Calculating E-field')
[Exy, Exz, Ezy, x, y, z, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = ...
    ArbitraryAntenna.EField(Center, k, Jface, -2, 2, -2, 2, -2, 2, 100, Area);
toc;

sx = size(J);
sy = size(Jface);
a = max(sx(1),sy(1));
J = [[J;zeros(abs([a 0]-sx))],[Jface;zeros(abs([a,0]-sy))]];

sx = size(center);
sy = size(Center);
a = max(sx(1),sy(1));
center = [[center ;zeros(abs([a 0]-sx))],[Center;zeros(abs([a,0]-sy))]];
ExyCrossX= [ExyCrossX Exy(:,50)];
ExzCrossZ = [ExzCrossZ  Exz(:,50)];

end
save('ConvFast', 'J', 'center', 'ExyCrossX', 'ExzCrossZ');
