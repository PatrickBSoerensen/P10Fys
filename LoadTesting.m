%% konstants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
%% load STL file into matlab
stl = stlread('AntBinMesh.stl');
%% faces and unique vertices
tic;
[p, t] = ArbitraryAntenna.RemoveEqualPoints(stl);
t = sort(t,2);
toc;
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
[Z, b, x] = ArbitraryAntenna.MoM(p, t, EdgeList, BasisNumber, BasisLA, A, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri);
toc;
            
    %% Emission
%     r = linspace(0,10,100);
%     E = zeros(100,100);
%              for i=1:length(r)
%                 
%                 g = exp(-1i.*k.*r)./(4.*pi.*r);
%                 G = g.*((1+1i./(k*r)-1./((k*r).^2)) - ...
%                             ((rx(i,:)).^2)./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
%                 E(i,:) = E(i,:)+G.*ant.Jthe(i);
%              end

