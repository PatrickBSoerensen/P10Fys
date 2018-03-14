clc;
%% Setting up constants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
%% Creation of objects
length = lambda/2;
LinSeg = linspace(2,52,51);
CircSeg = linspace(2,52,51);
Emax = zeros(51, 51);
parfor lin=1:51
for circ=1:51
%Antenna(length, pointsline, pointscircle, radii, centre, generator)
ant1 = Antenna(length, LinSeg(lin), CircSeg(circ), 0.0031, [0,0], 1);
% MoM solver object 
MoMobj = MoM(ant1);
% Area creation, where the antenna is placed
%Area(SingularityProtection, xsteps, zsteps, xmin, xmax, zmin, zmax, mu0)
FirstTestZone = Area(0.000001, 500, 500, -5, 5, -5, 5, mu0);
%% looping through alpha
for alpha=0:2
    alpha
    [MoMobj, ant1, FirstTestZone] = mombasis(MoMobj, ant1, FirstTestZone, alpha, k, w, pi/2, 0, 0, mu0);
end
Emax(lin,circ)=max(max(abs(FirstTestZone.Ethethe)));
end
end