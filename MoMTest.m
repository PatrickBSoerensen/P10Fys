clc;
%% Setting up constants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
% Creation of antenna, length defined for use in area setup
length = 0.995;
ant = Antenna(length, 40, 40, 0.0031, [0,0]);
MoMobj = MoM(0,0);
FirstTestZone = area(0.0001, 800, 800, -20, 20, -20, 20, mu0);

for alpha=0:2
    alpha
    ant = momself(MoMobj, ant, alpha, k, pi/2, 0);
    FirstTestZone = emission(MoMobj, ant, FirstTestZone, alpha, k, w, 0);
end

figure(1)
pcolor(FirstTestZone.z, FirstTestZone.x, abs(real(FirstTestZone.Ethethe)))
shading interp
colorbar
% caxis([0 0.2])