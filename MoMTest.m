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
length = 0.995;
ant = Antenna(length, 80, 40, 0.0031, [0,0]);
ant1 = Antenna(length, 60, 40, 0.0031, [0,0]);
ant2 = Antenna(length, 40, 40, 0.0031, [0,0]);
% MoM solver object 
MoMobj = MoM();
% Area creation, where the antenna is placed
FirstTestZone = Area(0.0001, 800, 800, -20, 20, -20, 20, mu0);
%% looping through alpha
for alpha=0:2
    alpha
    [ant1, FirstTestZone] = mom2on1(MoMobj, ant1, ant1, FirstTestZone, alpha, k, w, pi/2, 0, 0, mu0);
end
%% Plots
figure(2)
pcolor(FirstTestZone.z, FirstTestZone.x, abs(real(FirstTestZone.Ethethe)))
shading interp
colorbar
caxis([0 1*10^14])
