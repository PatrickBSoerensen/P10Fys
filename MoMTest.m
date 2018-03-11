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
ant1 = Antenna(length, 50, 20, 0.0031, [0,0], 1);
% MoM solver object 
MoMobj = MoM(ant1);
% Area creation, where the antenna is placed
FirstTestZone = Area(1, 1000, 1000, -5, 5, -5, 5, mu0);
%% looping through alpha
for alpha=0:2
    alpha
    [MoMobj, ant1, FirstTestZone] = mombasis(MoMobj, ant1, FirstTestZone, alpha, k, w, pi/2, 0, 0, mu0);
end
%% Plots
figure(1)
pcolor(FirstTestZone.z, FirstTestZone.x, abs(real(FirstTestZone.Ethethe)))
shading interp
colorbar
caxis([0 6*10^(-4)])
rectangle('Position',[-ant1.Radii -ant1.Length/2 2*ant1.Radii ant1.Length],'Curvature',1);%Antenna
figure(2)
%Middle segment
plot(abs(ant1.Jthe(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2)), 'k-*')
figure(3)
hold on
%Lower circ
plot(abs(ant1.Jthe(1:ant1.PointsCircle-1)), 'b-*')
%Upper circ
plot(abs(ant1.Jthe(ant1.PointsLine+ant1.PointsCircle-2:end)), 'r-*')
figure(4)
plot(ant1.E0)