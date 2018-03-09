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
ant1 = Antenna(length, 200, 100, 0.0031, [0,0], 1);
ant2 = Antenna(length, 10, 10, 0.0031, [0,0], 1);
% MoM solver object 
MoMobj = MoM();
% Area creation, where the antenna is placed
FirstTestZone = Area(1, 1000, 1000, -10, 10, -10, 10, mu0);
SecondTestZone = Area(0, 400, 400, -20, 20, -20, 20, mu0);
%% looping through alpha
for alpha=0:2
    alpha
    [ant1, FirstTestZone] = mombasis(MoMobj, ant1, FirstTestZone, alpha, k, w, pi/2, 0, 0, mu0);
%     [FirstTestZone] = emissionNew(MoMobj, ant1, FirstTestZone, alpha, k, w, 0);
%       [ant2, SecondTestZone] = mombasis(MoMobj, ant2, SecondTestZone, alpha, k, w, pi/2, 0, 0, mu0);
%   [ant2, SecondTestZone] = mom2on1(MoMobj, ant2, ant2, SecondTestZone, alpha, k, w, pi/2, 0, 0, mu0);
end
%% Plots
figure(1)
pcolor(FirstTestZone.z, FirstTestZone.x, abs(real(FirstTestZone.Ethethe)))
shading interp
% caxis([0 2*10^(-3)])
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