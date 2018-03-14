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
%5/2, 3/2, 1*, 2*
lambda=length/4;
f=c/lambda;
w=2*pi*f;
k=w/c;
%Antenna(length, pointsline, pointscircle, radii, centre, generator)
ant1 = Antenna(length, 20, 5, 0.0031, [0,0], 1);
%MoM solver object 
MoMobj = MoM(ant1);
%Area creation, where the antenna is placed
%Area(SingularityProtection, xsteps, zsteps, xmin, xmax, zmin, zmax, mu0)
FirstTestZone = Area(0.00001, 500, 500, -5, 5, -5, 5, mu0);
%% looping through alpha
for alpha=0:2
    alpha
    [MoMobj, ant1, FirstTestZone] = mombasis(MoMobj, ant1, FirstTestZone, alpha, k, w, pi/2, 0, 0, mu0);
end
% Plots
figure(1)
pcolor(FirstTestZone.x, FirstTestZone.z, abs(real(FirstTestZone.Ethethe)))
shading interp
colorbar
% caxis([0 2*10^(-4)])
rectangle('Position',[-ant1.Radii -ant1.Length/2 2*ant1.Radii ant1.Length],'Curvature',1);%Antenna
figure(2)
%Middle segment
plot(abs(ant1.Jthe(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2)), 'k-*')
figure(3)
hold on
%Lower circ
plot(abs(ant1.Jthe(1:ant1.PointsCircle-1)), 'b-*')
%Upper circ
plot(abs(ant1.Jthe(ant1.PointsLine+ant1.PointsCircle-1:end)), 'r-*')
figure(4)
plot(ant1.E0)