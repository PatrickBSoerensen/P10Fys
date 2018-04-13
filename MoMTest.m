%% Setting up constants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
%% Creation of objects
length = 0.1;
%2/5, 2/3, 1*, 2*
lambda=2/3*length;
f=c/lambda;
w=2*pi*f;
k=w/c;
%Antenna(length, pointsline, pointscircle, radii, centre, generator)
ant = Antenna(length, 20, 20, 0.003, [0,0], 1);
ant2D = TwoDAntenna(length, 20, 20, 0.0031, [0,0], 1);
%MoM solver object 
MoMobj = MoM(ant);
%Area creation, where the antenna is placed
%Area(SingularityProtection, xsteps, zsteps, xmin, xmax, zmin, zmax, mu0)
FirstTestZone = Area(0, 1000, 1000, -5, 5, -5, 5, mu0);
%% Removing the circle section for testing
% ant1.CoordTest=ant1.CoordTest(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2,:);
% ant1.T1 = ant1.T1(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2);
% ant1.T2 = ant1.T2(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2);
% ant1.gammaTest = ant1.gammaTest(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2,:);
% ant1.tHatTest = ant1.tHatTest(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2,:);
% ant1.zHatTest = ant1.zHatTest(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2,:);
% ant1.T1D = ant1.T1D(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2);
% ant1.T2D = ant1.T2D(ant1.PointsCircle:ant1.PointsLine+ant1.PointsCircle-2);
%% looping through alpha
for alpha=0:2
    tic;
    alpha
    [MoMobj, ant, FirstTestZone] = momBOR(MoMobj, ant, FirstTestZone, alpha, k, w, pi/2, 0);
    toc;
end
FirstTestZone.Ethethe = FirstTestZone.Ethethe/max(max(FirstTestZone.Ethethe));
%% Plots
figure(2)
pcolor(FirstTestZone.z, FirstTestZone.x, abs(real(FirstTestZone.Ethethe.')))
shading interp
colorbar
% caxis([0 4*10^(4)])
caxis([0 0.1])
rectangle('Position',[-ant.Radii -ant.Length/2 2*ant.Radii ant.Length],'Curvature',1);%Antenna
% figure(3)
% %Middle segment
% plot(abs(ant.Jthe(ant.PointsCircle:ant.PointsLine+ant.PointsCircle-2)), 'k-*')
% figure(4)
% hold on
% %Lower circ
% plot(abs(ant.Jthe(1:ant.PointsCircle-1)), 'b-*')
% %Upper circ
% plot(abs(ant.Jthe(ant.PointsLine+ant.PointsCircle-1:end)), 'r-*')