clc;close all;
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;

length = 0.995;
segmentslin = 10;
segmentscirc = 20;
radii = 0.2;

new = Antenna(length, segmentslin, segmentscirc, radii);
greenres = Green(new, k);
r = sqrt(new.Lin-circshift(new.Lin,-1)+(new.Radii).^2);
green = (exp((1i.*k.*r))./4.*pi.*r).*(1+(1i./r.*k)-1./(r.*k).^2 ...
-(new.Lin-circshift(new.Lin,-1)).^2./r.^2.*(1+3.*1i./(k.*r)-3./(k.*r).^2));
