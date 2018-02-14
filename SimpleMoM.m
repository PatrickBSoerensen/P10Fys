clc;
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
centrum = (0:0);

new = Antenna(length, segmentslin, segmentscirc, radii, centrum);
greenres = Green(new, k);
r = sqrt(new.Lin-circshift(new.Lin,-1)+(new.Radii).^2);
green = (exp((1i.*k.*r))./4.*pi.*r).*(1+(1i./r.*k)-1./(r.*k).^2 ...
-(new.Lin-circshift(new.Lin,-1)).^2./r.^2.*(1+3.*1i./(k.*r)-3./(k.*r).^2));

gam = cos(tHat*zHat)^(-1);

% for i=N:TotalElements
Func1 = @(phimark) cos(alpha*phimark)*exp(-1i*k*R)/R;
Func2 = @(phimark) cos(phimark)*cos(alpha*phimark)*exp(-1i*k*R)/R;
Func3 = @(phimark) sin(phimark)*sin(alpha*mark)*exp(-1i*k*R)/R;
G1 = integral(Func1, 0, pi);
G2 = integral(Func2, 0, pi);
G3 = integral(Func3, 0, pi);

% Ztt = T(i)*T*(sin(gam(i))*sin(gam)*G2+cos(gam(i))*cos(gam)*G1)-1/k^2*TD(i)*TD*G1;
% Zto = T(i)*T*sin(gam(i))*G3+1/k^2*alpha/rho*TD(i)*TD*G1;
% Zot = T(i)*T*sin(gam)*G3+1/k^2*alpha/rho(i)*TD(i)*TD*G1;
% Zoo = T(i)*T*(G2-1/k^2*alpha^2/(rho(i)*rho)*G1);

