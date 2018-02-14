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

% greenres = Green(new, k);
r = sqrt(new.Lin-circshift(new.Lin,-1)+(new.Radii).^2);
green = (exp((1i.*k.*r))./4.*pi.*r).*(1+(1i./r.*k)-1./(r.*k).^2 ...
-(new.Lin-circshift(new.Lin,-1)).^2./r.^2.*(1+3.*1i./(k.*r)-3./(k.*r).^2));
<<<<<<< HEAD

coord = CreateCoord(new);
tHat(:,1) = -new.Radii*sin(coord(:,3));
tHat(:,2) = 0;
tHat(:,3) = coord(:,1);
zHat = pi/2;
gam = cos(tHat*zHat)^(-1);
=======
alpha = 1;
f = coord(:,1:2);
f1 = @(z)(z-circshift(coord(:,1),1))/(coord(:,1)-circshift(coord(:,1),1));
f2 = @(z)(circshift(coord(:,1),-1)-z)/(circshift(coord(:,1),-1)-coord(:,1));
T1 = @(z, phi) f1.*exp(1i.*alpha.*phi)*tHat;
T2 = @(z, phi) f2.*exp(1i.*alpha.*phi)*zHat;
TD = 4;


% 2*new.SegmentsCircle+new.SegmentsLine
for i=1:new.SegmentsCircle+new.SegmentsLine
    R = sqrt((coord(i,1)-coord(:,1)).^2+(coord(i,2)-coord(:,2)).^2-2.*coord(i,2).*coord(i,2).*cos(0-0));
    
    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R)./R;
    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R)./R;
    Func3 = @(phimark) sin(phimark).*sin(alpha.*mark).*exp(-1i.*k.*R)./R;
    
    G1 = @(y1)arrayfun(@(phimark)coord(i,3).*coord(:,3).*integral(Func1, 0, pi),y1);
    G2 = @(y1)arrayfun(@(phimark)coord(i,3).*coord(:,3).*integral(Func2, 0, pi),y1);
    G3 = @(y1)arrayfun(@(phimark)coord(i,3).*coord(:,3).*integral(Func3, 0, pi),y1);
    
    Zoo = @(z)(T1(i)+T2(i))*T1*(G2-1/k^2*alpha^2/(rho(i)*rho)*G1);

    
end
