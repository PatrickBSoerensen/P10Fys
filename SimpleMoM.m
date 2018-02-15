clc;
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
r=1;
new = Antenna(0.995, 10, 20, 0.2, [0,0]);

green = (exp((1i.*k.*r))./4.*pi.*r).*(1+(1i./r.*k)-1./(r.*k).^2 ...
-(new.Lin-circshift(new.Lin,-1)).^2./r.^2.*(1+3.*1i./(k.*r)-3./(k.*r).^2));

coord = CreateCoord(new);
testing = linspace(-pi/2, 0, new.SegmentsCircle);
tHat(1:new.SegmentsCircle,1) ... 
    = -new.Radii.*sin(linspace(-pi/2, 0, new.SegmentsCircle));%x coord
tHat(new.SegmentsCircle+1:new.SegmentsCircle+new.SegmentsLine-2,1) ... 
    = 0;%x coord
tHat(new.SegmentsCircle+new.SegmentsLine-1: ... 
    2*new.SegmentsCircle+new.SegmentsLine-2,1)...
    = -new.Radii.*sin(linspace(0, pi/2, new.SegmentsCircle));%x coord
tHat(:,2) = 0;%y coord
tHat(:,3) = 1;%z coord

tHat(1:new.SegmentsCircle,3) ... 
    = new.Radii.*cos(linspace(-pi/2, 0, new.SegmentsCircle));%z coord
tHat(new.SegmentsCircle+1:new.SegmentsCircle+new.SegmentsLine-2,3) ... 
    = 1;%z coord
tHat(new.SegmentsCircle+new.SegmentsLine-1: ... 
    2*new.SegmentsCircle+new.SegmentsLine-2,3)...
    = new.Radii.*cos(linspace(0, pi/2, new.SegmentsCircle));%z coord
tHat = tHat./sqrt(tHat(:,1).^2+tHat(:,2).^2+tHat(:,3).^2);
zHat = tHat;

zHat(:,1) = 0; %x coord
zHat(:,2) = 0; %y coord
zHat(:,3) = 1; %z coord

gamma = acos(dot(tHat,zHat,2));

alpha = 1;
beta = 2;

f = coord(:,1:2);
%Should be changed to evaluate in centre points
f1 = @(z)(z-circshift(coord(:,1),1))/(coord(:,1)-circshift(coord(:,1),1));
f2 = @(z)(circshift(coord(:,1),-1)-z)/(circshift(coord(:,1),-1)-coord(:,1));

T1 = @(z, phi) f1.*exp(1i.*alpha.*phi)*tHat;
T2 = @(z, phi) f2.*exp(1i.*alpha.*phi)*zHat;
T3 = @(z, phi) f1.*exp(-1i.*beta.*phi)*tHat;
T4 = @(z, phi) f2.*exp(-1i.*beta.*phi)*zHat;

TD = 4;
wo = linspace(-pi/2, 0, new.SegmentsCircle);
TDtbm = 1/(new.Radii*cos(wo)).*-new.Radii*sin(wo).*f1.*exp(1i*beta*phi);
TDtan = 1/(new.Radii*cos(wo)).*-new.Radii*sin(wo).*f1.*exp(1i*alpha*phimark);
TDpbm = -f1./(new.Radii*cos(wo)).*-new.Radii*sin(wo).*1i.*beta.*exp(1i*beta*phi);
TDpan = 1i.*alpha.*f1/(new.Radii*cos(wo)).*exp(1i*alpha*phimark);

Z = [];

% 2*new.SegmentsCircle+new.SegmentsLine
for i=1:new.SegmentsCircle+new.SegmentsLine
    R = sqrt((coord(i,1)-coord(:,1)).^2+(coord(i,2)-coord(:,2)).^2-2.*coord(i,2).*coord(i,2).*cos(0-0));
    
    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R)./R;
    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R)./R;
    Func3 = @(phimark) sin(phimark).*sin(alpha.*mark).*exp(-1i.*k.*R)./R;
    
    G1 = @(y1)arrayfun(@(phimark)coord(i,3).*coord(:,3).*integral(Func1, 0, pi),y1);
    G2 = @(y1)arrayfun(@(phimark)coord(i,3).*coord(:,3).*integral(Func2, 0, pi),y1);
    G3 = @(y1)arrayfun(@(phimark)coord(i,3).*coord(:,3).*integral(Func3, 0, pi),y1);
    
    Ztt = @(z)(T1(i)+T2(i))*T1*(sin(gamma(i))*sin(gamma)*G2+cos(gamma(i))*cos(gamma)*G1)-1/k^2*TD(i)*TD*G1;
    Zto = @(z)(T1(i)+T2(i))*T1*sin(gamma(i))*G3+1/k^2*alpha/rho*TD(i)*TD*G1;
    Zot = @(z)(T1(i)+T2(i))*T1*sin(gamma)*G3+1/k^2*alpha/rho(i)*TD(i)*TD*G1;
    Zoo = @(z)(T1(i)+T2(i))*T1*(G2-1/k^2*alpha^2/(rho(i)*rho)*G1);

    
end
