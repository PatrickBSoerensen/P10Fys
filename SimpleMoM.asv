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

phiHat(:,1) = 0; %x coord
phiHat(:,2) = 1; %y coord
phiHat(:,3) = 0; %z coord

gamma = acos(dot(tHat,zHat,2));

alpha = 1;
beta = 2;

f = coord(:,1:2);
%Basis triangle functions
%Should be changed to evaluate in centre points
testingpoints = CreateTestCoord(new);
testend = size(testingpoints(:,1));
T1 = (testingpoints(:,1)-coord(1:testend,1))./coord(1:testend,1);
T1 = [0;T1];
T2 = (coord(2:size(coord),1)-testingpoints(:,1))./coord(1:testend,1);
T2 = [T2;0];
T1D = 1/coord(1:testend,1);
T1D = [T1D,0];
T2D = -1/coord(1:testend,1);
T2D = [T2D,0];

FT = [T1,T2];%Full T functions

ftan = @(z, phi) T1.*exp(1i.*alpha.*phi)*tHat;%T, alpha, n. Expansions function
fpan = @(z, phi) T2.*exp(1i.*alpha.*phi)*zHat;%Phi, alpha, n. Expansions function
ftbn = @(z, phi) T1.*exp(-1i.*beta.*phi)*tHat;%T, beta, n. Test function
Tpbn = @(z, phi) T2.*exp(-1i.*beta.*phi)*zHat;%Phi, beta, n. Test function

N = 2*new.SegmentsCircle+new.SegmentsLine-2;
Z = zeros(2*N,2*N);
btthe = (1:N);
btphi = (1:N);
bphithe = (1:N);
bphiphi = (1:N);

for i=1:N
    for j=1:N
    if i==j
       R = @(phimark) sqrt((coord(i,3)/4)^.2 ... 
       +2*coord(i,2).^2.*(1-cos(phimark)));
    else
       R = @(phimark) sqrt((coord(i,1)-coord(j,1)).^2 ...
       +(coord(i,2)-coord(j,2)).^2-2.*coord(i,2).*coord(j,2).*cos(1-phimark));    
    end
    
    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
    Func3 = @(phimark) sin(phimark).*sin(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
    
    G1 = coord(i,3).*coord(j,3).*integral(Func1, 0, pi);
    G2 = coord(i,3).*coord(j,3).*integral(Func2, 0, pi);
    G3 = coord(i,3).*coord(j,3).*integral(Func3, 0, pi);
       
    Ztt = (T1(i)+T2(i)).*(T1(j)+T2(j)).*(sin(gamma(i)).*sin(gamma(j)).*G2+cos(gamma(i)).*cos(gamma(j)).*G1)-1./k.^2.*(T1D(i)+T2D(i)).*(T1D(j)+T2D(j)).*G1;
    Zto = 1i.*((T1(i)+T2(i)).*(T1(j)+T2(j)).*sin(gamma(i)).*G3+1./k.^2.*alpha./coord(j,2).*(T1D(i)+T2D(i)).*(T1(j)+T2(j)).*G1);
    Zot = 1i.*((T1(i)+T2(i)).*(T1(j)+T2(j)).*sin(gamma(j)).*G3+1./k.^2.*alpha./coord(i,2).*(T1(i)+T2(i)).*(T1D(j)+T2D(j)).*G1);
    Zoo = -(T1(i)+T2(i)).*(T1(j)+T2(j)).*(G2-1/k.^2.*alpha.^2/(coord(i,2).*coord(j,2)).*G1);
   
    Z(i,j) = Ztt;
    Z(i+N,j) = Zto;
    Z(i,j+N) = Zot;
    Z(i+N,j+N) = Zoo;
    end
    thetai = pi/2;
    J0 = besselj(0, k*coord(i,2)*sin(thetai));
    J1 = besselj(1, k*coord(i,2)*sin(thetai));
    J2 = besselj(2, k*coord(i,2)*sin(thetai));
    
    btthe(i) = pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
        *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
        *sin(gamma(i)*1i*(J2-J0)-2*sin(thetai)*cos(gamma(i))*J1));
    
    bphithe(i) = -pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
        *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
        *(J2+J0));
    
    btphi(i) = pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
        *exp(1i*k*coord(i,1)*sin(thetai))*(sin(gamma(i))...
        *(J2-J0));
    
    bphiphi(i) = pi*1i^(alpha+1)*(T1(i)+T2(i))*coord(i,3)...
        *exp(1i*k*coord(i,1)*sin(thetai))*(J2-J0);
end
invZ = Z^(-1);
bthe = [btthe, bphithe];
xtthe = invZ*bthe.';
xphithe = invZ*bphithe.';
bphi = [btphi,bphiphi];
xtphi = invZ*bphi.';
xphiphi = invZ*bphiphi.';
ftn = T1.*exp(1i.*alpha.*0);%T, alpha, n. Expansions function
fpn = T2.*exp(1i.*alpha.*0);%Phi, alpha, n. Expansions function
J1 = sum(xtthe(1:48).*ftn+xtthe(49:96).*fpn);
Jthe = xtthe*ftn+2*xtthe;
Jphi

% b calculations for non plane waves 
% bt = -1i./(w.*mu0).*(T1.*T2).*coord(i,3).*integral((@(phi) exp(-1i.*alpha.*phi)),0,2*pi).*tHat.*E(r);
% bo = -1i./(w.*mu0).*(T1.*T2).*coord(i,3).*integral((@(phi) exp(-1i.*alpha.*phi)),0,2*pi).*phiHat.*E(r);
