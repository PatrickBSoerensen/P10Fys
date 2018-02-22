clc;
%% Setting up constants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
% Creation of antenna, length defined for use in area setup
length = 0.995;
new = Antenna(length, 30, 22, 0.01, [0,0]);

coord = new.Coord;
%Unit vectors should be part of Antenna class
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
%Gamma should also be part of Antenna as it is specific
gamma = acos(dot(tHat,zHat,2));
%N should be Antenna property
N = 2*new.SegmentsCircle+new.SegmentsLine-2;
%Basis triangle functions
%Should be made Antenna property
testingpoints = new.CoordTest;
T1 = (testingpoints(:,1)-coord(1:N-1,1))./coord(1:N-1,1);
T1 = [0;T1];
T2 = (coord(2:N,1)-testingpoints(:,1))./coord(1:N-1,1);
T2 = [T2;0];
T1D = 1./coord(1:N-1,3);
T1D = [T1D;0];
T2D = -1./coord(1:N-1,3);
T2D = [T2D;0];
%Setting size of Z matrix and b vectors
Z = zeros(2*N,2*N);
btthe = (1:N);
btphi = (1:N);
bphithe = (1:N);
bphiphi = (1:N);
btthe0 = (1:N);
btphi0 = (1:N);
bphithe0 = (1:N);
bphiphi0 = (1:N);
%% Setting field calculation
xsteps = 100;
zsteps = 100;
x = linspace(-length*2, length*2, xsteps);
z = linspace(-length*2, length*2, zsteps);
rz = (z-coord(:,1));
rx = (x-coord(:,2)-0.001);
%% calculating for alpha 0
alpha = 0
for i=1:N
    for j=1:N
        if i==j
            R = @(phimark) sqrt((coord(i,3)/4)^.2 ... 
            +2*coord(i,2).^2.*(1-cos(phimark)));
        else
            R = @(phimark) sqrt((coord(i,1)-coord(j,1)).^2 ...
            +(coord(i,2)-coord(j,2)).^2-2.*coord(i,2).*coord(j,2).*(1-cos(phimark)));    
        end
    
        Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        Func3 = @(phimark) sin(phimark).*sin(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        %Should possibly integrate to 2*pi
        G1 = coord(i,3).*coord(j,3).*integral(Func1, 0, pi);
        G2 = coord(i,3).*coord(j,3).*integral(Func2, 0, pi);
        G3 = coord(i,3).*coord(j,3).*integral(Func3, 0, pi);
        
        Ztt = (T1(i)+T2(i)).*(T1(j)+T2(j)).*(sin(gamma(i)).*sin(gamma(j)).*G2+cos(gamma(i)).*cos(gamma(j)).*G1)-1./k.^2.*(T1D(i)+T2D(i)).*(T1D(j)+T2D(j)).*G1;
        Ztphi = 1i*(sin(gamma(i))*((T1(i)+T2(i))*(T1(j)+T2(j)))*G3+(1/k^2)*(alpha/coord(j,2))*((T1D(i)+T2D(i))*(T1(j)*T2(j)))*G1);
        Zphit = 1i*(((T1(i)+T2(i)).*(T1(j)+T2(j)))*sin(gamma(j))*G3+(1/k^2)*(alpha/coord(i,2))*(((T1(i)+T2(i))*(T1D(j)+T2D(j))).*G1));
        Zphiphi = -(T1(i)+T2(i)).*(T1(j)+T2(j)).*(G2-1/k.^2.*alpha.^2/(coord(i,2).*coord(j,2)).*G1);
        
        Z(i,j) = Ztt;
        Z(i+N,j) = Ztphi;
        Z(i,j+N) = Zphit;
        Z(i+N,j+N) = Zphiphi;
    end
    %Indfaldsvinklen af plan bølgen pi/2
    thetai = pi/2;
    
    J0 = besselj(alpha-1, k*coord(i,2)*sin(thetai));
    J1 = besselj(alpha, k*coord(i,2)*sin(thetai));
    J2 = besselj(alpha+1, k*coord(i,2)*sin(thetai));
    %Zero order b
    btthe0(i) = pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
    *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
    *sin(gamma(i))*1i*(J2-J0)-2*sin(thetai)*cos(gamma(i))*J1);
    
    bphithe0(i) = -pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
    *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
    *(J2+J0));
    
    btphi0(i) = pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
    *exp(1i*k*coord(i,1)*sin(thetai))*(sin(gamma(i))...
    *(J2-J0));
    
    bphiphi0(i) = pi*1i^(alpha+1)*(T1(i)+T2(i))*coord(i,3)...
    *exp(1i*k*coord(i,1)*sin(thetai))*(J2-J0);
end

invZ = Z^(-1);
bthe = [btthe0, bphithe0];
xthe = invZ*bthe.';
xtthe = xthe(1:N);
xphithe = xthe(N+1:2*N);
bphi = [btphi0,bphiphi0];
xphi = invZ*bphi.';
xtphi = xphi(1:N);
xphiphi = xphi(N+1:2*N);
ftn = (T1+T2).*exp(1i.*alpha.*0);%T, alpha, n. Expansions function
fpn = (T1+T2).*exp(1i.*alpha.*0);%Phi, alpha, n. Expansions function
xNulAlphaThe = xtthe.*ftn;% I tHat retning
xNulAlphaPhi = xtphi.*fpn;% I phiHat retning

Jthe=xtthe.*ftn;
Jphi=xtphi.*fpn;

for i=1:N
        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);

        Ethethe = B/2 * xNulAlphaThe(i) * btthe0(i);
        Ephithe = 0;
        Ethephi = 0;
        Ephiphi = B/2 * xNulAlphaPhi(i) * btphi0(i);
        
        rx = (x+coord(:,2)+0.001);
        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);
        
        Ethethe = B/2 * xNulAlphaThe(i) * btthe0(i)+Ethethe;
        Ephiphi = B/2 * xNulAlphaPhi(i) * btphi0(i)+Ephiphi;
        
end

%% Actual script
for alpha=1:2
    alpha
    %beta should propably also turn into a loop
    beta = alpha;
    ftan = @(z, phi) T1.*exp(1i.*alpha.*phi)*tHat;%T, alpha, n. Expansions function
    fpan = @(z, phi) T2.*exp(1i.*alpha.*phi)*zHat;%Phi, alpha, n. Expansions function
    ftbn = @(z, phi) T1.*exp(-1i.*beta.*phi)*tHat;%T, beta, n. Test function
    Tpbn = @(z, phi) T2.*exp(-1i.*beta.*phi)*zHat;%Phi, beta, n. Test function

    for i=1:N
        for j=1:N
        if i==j
            R = @(phimark) sqrt((coord(i,3)/4)^.2 ... 
            +2*coord(i,2).^2.*(1-cos(phimark)));
        else
            R = @(phimark) sqrt((coord(i,1)-coord(j,1)).^2 ...
            +(coord(i,2)-coord(j,2)).^2-2.*coord(i,2).*coord(j,2).*(1-cos(phimark)));
        end
        
        Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        Func3 = @(phimark) sin(phimark).*sin(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        %Should possibly integrate to 2*pi
        G1 = coord(i,3).*coord(j,3).*integral(Func1, 0, pi);
        G2 = coord(i,3).*coord(j,3).*integral(Func2, 0, pi);
        G3 = coord(i,3).*coord(j,3).*integral(Func3, 0, pi);
       
        Ztt = (T1(i)+T2(i)).*(T1(j)+T2(j)).*(sin(gamma(i)).*sin(gamma(j)).*G2+cos(gamma(i)).*cos(gamma(j)).*G1)-1./k.^2.*(T1D(i)+T2D(i)).*(T1D(j)+T2D(j)).*G1;
        Ztphi = 1i*(sin(gamma(i))*((T1(i)+T2(i))*(T1(j)+T2(j)))*G3+(1/k^2)*(alpha/coord(j,2))*((T1D(i)+T2D(i))*(T1(j)*T2(j)))*G1);
        Zphit = 1i*(((T1(i)+T2(i)).*(T1(j)+T2(j)))*sin(gamma(j))*G3+(1/k^2)*(alpha/coord(i,2))*(((T1(i)+T2(i))*(T1D(j)+T2D(j))).*G1));
        Zphiphi = -(T1(i)+T2(i)).*(T1(j)+T2(j)).*(G2-1/k.^2.*alpha.^2/(coord(i,2).*coord(j,2)).*G1);
   
        Z(i,j) = Ztt;
        Z(i+N,j) = Ztphi;
        Z(i,j+N) = Zphit;
        Z(i+N,j+N) = Zphiphi;
        end
        
        J0 = besselj(alpha-1, k*coord(i,2)*sin(thetai));
        J1 = besselj(alpha, k*coord(i,2)*sin(thetai));
        J2 = besselj(alpha+1, k*coord(i,2)*sin(thetai));
    
        btthe(i) = pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
        *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
        *sin(gamma(i))*1i*(J2-J0)-2*sin(thetai)*cos(gamma(i))*J1);
    
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
    xthe = invZ*bthe.';
    xtthe = xthe(1:N);
    xphithe = xthe(N+1:2*N);
    bphi = [btphi,bphiphi];
    xphi = invZ*bphi.';
    xtphi = xphi(1:N);
    xphiphi = xphi(N+1:2*N);
    ftn = (T1+T2).*exp(1i.*alpha.*0);%T, alpha, n. Expansions function
    fpn = (T1+T2).*exp(1i.*alpha.*0);%Phi, alpha, n. Expansions function
    phi=0;

    Jthe = Jthe+2*(xtthe.*ftn.*cos(alpha.*phi)+xphithe.*ftn.*sin(alpha.*phi));
    Jphi = Jphi+2*(xtphi.*fpn.*sin(alpha.*phi)+xphiphi.*fpn.*cos(alpha.*phi));
    phiS = 0;
    
    for i=1:N
        rx = (x+coord(:,2)+0.001);
        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);

        Ethethe = B*(xtthe(i)*btthe(i)+xphithe(i)*bphithe(i))*cos(alpha*phiS)+Ethethe;
        Ephithe = 1i*B*(xtthe(i)*btphi(i)+xphithe(i)*bphiphi(i))*sin(alpha*phiS)+Ephithe;
        Ethephi = 1i*B*(xtphi(i)*btthe(i)+xphiphi(i)*bphithe(i))*sin(alpha*phiS)+Ethephi;
        Ephiphi = B*(xtphi(i)*btphi(i)+xphiphi(i)*bphiphi(i))*cos(alpha*phiS)+Ephiphi;
        
        rx = (x-coord(:,2)-0.001);

        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);

        Ethethe = B*(xtthe(i)*btthe(i)+xphithe(i)*bphithe(i))*cos(alpha*phiS)+Ethethe;    
        Ephithe = 1i*B*(xtthe(i)*btphi(i)+xphithe(i)*bphiphi(i))*sin(alpha*phiS)+Ephithe;
        Ethephi = 1i*B*(xtphi(i)*btthe(i)+xphiphi(i)*bphithe(i))*sin(alpha*phiS)+Ethephi;
        Ephiphi = B*(xtphi(i)*btphi(i)+xphiphi(i)*bphiphi(i))*cos(alpha*phiS)+Ephiphi;
    end
end
figure(1)
pcolor(abs(Ethethe))
shading interp
% figure(2)
% pcolor(real(Ethephi))
% shading interp
% figure(3)
% pcolor(real(Ephithe))
% shading interp    
figure(4)
pcolor(abs(Ephiphi))
shading interp
