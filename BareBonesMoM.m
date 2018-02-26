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
new = Antenna(length, 40, 40, 0.01, [0,0]);

coord = new.Coord;
tHat = new.tHat;
zHat = new.zHat;
%Gamma should also be part of Antenna as it is specific
gamma = new.gamma;
%N should be Antenna property
N = new.Segments;
%Basis triangle functions
%Should be made Antenna property
testingpoints = new.CoordTest;
T1 = (testingpoints(:,1)-coord(1:N-1,1))./coord(1:N-1,3);
T1 = [0;T1];
T2 = (coord(2:N,1)-testingpoints(:,1))./coord(1:N-1,3);
T2 = [T2;0];
T1D = 1./coord(1:N-1,3);
T1D = [T1D;0];
T2D = -1./coord(1:N-1,3);
T2D = [T2D;0];
%Setting size of Z matrix and b vectors
Z = zeros(N,N);
btthe = (1:N);
btthe0 = (1:N);
Ethethe = 0;
%% Setting up field calculations
SingularityProtection = 0.0001;
xsteps = 100;
zsteps = 100;
x = linspace(-20, 20, xsteps);
z = linspace(-20, 20, zsteps);
rz = (z-coord(:,1));
rx = (x-coord(:,2)-SingularityProtection);
%% calculating for alpha 0
alpha = 0
for i=1:N
    for j=1:N
        if i==j
            R = @(phimark) sqrt((coord(i,3)/4)^.2 ... 
            +2*coord(i,2).^2.*(1-cos(phimark)));
        else
            R = @(phimark) sqrt((coord(i,1)-coord(j,1)).^2 ...
            +(coord(i,2)-coord(j,2)).^2+2.*coord(i,2).*coord(j,2).*(1-cos(phimark)));
        end
        
        Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        
        %Should possibly integrate to 2*pi
        G1 = coord(i,3).*coord(j,3).*integral(Func1, 0, pi);
        G2 = coord(i,3).*coord(j,3).*integral(Func2, 0, pi);
        
        Z(i,j) = (T1(i)+T2(i)).*(T1(j)+T2(j)).*(sin(gamma(i)).*sin(gamma(j)).*G2+cos(gamma(i)).*cos(gamma(j)).*G1)-1./k.^2.*(T1D(i)+T2D(i)).*(T1D(j)+T2D(j)).*G1;
    end
    %Indfaldsvinklen af plan bølgen pi/2
    thetai = pi/2;
    
    J0 = besselj(alpha-1, k*coord(i,2)*sin(thetai));
    J1 = besselj(alpha, k*coord(i,2)*sin(thetai));
    J2 = besselj(alpha+1, k*coord(i,2)*sin(thetai));
    %Zero order b
    btthe0(i) = -1i/(w*mu0)*pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
    *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
    *sin(gamma(i))*1i*(J2-J0)-2*sin(thetai)*cos(gamma(i))*J1);
end

invZ = Z^(-1);
xthe = invZ*btthe0.';
ftn = (T1+T2).*exp(1i.*alpha.*0)./coord(:,2);%T, alpha, n. Expansions functio
Jthe=xthe.*ftn;

for i=1:N
        rx = (x-coord(:,2)-SingularityProtection);
        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);

        Ethethe = B/2 * Jthe(i) * btthe0(i)+Ethethe;
        
        rx = (x+coord(:,2)+SingularityProtection);
        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);
        
        Ethethe = B/2 * Jthe(i) * btthe0(i)+Ethethe;
end

%% Actual script
for alpha=1:2
    alpha
    %beta should propably also turn into a loop
    beta = alpha;
    ftan = @(phi) T1.*exp(1i.*alpha.*phi)*tHat;%T, alpha, n. Expansions function
    fpan = @(phi) T2.*exp(1i.*alpha.*phi)*zHat;%Phi, alpha, n. Expansions function
   
    for i=1:N
        for j=1:N
        if i==j
            R = @(phimark) sqrt((coord(i,3)/4)^.2 ... 
            +2*coord(i,2).^2.*(1-cos(phimark)));
        else                        
            R = @(phimark) sqrt((coord(i,1)-coord(j,1)).^2 ...
            +(coord(i,2)-coord(j,2)).^2+2.*coord(i,2).*coord(j,2).*(1-cos(phimark)));
        end
        
        Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
        %Should possibly integrate to 2*pi
        G1 = coord(i,3).*coord(j,3).*integral(Func1, 0, pi);
        G2 = coord(i,3).*coord(j,3).*integral(Func2, 0, pi);
       
        Z(i,j) = (T1(i)+T2(i)).*(T1(j)+T2(j)).*(sin(gamma(i)).*sin(gamma(j)).*G2+cos(gamma(i)).*cos(gamma(j)).*G1)-1./k.^2.*(T1D(i)+T2D(i)).*(T1D(j)+T2D(j)).*G1;
       
        end
        
        J0 = besselj(alpha-1, k*coord(i,2)*sin(thetai));
        J1 = besselj(alpha, k*coord(i,2)*sin(thetai));
        J2 = besselj(alpha+1, k*coord(i,2)*sin(thetai));
    
        btthe(i) = -1i/(w*mu0)*pi*1i^(alpha)*(T1(i)+T2(i))*coord(i,3)...
        *exp(1i*k*coord(i,1)*cos(thetai))*(cos(thetai)...
        *sin(gamma(i))*1i*(J2-J0)-2*sin(thetai)*cos(gamma(i))*J1);
    end
    
    invZ = Z^(-1);
    xthe = invZ*btthe.';
    phi=0;
    ftn = (T1+T2).*exp(1i.*alpha.*phi)./coord(:,2);%T, alpha, n. Expansions function
    
    
    Jthe = Jthe+2*(xthe.*ftn.*cos(alpha.*phi));
    phiS = 0;
    Jthe(1,1) = 0;
    Jthe(end,1) = 0;
    
    for i=1:N
        rx = (x+coord(:,2)+SingularityProtection);
        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);

        Ethethe = B*(xthe(i)*btthe(i))+Ethethe;
        
        rx = (x-coord(:,2)-SingularityProtection);

        r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
        B = -(1i*w*mu0)/(2*pi)*(exp(-1i*k*r)./r);

        Ethethe = B*(xthe(i)*btthe(i))+Ethethe;
    end
end
figure(1)
pcolor(abs(Ethethe))
shading interp
