set(0,'defaulttextinterpreter','latex')
%% ROD NUMBER 1 IS ALWAYS THE DRIVER.

b=[ ];
L=[ ];
delta=[ ];

antalstang=1; %9

%Lengths
L(1)=0.955; %Length of driver
L(2)=1.038; %Length of reflector
L(3)=0.956; %Length of director 1
L(4)=0.932; %etc
L(5)=0.916;
L(6)=0.906;
L(7)=0.897;
L(8)=0.891;
L(9)=0.887;

%Positions
b(1)=(312-312)*10^(-3); %Driver
b(2)=(0-312)*10^(-3); %Reflector
b(3)=(447-312)*10^(-3); %Director 1
b(4)=(699-312)*10^(-3); %Director 2
b(5)=(1050-312)*10^(-3); %etc
b(6)=(1482-312)*10^(-3);
b(7)=(1986-312)*10^(-3);
b(8)=(2553-312)*10^(-3);
b(9)=(3168-312)*10^(-3);

opl=100;

delta=L./N; %Length of segments

Eiz=1;

SE=N/10; 
a=3.1*10^(-3); %Radius of the antenna

eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

N=20; %Number of
f=146.5*10^6; %Frequency
lambda=c/f;
w=2*pi*f;
k=w/c;

f0=@(z) (0.5-z).*(1.0-z)/(0.5-0)/(1.0-0);
f1=@(z) (0-z).*(1.0-z)/(0-0.5)/(1.0-0.5);
f2=@(z) (0-z).*(0.5-z)/(0-1)/(0.5-1);

%% Setting up positions of segments
eta=1;
for n=1:antalstang
 for j=1:N
    for m=0:2
        zv(eta,m+1)=-L(n)/2+delta(n)*(j-1)+delta(n)*m/2;
        zv(eta,4)=b(n);
    end
eta=eta+1;
end
end

DM=eye(N);
DM=circshift(DM,-1);
DMTOT=zeros(N*antalstang,N*antalstang);

for j=1:antalstang
 DMTOT((j-1)*N+1:j*N,(j-1)*N+1:j*N)=DM;
end
%% Setting up incident field

DSE=SE*delta(1); %Area hit by incident E-field

for j=1:N
for m=0:1
z=zv(j,m+1);
E0=0;
if z>=-DSE/2 && z<=0
E0=1i/(w*mu0)*Eiz*(z+(DSE/2))/(DSE/2);
end
if z>0 && z<=DSE/2
E0=1i/(w*mu0)*Eiz*(1-z/(DSE/2));
end
E0v(j,m+1)=E0;
end
end

nulvektor=zeros(N*(antalstang-1),1);

bv=[E0v(:,1);nulvektor(:,1);E0v(:,2);nulvektor(:,1)];
bv(1)=0; bv(2*N)=0;
%% Constructing matrices

B=zeros(2*antalstang*N,2*antalstang*N);

for m=0:2
for mm=0:1
AM=zeros(N,N);
 for n=1:antalstang*N
zns=zv(n,1);
zne=zv(n,3);
for nm=1:antalstang*N
afstand=zv(n,4)-zv(nm,4);
if abs(afstand) < 10^-5, afstand=a; end
zm=zv(nm,mm+1);
R=@(z)sqrt((z-zm).^2+afstand^2);
g=@(z) exp(1i*k*R(z))./(4*pi*R(z));
Gzz=@(z) g(z).*((1+1i./(k*R(z))-1./((k*R(z)).^2)) - ...
((zm-z).^2)./(R(z).^2).*(1+3i./(k*R(z))-3./((k*R(z)).^2)));
if m==0, f=@(z) f0((z-zns)/(zne-zns)); end
if m==1, f=@(z) f1((z-zns)/(zne-zns)); end
if m==2, f=@(z) f2((z-zns)/(zne-zns)); end
A=quadgk(@(z)(Gzz(z).*f(z)),zns,zne);
AM(nm,n)=A;
end
end
if m==0 && mm==0
B(1:antalstang*N,1:antalstang*N)...
=B(1:antalstang*N,1:antalstang*N)+AM;
end
if m==2 && mm==0
B(1:antalstang*N,1:antalstang*N)...
=B(1:antalstang*N,1:antalstang*N)+AM*DMTOT;
end
if m==1 && mm==0
B(1:antalstang*N,antalstang*N+1:2*antalstang*N)...
=B(1:antalstang*N,antalstang*N+1:2*antalstang*N)+AM;
end
if m==0 && mm==1
B(antalstang*N+1:2*antalstang*N,1:antalstang*N)...
=B(antalstang*N+1:2*antalstang*N,1:antalstang*N)+AM;
end
if m==2 && mm==1
B(antalstang*N+1:2*antalstang*N,1:antalstang*N)...
=B(antalstang*N+1:2*antalstang*N,1:antalstang*N)+AM*DMTOT;
end
if m==1 && mm==1
B(antalstang*N+1:2*antalstang*N,antalstang*N+1:2*antalstang*N)...
=B(antalstang*N+1:2*antalstang*N,antalstang*N+1:2*antalstang*N)+AM;
end
end
end
%% Plotting current
dummy1=0;
for j=1:antalstang
B(dummy1+1,:)=zeros(1,2*N*antalstang);
B(1+dummy1,1+dummy1)=1;
dummy1=dummy1+N;
end
av=(B^(-1))*bv;

dummy1=0;
dummy2=N*antalstang;
 for j=1:antalstang
av0(:,j)=av(dummy1+1:dummy1+N);
av1(:,j)=av(dummy2+1:dummy2+N);
dummy1=dummy1+N;
dummy2=dummy2+N;
end

eta=1;
for j1=1:antalstang
zout=[]; Iout=[];
for j=1:N
zout=[zout zv(eta,1) zv(eta,2)];
Iout=[Iout av0(j,j1) av1(j,j1)];
eta=eta+1;
end
zout=[zout zv(eta-1,3)];
Iout=[Iout 0];
figure(j1)
plot(zout,abs(Iout),'k*-')
zouts(:,:,j1)=zout;
Iouts(:,:,j1)=Iout;
if j1==1
title('Driver')
elseif j1==2
title('Reflector') %Change this if no reflector
else
title(sprintf('Director Number %d',j1-2));
end
xlabel('$z$-axis [m]');
ylabel('Absolute Value of the Current');
end
%% Plotting E-field

zp=linspace(-20,20,opl);
pp=linspace(-20,20,opl);
E=zeros(opl,opl);
for j=1:opl %Plotting E using I
zk=zp(j);
for m=1:opl
for eta=1:antalstang
pk=pp(m);
r=@(z) sqrt((z-zk).^2+(pk-b(eta)).^2);
g=@(z) exp(1i*k*r(z))./(4*pi*r(z));
G=@(z) g(z).*((1+1i./(k*r(z))-1./((k*r(z)).^2)) - ...
((zk-z).^2)./(r(z).^2).*(1+3i./(k*r(z))-3./((k*r(z)).^2)));
E(m,j)=E(m,j)+sum(G(zouts(:,:,eta)).*Iouts(:,:,eta)*delta(eta));
end
end
end
figure(antalstang+1)
pcolor(zp,pp,abs(real(E))), shading interp, colorbar, caxis([0 2*10^(-5)]), axis image

title('Far field from a half-wave dipole for $\theta$');
xlabel('$z$-axis [m]');
ylabel('$y$-axis [m]');

%% Far Field Plots
%In order to plot 3D the Antenna Toolbox for matlab is needed
str=12*lambda; %Distance chosen to be far field

%% Plotting xz-plane
theta=linspace(0,2*pi,opl);
E=zeros(1,opl);
for j=1:opl
for j2=1:antalstang
v=theta(j);
 if j2==2, Ry=a; else, Ry=b(j2); end
thetam=@(z) pi/2-atan(z/Ry);
r=@(z) z.*cos(v)+Ry*sin(v);
konst=exp(1i*k*str)/(4*pi*str);
g=@(z) exp(-1i.*k.*r(z));
E(j)=E(j)+sum(-sin(v).*konst.*g(zouts(:,:,j2)).*Iouts(:,:,j2)*delta(j2));
end
end
figure(antalstang+2)
polarplot(theta,abs(E)./max(abs(E)))
title('Far field for a Half-wave dipole for $\theta$');
%% Plotting xy-plane

phi=linspace(0,2*pi,opl);
theta=pi/2;
E=zeros(1,opl);
for j=1:opl
for j2=1:antalstang
v=phi(j);
if j2==2, Ry=a; else Ry=b(j2); end
r=@(z) z.*cos(theta)+Ry*sin(theta)*sin(v);
konst=exp(1i*k*str)/(4*pi*str);
g=@(z) exp(-1i.*k.*r(z));
E(j)=E(j)+sum(-sin(theta).*konst.*g(zouts(:,:,j2)).*Iouts(:,:,j2)*delta(j2));
end
end
figure(antalstang+3)
polarplot(phi,abs(E)./max(abs(E)))
title('Far field for a Half-wave dipole for $\phi$');
%% Far Field 3d

opl1=73;
opl2=37;
phi=linspace(0,2*pi,opl1);
theta=linspace(0,pi,opl2);
E=zeros(opl1,opl2);
for j=1:opl1
for j12=1:opl2
 for j2=1:antalstang
v=phi(j);
v2=theta(j12);
if j2==2, Ry=a; else Ry=b(j2); end
 thetam=@(z) pi/2-atan(z/Ry);
r=@(z) z.*cos(v2)+Ry*sin(v2)*sin(v);
konst=exp(1i*k*str)/(4*pi*str);
g=@(z) exp(-1i.*k.*r(z));
E(j,j12)=E(j,j12)+sum(-sin(v2).*konst.*g(zouts(:,:,j2)).*Iouts(:,:,j2)*delta(j2));
 end
end
end
figure(antalstang+4)
patternCustom(abs(E)./(max(max(abs(E)))),360.*theta./(2*pi),(360.*phi./(2*pi))')
title('Far Field 3D plot');