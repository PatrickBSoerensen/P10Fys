wavelength = 700; %[nm]
eps2 = 12;
eps1 = 1;

a = 200;
steps = 256;
hmax = 2*pi*a/steps;
A=[];
for j=1:steps
    phi1=(j-1)*2*pi/steps;
    phi2=j*2*pi/steps;
    row=[2 a*cos(phi1) a*cos(phi2) a*sin(phi1) a*sin(phi2) 1 0];
    A(j,:)=row;
end

[p e t]=initmesh(A','hmax',hmax)
h=pdemesh(p,e,t)

size(p)
size(e)
size(t)

xv=p(1,:);
yv=p(1,:);
N=length(xv); % Number of points

k0=2*pi/wavelength;
n1=sqrt(eps1);
n2=sqrt(eps2);

%************************
%* Setup incident field *
%************************
E0v=transpose(exp(1i*k0*n1*xv));

%*************************************************
%* Construct matrix with coupling between points *
%*************************************************
gmat=zeros(N,N);
for i=1:N,
    xi=p(1,i); yi=p(2,i);
    for j=1:i,
        xj=p(1,j); yj=p(2,j);
        r_ij=sqrt((xi-xj)^2+(yi-yj)^2);
        gmat(i,j)=1i/4*besselh(0,1,k0*n1*r_ij)*k0^2*(eps2-eps1);
    end
end

%*******************************
%* Setup area of all triangles *
%*******************************
TriangleAreav=zeros(1,length(t));
for k=1:length(t),
   x1=p(1,t(1,k)); y1=p(2,t(1,k));
   x2=p(1,t(2,k)); y2=p(2,t(2,k)); 
   x3=p(1,t(3,k)); y3=p(2,t(3,k));
   dA_dudv=abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
   TriangleAreav(k)=dA_dudv/2;
end

%********************
%* Construct matrix *
%********************
M=zeros(N,N);
for k=1:length(t),
   k 
   % The three corners of the triangle
   i1=t(1,k); i2=t(2,k); i3=t(3,k); 
   x1=p(1,i1); y1=p(2,i1);
   x2=p(1,i2); y2=p(2,i2); 
   x3=p(1,i3); y3=p(2,i3);
   Delta=TriangleAreav(k);
   
   xf=@(u,v) x1+(x2-x1)*u+(x3-x1)*v;
   yf=@(u,v) y1+(y2-y1)*u+(y3-y1)*v;
   
   phi1=@(u,v) 1-u-v;
   phi2=@(u,v) u;
   phi3=@(u,v) v;
   
   Ak=zeros(N,3);
   for i=1:N, 
      if i==i1 || i==i2 || i==i3,
        %point i is on triangle k and rigorous integration is necessary
        xi=p(1,i); yi=p(2,i);
        f1=@(u,v) 1i/4*besselh(0,1,k0*n1*sqrt((xi-xf(u,v)).^2+...
            (yi-yf(u,v)).^2))*k0^2*(eps2-eps1).*phi1(u,v)*Delta*2;
        f2=@(u,v) 1i/4*besselh(0,1,k0*n1*sqrt((xi-xf(u,v)).^2+...
            (yi-yf(u,v)).^2))*k0^2*(eps2-eps1).*phi2(u,v)*Delta*2;
        f3=@(u,v) 1i/4*besselh(0,1,k0*n1*sqrt((xi-xf(u,v)).^2+...
            (yi-yf(u,v)).^2))*k0^2*(eps2-eps1).*phi3(u,v)*Delta*2;
        vmax=@(u) 1-u;
        Ak(i,1)=integral2(f1,0,1,0,vmax);
        Ak(i,2)=integral2(f2,0,1,0,vmax);
        Ak(i,3)=integral2(f3,0,1,0,vmax);
      else
          if i>=i1, f1=gmat(i,i1); else f1=gmat(i1,i); end
          if i>=i2, f2=gmat(i,i2); else f2=gmat(i2,i); end
          if i>=i3, f3=gmat(i,i3); else f3=gmat(i3,i); end          
          Ak(i,1)=Delta/12*(f1*2+f2+f3);
          Ak(i,2)=Delta/12*(f1+f2*2+f3);
          Ak(i,3)=Delta/12*(f1+f2+f3*2);
      end
   end
%   Pk=zeros(3,N);
%   Pk(1,t(1,k))=1;
%   Pk(2,t(2,k))=1;   
%   Pk(3,t(3,k))=1;
%   if k==1,
%      M=Ak*Pk; 
%   else
%      M=M+Ak*Pk; 
%   end
    M(:,i1)=M(:,i1)+Ak(:,1);
    M(:,i2)=M(:,i2)+Ak(:,2);
    M(:,i3)=M(:,i3)+Ak(:,3);
end

Ev=((eye(N)-M)^(-1))*E0v;
pdesurf(p,t,abs(Ev)), colormap(gray(50)), colorbar, axis image

%************************************
%* Plot field outside the scatterer *
%************************************

% Generate mesh for the region outside the scatterer
A=[];
hmax2 = 2*pi*a/steps*2;
for j=1:steps,
    phi1=(j-1)*2*pi/steps;
    phi2=j*2*pi/steps;
    row=[2 a*cos(phi1) a*cos(phi2) a*sin(phi1) a*sin(phi2) 0 1];
    A(j,:)=row;
end

%hmax=40;
count=steps;
Lbox=1200;

Nx=ceil(Lbox/hmax2); Ny=Nx;
for j=1:Nx,
    x1=-Lbox/2+Lbox/Nx*(j-1);
    y1=-Lbox/2;
    x2=-Lbox/2+Lbox/Nx*j;
    y2=-Lbox/2;
    row=[2 x1 x2 y1 y2 1 0];
    A(j+count,:)=row;
end
count=count+Nx;
for j=1:Ny,
    y1=-Lbox/2+Lbox/Ny*(j-1);
    x1=Lbox/2;
    y2=-Lbox/2+Lbox/Ny*j;
    x2=Lbox/2;
    row=[2 x1 x2 y1 y2 1 0];
    A(j+count,:)=row;
end
count=count+Ny;
for j=1:Nx,
    x1=Lbox/2-Lbox/Nx*(j-1);
    y1=Lbox/2;
    x2=Lbox/2-Lbox/Nx*j;
    y2=Lbox/2;
    row=[2 x1 x2 y1 y2 1 0];
    A(j+count,:)=row;
end
count=count+Nx;
for j=1:Ny,
    y1=Lbox/2-Lbox/Ny*(j-1);
    x1=-Lbox/2;
    y2=Lbox/2-Lbox/Ny*j;
    x2=-Lbox/2;
    row=[2 x1 x2 y1 y2 1 0];
    A(j+count,:)=row;
end
steps=steps+Nx;

[p2 e2 t2]=initmesh(A','hmax',hmax2)
%h=pdemesh(p,e,t)

size(p2)
size(e2)
size(t2)

%* Calculate field on outer mesh *
Ev2=[];
length(t2)
for i=1:length(p2),
    i
    xi=p2(1,i); yi=p2(2,i);
    E0=exp(1i*k0*n1*xi);
    E=E0;
    % calculate coupling from point i in the outside mesh to all points in
    % the inside mesh
    xjv=p(1,:); yjv=p(2,:);
    r_ijv=sqrt((xi-xjv).^2+(yi-yjv).^2);
    g_ijv=1i/4*besselh(0,1,k0*n1*r_ijv)*k0^2*(eps2-eps1);
    if min(r_ijv)<1e-3,
       %point coincides with edge of scatterer
       for j=1:length(p),
          if r_ijv(j)<1e-3,
              Ev2(i)=Ev(j);
          end 
       end
    else
        %Integral over scatterer
        for k=1:length(t),
            i1=t(1,k); i2=t(2,k); i3=t(3,k); 
            x1=p(1,i1); y1=p(2,i1); E1=Ev(i1);
            x2=p(1,i2); y2=p(2,i2); E2=Ev(i2);
            x3=p(1,i3); y3=p(2,i3); E3=Ev(i3);
            Delta=TriangleAreav(k);
            f1=g_ijv(i1);
            f2=g_ijv(i2);
            f3=g_ijv(i3);
            E=E+Delta/12*((f1*E1*2+f2*E1+f3*E1)+(f1*E2+f2*E2*2+f3*E2)+...
                (f1*E3+f2*E3+f3*E3*2));
        end
        Ev2(i)=E;
    end
end
        
%hold on
%pdesurf(p2,t2,abs(Ev2')), colormap(gray(50)), colorbar, axis image
%hold off
    
figure
pdeplot(p,e,t,'xydata',abs(Ev)), colormap(gray(50)), colorbar, axis image
hold on
pdeplot(p2,e2,t2,'xydata',abs(Ev2)), colormap(gray(50)), colorbar 
axis image
hold off
set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',28)
xlabel('x [nm]','FontName','Times New Roman','fontsize',28)
ylabel('y [nm]','FontName','Times New Roman','fontsize',28)
set(colorbar,'FontName','Times New Roman','fontsize',28)
colormap(gray(100))

figure
pdeplot(p,e,t,'xydata',abs(Ev)), colormap(gray(50)), colorbar, axis image
set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',28)
xlabel('x [nm]','FontName','Times New Roman','fontsize',28)
ylabel('y [nm]','FontName','Times New Roman','fontsize',28)
set(colorbar,'FontName','Times New Roman','fontsize',28)
colormap(gray(100))
axis([-600 600 -600 600])

figure
pdeplot(p2,e2,t2,'xydata',abs(Ev2)), colormap(gray(50)), colorbar
axis image
set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',28)
xlabel('x [nm]','FontName','Times New Roman','fontsize',28)
ylabel('y [nm]','FontName','Times New Roman','fontsize',28)
set(colorbar,'FontName','Times New Roman','fontsize',28)
colormap(gray(100))
axis([-600 600 -600 600])

figure
x1=-200:200;
y=0;
Eplotx1=tri2grid(p,t,Ev,x1,y);
x2=200:600;
Eplotx2=tri2grid(p2,t2,transpose(Ev2),x2,y);
x3=-600:-200;
Eplotx3=tri2grid(p2,t2,transpose(Ev2),x3,y);
plot(x1,real(Eplotx1),'k','linewidth',2)
hold on
plot(x2,real(Eplotx2),'k','linewidth',2)
plot(x3,real(Eplotx3),'k','linewidth',2)
hold off
xv=[x3 x1 x2];
yv=[Eplotx3 Eplotx1 Eplotx2];
yvr=real(yv);
yvi=imag(yv);
M=[xv' yvr' yvi'];
save Result_cross_sect_tri_a200_eps12_hmax5.m -ascii M


