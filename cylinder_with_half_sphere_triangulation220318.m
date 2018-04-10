clear all; close all

% Radius of cylinder and half-spheres
a = 5; %[nm]
steps = 2; % Discretization along theta and phi for a quarter half-sphere

% Length of cylinder
H = 10;
stepsH = 2; %Discretization along the cylinder height

add_triangle_normal_vector = false;
add_element_type_and_details = true;

add_points_for_second_order_interpolation = true;


%figure(5)
%**********************************
%* Mesh for a quarter half-sphere *
%**********************************
hmax = pi/2*a/steps*1.3;
A=[];
for j=1:steps
    t1a=pi/2*a*(j-1)/steps;
    t2a=-a*pi/4*sin(t1a/a);
    t1b=pi/2*a*(j)/steps;
    t2b=-a*pi/4*sin(t1b/a);
    row=[2 t2a t2b t1a t1b 1 0];
    A(j,:)=row;
%    plot(t2a,t1a,'k*'); hold on
end
for j=1:steps
    t1a=pi/2*a; t1b=t1a;
    t2a=-a*pi/4+a*pi/2*(j-1)/steps;
    t2b=-a*pi/4+a*pi/2*j/steps;
    row=[2 t2a t2b t1a t1b 1 0];
    A(j+steps,:)=row;
%    plot(t2a,t1a,'r*');    
end
for j=1:steps
    t1a=pi/2*a*(steps-j+1)/steps;
    t2a=a*pi/4*sin(t1a/a);
    t1b=pi/2*a*(steps-j)/steps;
    t2b=a*pi/4*sin(t1b/a);
    row=[2 t2a t2b t1a t1b 1 0];
    A(j+2*steps,:)=row;
%    plot(t2a,t1a,'b*');   
end

figure(1)
[p, e, t]=initmesh(A','hmax',hmax);
h=pdemesh(p,e,t);
size(p)
size(e)
size(t)
set(h,'linewidth',2,'color','k')
axis image
xlabel('t_2','fontsize',28,'FontName','Times New Roman')
ylabel('t_1','fontsize',28,'FontName','Times New Roman')
set(gca,'fontsize',28,'FontName','Times New Roman')
set(gca,'linewidth',1)

%hold on
%plot(p(1,:),p(2,:),'k*')
%hold off

%*******************************
%* Mesh for a quarter cylinder *
%*******************************
if H>0
A=[];
for j=1:steps
    t1a=0;
    t2a=pi/2*a*(j-1)/steps;
    t1b=0;
    t2b=pi/2*a*(j)/steps;
    row=[2 t2a t2b t1a t1b 1 0];
    A(j,:)=row;
end
for j=1:stepsH
    t1a=H*(j-1)/stepsH; t1b=H*j/stepsH;
    t2a=pi/2*a;
    t2b=pi/2*a;
    row=[2 t2a t2b t1a t1b 1 0];
    A(j+steps,:)=row;
end
for j=1:steps
    t1a=H; 
    t2a=pi/2*a*(steps-j+1)/steps;
    t1b=H;
    t2b=pi/2*a*(steps-j)/steps;
    row=[2 t2a t2b t1a t1b 1 0];
    A(j+steps+stepsH,:)=row;
end
for j=1:stepsH
    t1a=H*(stepsH-j+1)/stepsH; t1b=H*(stepsH-j)/stepsH;
    t2a=0;
    t2b=0;
    row=[2 t2a t2b t1a t1b 1 0];
    A(j+2*steps+stepsH,:)=row;
end
%hmax=H/stepsH;
figure(2)
[pc, ec, tc]=initmesh(A','hmax',hmax);
h=pdemesh(pc,ec,tc);
set(h,'linewidth',2,'color','k')
axis image


% Convert to 3D
p3=zeros(3,length(p));
for j=1:length(p)
    t2=p(1,j); t1=p(2,j);
    z=a*cos(t1/a);
    x=a*sin(t1/a)*cos(t2/(a*sin(t1/a))-pi/4-pi/2);
    y=a*sin(t1/a)*sin(t2/(a*sin(t1/a))-pi/4-pi/2);
    if t1==0
        x=0; y=0;
    end
    p3(1,j)=x;
    p3(2,j)=y;
    p3(3,j)=z;
end
end

% Construct half-sphere surfaces
p3=[];
t3=[];
counter = 0;
for jflip=1:2
    for jphi=1:4
        phi=pi/2*(jphi-1);
        for j=1:length(p)
            xc=0; yc=0; zc=0;
            t2=p(1,j); t1=p(2,j);
            if jflip==1
                z=a*cos(t1/a)+H;
                zc=H;
            else
                z=-a*cos(t1/a);
            end
            x=a*sin(t1/a)*cos(t2/(a*sin(t1/a))-pi/4-pi/2+phi);
            y=a*sin(t1/a)*sin(t2/(a*sin(t1/a))-pi/4-pi/2+phi);
            if t1==0
                x=0; y=0;
            end
            p3(1,j+counter)=x;
            p3(2,j+counter)=y;
            p3(3,j+counter)=z;
        end
        % 1 : means part of spherical surface, xc, yc, zc: center of sphere
        % The normal vector is assumed to be in the direction from the
        % sphere center to the point on the surface
        td=zeros(7,length(t));
        td(1,:)=1; td(2,:)=xc; td(3,:)=yc; td(4,:)=zc; 
        td(5,:)=0; td(6,:)=0; td(7,:)=0;
        t2=[t(1,:)+counter;t(2,:)+counter;t(3,:)+counter;td]; 
        t3=[t3 t2];
        counter=counter+length(p);
    end
end

if H>0
% Add cylinder parts
for jphi=1:4
    xc=0; yc=0; zc=0;
    nx=0; ny=0; nz=1;
    phi=pi/2*(jphi-1);
    for j=1:length(pc)
        t1=pc(1,j); t2=pc(2,j);
        z=t2;
        x=a*cos(t1/a+phi);
        y=a*sin(t1/a+phi);
        p3(1,j+counter)=x;
        p3(2,j+counter)=y;
        p3(3,j+counter)=z;
    end
    % 2 : means part of cylinder surface, xc, yc, zc: point on cylinder
    % axis. nx, ny, nz: direction of cylinder axis.
    
    td=zeros(7,length(tc));
    td(1,:)=2; td(2,:)=xc; td(3,:)=yc; td(4,:)=zc; 
    td(5,:)=nx; td(6,:)=ny; td(7,:)=nz;
    
    t2=[tc(1,:)+counter;tc(2,:)+counter;tc(3,:)+counter; td];
    t3=[t3 t2];
    counter=counter+length(pc);
end
end

t=t3;

pfin=p3; tfin=t3;
% Remove points that appear twice
counter=1;
j2=1;
while j2<length(pfin)
  x=pfin(1,j2); y=pfin(2,j2); z=pfin(3,j2); 
  found=false;
  jfound=0;
   for j1=1:(j2-1)
      if abs(pfin(1,j1)-x)<1e-7
          if abs(pfin(2,j1)-y)<1e-7
              if abs(pfin(3,j1)-z)<1e-7
                 found=true; 
                 jfound=j1;
              end
          end
      end
   end
   if found==false
       % do nothing
       j2=j2+1;
   end
   if found==true
       % Remove point
       pd1=pfin(:,1:j2-1); pd2=pfin(:,j2+1:length(pfin));
       pfin=[pd1 pd2];
       for jt1=1:3
           for jt2=1:length(tfin)
               if tfin(jt1,jt2)==j2, tfin(jt1,jt2)=jfound; end
               if tfin(jt1,jt2)>j2, tfin(jt1,jt2)=tfin(jt1,jt2)-1; end
           end
       end
   end
end

p3=pfin; t=tfin;



if add_element_type_and_details == true
    % Do nothing
end
if add_triangle_normal_vector == true
   tfin2=tfin(1:6,1:length(tfin));
   for k=1:length(tfin)
      if tfin(4,k)==1 || tfin(4,k)==2
          %part of sphere
          xc=tfin(5,k); yc=tfin(6,k); zc=tfin(7,k);
          x1=pfin(1,tfin(1,k)); y1=pfin(2,tfin(1,k)); z1=pfin(3,tfin(1,k));
          x2=pfin(1,tfin(2,k)); y2=pfin(2,tfin(2,k)); z2=pfin(3,tfin(2,k));
          x3=pfin(1,tfin(3,k)); y3=pfin(2,tfin(3,k)); z3=pfin(3,tfin(3,k));
          nx=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
          ny=(z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
          nz=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
          dummy=sqrt(nx^2+ny^2+nz^2); 
          nx=nx/dummy; ny=ny/dummy; nz=nz/dummy;
          if nx*(x1-xc)+ny*(y1-yc)+nz*(z1-zc)>0
              tfin2(4,k)=nx;
              tfin2(5,k)=ny;
              tfin2(6,k)=nz;
          else
              tfin2(4,k)=-nx;
              tfin2(5,k)=-ny;
              tfin2(6,k)=-nz;
          end
      end
   end
   tfin=tfin2;
end




figure(3)
% Draw surface
for k=1:length(t)
   r1=p3(:,t(1,k)); x1=r1(1); y1=r1(2); z1=r1(3);
   r2=p3(:,t(2,k)); x2=r2(1); y2=r2(2); z2=r2(3);   
   r3=p3(:,t(3,k)); x3=r3(1); y3=r3(2); z3=r3(3);
   xv=[x1 x2 x3 x1];
   yv=[y1 y2 y3 y1];
   zv=[z1 z2 z3 z1];
   h=fill3(xv,yv,zv,'w'); set(h,'linewidth',2)   
   if k==1, hold on; end
end
hold off
axis image

set(gca,'linewidth',1)
xlabel('x','fontsize',24,'FontName','Times New Roman')
ylabel('y','fontsize',24,'FontName','Times New Roman')
zlabel('z','fontsize',24,'FontName','Times New Roman')
set(gca,'fontsize',24,'FontName','Times New Roman')

%save sphere_mesh_R5_240318Pv2exactsurf.m -ascii pfin
%save sphere_mesh_R5_240318Tv2exactsurf.m -ascii tfin


%********************************************************
%* Add points for second order polynomial interpolation *
%********************************************************
if add_points_for_second_order_interpolation == true && ...
        add_element_type_and_details == true
    p=pfin; t=tfin;
    p2=[]; t2=[];
    NP=length(p); NT=length(t);
    CP=1;
    td1=t(1:3,1:length(t));
    td2=t(4:10,1:length(t));
    t=[td1;(td1*0);td2];

    for k=1:length(t)
        x1=p(1,t(1,k)); y1=p(2,t(1,k)); z1=p(3,t(1,k));
        x2=p(1,t(2,k)); y2=p(2,t(2,k)); z2=p(3,t(2,k));
        x3=p(1,t(3,k)); y3=p(2,t(3,k)); z3=p(3,t(3,k));        
    
        x4=(x1+x2)/2; y4=(y1+y2)/2; z4=(z1+z2)/2;
        x5=(x1+x3)/2; y5=(y1+y3)/2; z5=(z1+z3)/2;
        x6=(x2+x3)/2; y6=(y2+y3)/2; z6=(z2+z3)/2;
        
        if t(7,k)==1
            xc=t(8,k); yc=t(9,k); zc=t(10,k);
            R=sqrt((x1-xc)^2+(y1-yc)^2+(z1-zc)^2);
            nx=(x4-xc); ny=(y4-yc); nz=(z4-zc);
            n=sqrt(nx^2+ny^2+nz^2); nx=nx/n; ny=ny/n; nz=nz/n;
            x4=xc+nx*R; y4=yc+ny*R; z4=zc+nz*R;
            nx=(x5-xc); ny=(y5-yc); nz=(z5-zc);
            n=sqrt(nx^2+ny^2+nz^2); nx=nx/n; ny=ny/n; nz=nz/n;
            x5=xc+nx*R; y5=yc+ny*R; z5=zc+nz*R;
            nx=(x6-xc); ny=(y6-yc); nz=(z6-zc);
            n=sqrt(nx^2+ny^2+nz^2); nx=nx/n; ny=ny/n; nz=nz/n;
            x6=xc+nx*R; y6=yc+ny*R; z6=zc+nz*R;
        end
        if k==1
            pd=[x4 x5 x6;y4 y5 y6;z4 z5 z6];
            p=[p pd];
            t(4,k)=length(p)-2; t(5,k)=length(p)-1; t(6,k)=length(p);
        else
            % Check if these points are already found
            found4=false; found5=false; found6=false; i4=0; i5=0; i6=0;
            for j=1:(length(p)-NP)
                if abs(p(1,j+NP)-x4)<1e-7 && abs(p(2,j+NP)-y4)<1e-7 && abs(p(3,j+NP)-z4)<1e-7
                    found4=true;
                    i4=j+NP;
                end
                if abs(p(1,j+NP)-x5)<1e-7 && abs(p(2,j+NP)-y5)<1e-7 && abs(p(3,j+NP)-z5)<1e-7
                    found5=true;
                    i5=j+NP;
                end
                if abs(p(1,j+NP)-x6)<1e-7 && abs(p(2,j+NP)-y6)<1e-7 && abs(p(3,j+NP)-z6)<1e-7
                    found6=true;
                    i6=j+NP;
                end
            end
            if found4==false
               pd=[x4;y4;z4];
               p=[p pd];
               t(4,k)=length(p);
            else
               t(4,k)=i4;
            end
            if found5==false
               pd=[x5;y5;z5];
               p=[p pd];
               t(5,k)=length(p);
            else
               t(5,k)=i5;
            end
            if found6==false
               pd=[x6;y6;z6];
               p=[p pd];
               t(6,k)=length(p);
            else
               t(6,k)=i6;
            end
        end
    end
    
    figure(4)
    % Draw surface
    for k=1:length(t)
       r1=p(:,t(1,k)); x1=r1(1); y1=r1(2); z1=r1(3);
       r2=p(:,t(2,k)); x2=r2(1); y2=r2(2); z2=r2(3);   
       r3=p(:,t(3,k)); x3=r3(1); y3=r3(2); z3=r3(3);
       r4=p(:,t(4,k)); x4=r4(1); y4=r4(2); z4=r4(3);
       r5=p(:,t(5,k)); x5=r5(1); y5=r5(2); z5=r5(3);
       r6=p(:,t(6,k)); x6=r6(1); y6=r6(2); z6=r6(3);
       
       if t(7,k)==1
           % part of spherical surface
           xc=t(8,k); yc=t(9,k); zc=t(10,k);
           R=5;%sqrt((x1-xc)^2+(y1-yc)^2+(z1-zc)^2);
           nx=(x4-xc); ny=(y4-yc); nz=(z4-zc); n=sqrt(nx^2+ny^2+nz^2);
           nx=nx/n; ny=ny/n; nz=nz/n;
           x4=xc+nx*R; y4=yc+ny*R; z4=zc+nz*R;
           nx=(x5-xc); ny=(y5-yc); nz=(z5-zc); n=sqrt(nx^2+ny^2+nz^2);
           nx=nx/n; ny=ny/n; nz=nz/n;
           x5=xc+nx*R; y5=yc+ny*R; z5=zc+nz*R;
           nx=(x6-xc); ny=(y6-yc); nz=(z6-zc); n=sqrt(nx^2+ny^2+nz^2);
           nx=nx/n; ny=ny/n; nz=nz/n;
           x6=xc+nx*R; y6=yc+ny*R; z6=zc+nz*R;
           
           %a few more points are added
           v=0;
           uv=linspace(0,1,11);
           xv1=x1+(x2-x1)*uv; yv1=y1+(y2-y1)*uv; zv1=z1+(z2-z1)*uv;
           xv2=x2+(x3-x2)*uv; yv2=y2+(y3-y2)*uv; zv2=z2+(z3-z2)*uv;
           xv3=x3+(x1-x3)*uv; yv3=y3+(y1-y3)*uv; zv3=z3+(z1-z3)*uv;
           xv=[xv1(1:10) xv2(1:10) xv3]; yv=[yv1(1:10) yv2(1:10) yv3]; zv=[zv1(1:10) zv2(1:10) zv3];
           for ju=1:length(xv)
              nx=xv(ju)-xc; ny=yv(ju)-yc; nz=zv(ju)-zc;
              n=sqrt(nx^2+ny^2+nz^2); nx=nx/n; ny=ny/n; nz=nz/n;
              xv(ju)=xc+nx*R; yv(ju)=yc+ny*R; zv(ju)=zc+nz*R;
           end
           xvd=[x1 x2 x3 x1]; yvd=[y1 y2 y3 y1]; zvd=[z1 z2 z3 z1];
           h=fill3(xvd,yvd,zvd,'w'); set(h,'linewidth',2,'EdgeColor','r')              
           
          h=fill3(xv,yv,zv,'g'); set(h,'linewidth',3)
       end
        % Part of Cylinder
       if t(7,k)==2
          xc=t(8,k); yc=t(9,k); zc=p(3,t(4,k));
          R = a;
          nx=(x4-xc); ny=(y4-yc); nz=(z4-zc); n=sqrt(nx^2+ny^2+nz^2);
          nx=nx/n; ny=ny/n; nz=nz/n;
          x4=xc+nx*R; y4=yc+ny*R; z4=zc+nz*R;
          xc=t(8,k); yc=t(9,k); zc=p(3,t(5,k));
          nx=(x5-xc); ny=(y5-yc); nz=(z5-zc); n=sqrt(nx^2+ny^2+nz^2);
          nx=nx/n; ny=ny/n; nz=nz/n;
          x5=xc+nx*R; y5=yc+ny*R; z5=zc+nz*R;
          xc=t(8,k); yc=t(9,k); zc=p(3,t(6,k));
          nx=(x6-xc); ny=(y6-yc); nz=(z6-zc); n=sqrt(nx^2+ny^2+nz^2);
          nx=nx/n; ny=ny/n; nz=nz/n;
          x6=xc+nx*R; y6=yc+ny*R; z6=zc+nz*R;
%           %a few more points are added to cyl
%            v=0;
%            uv=linspace(0,1,11);
%            xv1=x1+(x2-x1)*uv; yv1=y1+(y2-y1)*uv; zv1=z1+(z2-z1)*uv;
%            xv2=x2+(x3-x2)*uv; yv2=y2+(y3-y2)*uv; zv2=z2+(z3-z2)*uv;
%            xv3=x3+(x1-x3)*uv; yv3=y3+(y1-y3)*uv; zv3=z3+(z1-z3)*uv;
%            xv=[xv1(1:10) xv2(1:10) xv3]; yv=[yv1(1:10) yv2(1:10) yv3]; zv=[zv1(1:10) zv2(1:10) zv3];
%            for ju=1:length(xv)
%               nx=xv(ju)-xc; ny=yv(ju)-yc; nz=zv(ju)-zc;
%               n=sqrt(nx^2+ny^2+nz^2); nx=nx/n; ny=ny/n; nz=nz/n;
%               xv(ju)=xc+nx*R; yv(ju)=yc+ny*R; zv(ju)=zc+nz*R;
%            end
%         
        end
     %            xvd=[x1 x2 x3 x1]; yvd=[y1 y2 y3 y1]; zvd=[z1 z2 z3 z1];
    %       h=fill3(xvd,yvd,zvd,'w'); set(h,'linewidth',2,'EdgeColor','w')              
           
      %     h=fill3(xv,yv,zv,'w'); set(h,'linewidth',3)
         % h=plot3(xv,yv,zv,'k'); set(h,'linewidth',2)   
       if k==1
           hold on;
       end
       h=plot3(x1,y1,z1,'ko'); set(h,'markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
       h=plot3(x2,y2,z2,'ko'); set(h,'markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
       h=plot3(x3,y3,z3,'ko'); set(h,'markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
       h=plot3(x4,y4,z4,'ko'); set(h,'markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
       h=plot3(x5,y5,z5,'ko'); set(h,'markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
       h=plot3(x6,y6,z6,'ko'); set(h,'markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
    end
    hold off
    axis image

    set(gca,'linewidth',1)
    xlabel('x','fontsize',24,'FontName','Times New Roman')
    ylabel('y','fontsize',24,'FontName','Times New Roman')
    zlabel('z','fontsize',24,'FontName','Times New Roman')
    set(gca,'fontsize',24,'FontName','Times New Roman')
    
    %Output files
    save sphere_mesh_R5_280318Pv0exactsurf_quadratic.m -ascii p
    save sphere_mesh_R5_280318Tv0exactsurf_quadratic.m -ascii t
end



