% Code for calculating the indirect Greens tensor   
    k0 = 2;
    strip_length = 10;
    strip_width = 0.5;
    dx = 0.1;
    Nz = 10;
    Nx = 10;
    Ny = 10;
    epsL2 = 3.4;
    eps1 = 1;
    epsL3 = 2;
    layer_thickness = 20;
    k1i2=1;
    
    kz1f=@(krho) sqrt(k0^2*eps1-krho.*krho);
    kzL2f=@(krho) sqrt(k0^2*epsL2-krho.*krho);
    kzL3f=@(krho) sqrt(k0^2*epsL3-krho.*krho);
    
    rp1L2f=@(krho) (epsL2*kz1f(krho)-eps1*kzL2f(krho))./(epsL2*kz1f(krho)+eps1*kzL2f(krho));
    rs1L2f=@(krho) (kz1f(krho)-kzL2f(krho))./(kz1f(krho)+kzL2f(krho));    
    rpL2L3f=@(krho) (epsL3*kzL2f(krho)-epsL2*kzL3f(krho))./(epsL3*kzL2f(krho)+epsL2*kzL3f(krho));
    rsL2L3f=@(krho) (kzL2f(krho)-kzL3f(krho))./(kzL2f(krho)+kzL3f(krho));    
    rpf=@(krho) (rp1L2f(krho)+rpL2L3f(krho).*exp(2i*kzL2f(krho)*layer_thickness))...
        ./(1+rp1L2f(krho).*rpL2L3f(krho).*exp(2i*kzL2f(krho)*layer_thickness));
    rsf=@(krho) (rs1L2f(krho)+rsL2L3f(krho).*exp(2i*kzL2f(krho)*layer_thickness))...
        ./(1+rs1L2f(krho).*rsL2L3f(krho).*exp(2i*kzL2f(krho)*layer_thickness));
    
    ellipse_length = k0*5;
    ellipse_height = k0*0.2;
    krhof=@(alpha) (1+cos(alpha))*ellipse_length/2+1i*ellipse_height*sin(alpha);
    dkrho_dalphaf=@(alpha) -sin(alpha)*ellipse_length/2+1i*ellipse_height*cos(alpha);
    rhomin=0; rhomax=sqrt(strip_length^2+strip_width^2);
    rhov1=[0 dx dx*sqrt(2) dx*2 dx*sqrt(5) dx*sqrt(8)];
    if rhomax>dx*3
        rhov2=linspace(dx*3,rhomax,ceil((rhomax-dx*3)/dx)+1);
    else
        rhov2=[];
    end
    rhotabv=[rhov1 rhov2];
    Nrho=length(rhotabv);
    
    GItabzz=zeros(Nrho,2*Nz-1); GItabzr=GItabzz; 
    GItabrr=GItabzz; GItabpp=GItabzz;
    tab_z=GItabzz; tab_r=GItabzz;
     for jrho=1:Nrho
         rho=rhotabv(jrho);
         for jz=1:Nz*2-1
            z=jz*dx;
                
            dJmf=@(krho) -besselj(1,krho*rho)./(krho*rho);
            if jrho==1
                dJmf=@(krho) -0.5;
            end
            Jmmf=@(krho) -0.5*(besselj(0,krho*rho)-besselj(2,krho*rho));
            Gizzf=@(krho) 1i/(4*pi)*rpf(krho).*besselj(0,krho*rho).*(krho.^2)/k1i2.*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
            Gizrf=@(krho) 1/(4*pi)*rpf(krho).*krho.*kz1f(krho)/k1i2.*(-besselj(1,krho*rho)).*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
            Gippf=@(krho) 1i/(4*pi)*(rpf(krho).*dJmf(krho).*(kz1f(krho).^2)/k1i2-rsf(krho).*Jmmf(krho)).*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
            Girrf=@(krho) 1i/(4*pi)*(rpf(krho).*Jmmf(krho).*(kz1f(krho).^2)/k1i2-rsf(krho).*dJmf(krho)).*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);

            Gzz=0; Gzr=0; Grr=0; Gpp=0;  
            % Integral going into the complex plane to avoid poles
            Gzz=quadgk(@(alpha) Gizzf(krhof(alpha)).*dkrho_dalphaf(alpha),-pi,0,'MaxIntervalCount',10000);
            Gzr=quadgk(@(alpha) Gizrf(krhof(alpha)).*dkrho_dalphaf(alpha),-pi,0,'MaxIntervalCount',10000);
            Gpp=quadgk(@(alpha) Gippf(krhof(alpha)).*dkrho_dalphaf(alpha),-pi,0,'MaxIntervalCount',10000);
            Grr=quadgk(@(alpha) Girrf(krhof(alpha)).*dkrho_dalphaf(alpha),-pi,0,'MaxIntervalCount',10000);
                
            % The rest of the integrals
            Gzz=Gzz+quadgk(@(krho) Gizzf(krho),ellipse_length,inf,'MaxIntervalCount',10000);
            Gzr=Gzr+quadgk(@(krho) Gizrf(krho),ellipse_length,inf,'MaxIntervalCount',10000);
            Gpp=Gpp+quadgk(@(krho) Gippf(krho),ellipse_length,inf,'MaxIntervalCount',10000);
            Grr=Grr+quadgk(@(krho) Girrf(krho),ellipse_length,inf,'MaxIntervalCount',10000);                

            GItabzz(jrho,jz)=Gzz;
            GItabzr(jrho,jz)=Gzr;
            GItabpp(jrho,jz)=Gpp;
            GItabrr(jrho,jz)=Grr;
            tab_z(jrho,jz)=z;
            tab_r(jrho,jz)=rho;
         end
     end

     % Conversion to cartesian components
     
     for jx=-(Nx-1):(Nx-1)
        x=jx*dx
        for jy=-(Ny-1):(Ny-1)
            y=jy*dx;
            for jz=1:Nz*2-1
                z=jz*dx;

                rho=sqrt(x*x+y*y);
                phi=atan2(y,x);

               Gzz=interp2(tab_z,tab_r,GItabzz,z,rho,'spline'); 
               Gzr=interp2(tab_z,tab_r,GItabzr,z,rho,'spline'); 
               Grr=interp2(tab_z,tab_r,GItabrr,z,rho,'spline'); 
               Gpp=interp2(tab_z,tab_r,GItabpp,z,rho,'spline');                
               Grz=-Gzr;
               
               GIxx(jx+Nx,jy+Ny,jz)=((sin(phi))^2)*Gpp+((cos(phi))^2)*Grr;
               GIxy(jx+Nx,jy+Ny,jz)=sin(phi)*cos(phi)*(Grr-Gpp);
               GIxz(jx+Nx,jy+Ny,jz)=cos(phi)*Grz;
               GIyx(jx+Nx,jy+Ny,jz)=sin(phi)*cos(phi)*(Grr-Gpp);
               GIyy(jx+Nx,jy+Ny,jz)=((cos(phi))^2)*Gpp+((sin(phi))^2)*Grr;
               GIyz(jx+Nx,jy+Ny,jz)=sin(phi)*Grz;
               GIzx(jx+Nx,jy+Ny,jz)=cos(phi)*Gzr;
               GIzy(jx+Nx,jy+Ny,jz)=sin(phi)*Gzr;
               GIzz(jx+Nx,jy+Ny,jz)=Gzz;
            end
        end
     end
     
    
     