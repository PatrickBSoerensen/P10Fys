classdef Area
    %AREA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SingularityProtection;
        x;
        z;
        Ethethe;
        mu0;
    end
    
    methods
        function obj = Area(SingularityProtection, xsteps, zsteps, xmin, xmax, zmin, zmax, mu0)
            %AREA Construct an instance of this class
            %   Detailed explanation goes here
            obj.SingularityProtection = SingularityProtection;
            obj.x = linspace(xmin, xmax, xsteps);
            obj.z = linspace(zmin, zmax, zsteps);
            
            obj.Ethethe = 0;
            
            obj.mu0 = mu0;
        end
        
        function ThreeDplot(obj, opl1, opl2, ant, k, str)
            %% Far Field 3d
        phi=linspace(0,2*pi,opl1);
        theta=linspace(0,pi,opl2);
        E=zeros(opl1,opl2);
        for j=1:opl1
            for j12=1:opl2
                    v=phi(j);
                    v2=theta(j12);
                    Ry=ant.Centre(1);
                    r=@(z) z.*cos(v2)+Ry*sin(v2)*sin(v);
                    konst=exp(1i*k*str)/(4*pi*str);
                    g=@(z) exp(-1i.*k.*r(z));
                    E(j,j12)=E(j,j12)+sum(-sin(v2).*konst.*g(ant.CoordTest(1:end-1,1)).*ant.Jthe.*ant.CoordTest(1:end-1,3));
                    E(j,j12)=E(j,j12)+sum(-sin(v2).*konst.*g(ant.CoordTest(2:end,1)).*ant.Jthe.*ant.CoordTest(2:end,3));
            end
        end
        figure(2)
        patternCustom(abs(E)./(max(max(abs(E)))),360.*theta./(2*pi),(360.*phi./(2*pi))')
        title('Far Field 3D plot');
        end
        
        function farfieldplotXZ(obj, ant, distance, k, w, nr)
        opl = length(obj.x);
        theta=linspace(0,2*pi,opl);
        E=zeros(1,opl);
        rz = (sin(theta)*distance-ant.CoordTest(:,1));
        rx = (sin(theta)*distance-(ant.CoordTest(:,2)+obj.SingularityProtection));
            for i=1:length(ant.T1)
                r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                
                B = -(1i.*w.*obj.mu0);
                g = B.*exp(-1i.*k.*r)./(4.*pi.*r);
                G = g.*((1+1i./(k*r)-1./((k*r).^2)) - ...
                            ((rx(i,:)).^2)./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
                E = E+G(i,:).*ant.Jthe(i).*ant.CoordTest(i,3);
            end
%             for j=1:opl
%             r = ant.CoordTest(:,1).*cos(theta(j))+ant.CoordTest(:,2)*sin(theta(j));
%             B = -(1i.*w.*obj.mu0);
%             g = B.*exp(-1i.*k.*r)./(4.*pi.*r);
%             G = g.*((1+1i./(k*r)-1./((k*r).^2)) - ...
%                 ((r).^2)./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
%             E(j)=E(j)+sum(G.*ant.Jthe.*ant.CoordTest(:,3));
%             end
        figure(nr)
        polarplot(theta,abs(E)./max(abs(E)))
        title('Far field for a Half-wave dipole for $\theta$');
        end
        
        function farfieldplotXY(obj, ant, distance, k)
            opl = length(obj.x);
            phi=linspace(0,2*pi,opl);
            theta=pi/2;
            E=zeros(1,opl);
            for j=1:opl
                    r= sqrt((ant.CoordTest(:,1)).^2+(ant.CoordTest(:,2)).^2);
                    konst=exp(1i*k*distance)/(4*pi*distance);
                    g= exp(-1i.*k.*r);
                    E(j)=E(j)+sum(-sin(theta).*konst.*g(1:end-1).*ant.Jthe.*ant.CoordTest(1:end-1,3));
                    E(j)=E(j)+sum(-sin(theta).*konst.*g(2:end).*ant.Jthe.*ant.CoordTest(2:end,3));
            end
            figure(9)
            polarplot(phi,abs(E)./max(abs(E)))
            title('Far field for a Half-wave dipole for $\phi$');
        end
    end
end

