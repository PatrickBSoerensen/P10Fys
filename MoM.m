classdef MoM
    %MOM holds methods for computing current density for an antenna, and
    %the resulting current density given other antennas influence.
    %Additionally it holds methods for computing the radiated E-field
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = MoM()
            %MOM Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function [ant1, area] = mom2on1(obj, ant1, ant2, area, alpha, k, w, thetaI, phi, phiS, mu)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for i=1:ant1.Segments
                for j=1:ant2.Segments
                    if i==j
                        R = @(phimark) sqrt((ant1.Coord(i,3)/4)^.2 ... 
                        +2*ant1.Coord(i,2).^2.*(1-cos(phimark)));
                    else
                        R = @(phimark) sqrt((ant1.Coord(i,1)-ant2.Coord(j,1)).^2 ...
                        +(ant1.Coord(i,2)-ant2.Coord(j,2)).^2 ...
                        +2.*ant1.Coord(i,2).*ant2.Coord(j,2).*(1-cos(phimark)));
                    end
        
                    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func3 = @(phimark) sin(phimark).*sin(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    
                    %Should possibly integrate to 2*pi
                    G1 = ant1.Coord(i,3).*ant2.Coord(j,3).*integral(Func1, 0, pi);
                    G2 = ant1.Coord(i,3).*ant2.Coord(j,3).*integral(Func2, 0, pi);
                    G3 = ant1.Coord(i,3).*ant2.Coord(j,3).*integral(Func3, 0, pi);
                    
                    %Ztt
                    ant1.Z(i,j) = (ant1.T1(i)+ant1.T2(i)).*(ant2.T1(j)+ant2.T2(j)).*(sin(ant1.gamma(i)).*sin(ant2.gamma(j)).*G2 ...  
                        +cos(ant1.gamma(i)).*cos(ant2.gamma(j)).*G1) ...
                        -1./k.^2.*(ant1.T1D(i)+ant1.T2D(i)).*(ant2.T1D(j)+ant2.T2D(j)).*G1;
                    %Ztphi
                    ant1.Z(i+ant1.Segments,j) = 1i*(sin(ant1.gamma(i))*((ant1.T1(i)+ant1.T2(i))*(ant2.T1(j)+ant2.T2(j)))*G3 ...
                        +(1/k^2)*(alpha/ant2.Coord(j,2))*((ant1.T1D(i)+ant1.T2D(i))*(ant2.T1(j)*ant2.T2(j)))*G1);
                    %Zphit
                    ant1.Z(i,j+ant1.Segments) = 1i*(((ant1.T1(i)+ant1.T2(i)).*(ant2.T1(j)+ant2.T2(j)))*sin(ant2.gamma(j))*G3 ... 
                        +(1/k^2)*(alpha/ant1.Coord(i,2))*(((ant1.T1(i)+ant1.T2(i))*(ant2.T1D(j)+ant2.T2D(j))).*G1));
                    %Zphiphi
                    ant1.Z(i+ant1.Segments,j+ant1.Segments) = -(ant1.T1(i)+ant1.T2(i)).*(ant2.T1(j)+ant2.T2(j)).*... 
                        (G2-1/k.^2.*alpha.^2/(ant1.Coord(i,2).*ant2.Coord(j,2)).*G1);
                end
        
                J0 = besselj(alpha-1, k*ant1.Coord(i,2)*sin(thetaI));
                J1 = besselj(alpha, k*ant1.Coord(i,2)*sin(thetaI));
                J2 = besselj(alpha+1, k*ant1.Coord(i,2)*sin(thetaI));
                %% planewave b equations
                if ant1.E0(i) ~= 0
                ant1.btTheta(i) = -ant1.E0(i)*1i/(w*mu)*pi*1i^(alpha)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3)...
                *exp(1i*k*ant1.Coord(i,1)*cos(thetaI))*(cos(thetaI)...
                *sin(ant1.gamma(i))*1i*(J2-J0)-2*sin(thetaI)*cos(ant1.gamma(i))*J1);
    
                ant1.bPhiTheta(i) = ant1.E0(i)*1i/(w*mu)*pi*1i^(alpha)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3)...
                *exp(1i*k*ant1.Coord(i,1)*cos(thetaI))*(cos(thetaI)...
                *(J2+J0));
    
                ant1.btPhi(i) = -ant1.E0(i)*1i/(w*mu)*pi*1i^(alpha)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3)...
                *exp(1i*k*ant1.Coord(i,1)*sin(thetaI))*(sin(ant1.gamma(i))...
                *(J2-J0));
    
                ant1.bPhiPhi(i) = -ant1.E0(i)*1i/(w*mu)*pi*1i^(alpha+1)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3)...
                *exp(1i*k*ant1.Coord(i,1)*sin(thetaI))*(J2-J0);
                else
                %% general b equations, missing hat vectors in integral
                ant1.btTheta(i) = -2*pi*1i/(w*mu)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3);
                ant1.bPhiTheta(i) = -2*pi*1i/(w*mu)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3);
                ant1.btPhi(i) = -2*pi*1i/(w*mu)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3);
                ant1.bPhiPhi(i) = -2*pi*1i/(w*mu)*(ant1.T1(i)+ant1.T2(i))*ant1.Coord(i,3);
                end
            end
            
            ant1.invZ = ant1.Z^(-1);
            bThe = [ant1.btTheta, ant1.bPhiTheta];
            xTheta = ant1.invZ*bThe.';
            bPhi = [ant1.btPhi, ant1.bPhiPhi];
            xPhi = ant1.invZ*bPhi.';
            
            ant1.xtTheta = xTheta(1:ant1.Segments);
            ant1.xPhiTheta = xTheta(ant1.Segments+1:2*ant1.Segments);    
            ant1.xtPhi = xPhi(1:ant1.Segments);
            ant1.xPhiPhi = xPhi(ant1.Segments+1:2*ant1.Segments);
            
            ftn = ant1.tHat(i,:).*(ant1.T1+ant1.T2).*exp(1i.*alpha.*phi)./ant1.Coord(:,2);%T, alpha, n. Expansions function
            fpn = [0 1 0].*(ant1.T1+ant1.T2).*exp(1i.*alpha.*phi)./ant1.Coord(:,2);%Phi, alpha, n. Expansions function
            
            if alpha == 0
                ant1.Jthe=ant1.xtTheta.*ftn;
                ant1.Jphi=ant1.xtPhi.*fpn;
            else
                ant1.Jthe = ant1.Jthe+2*(ant1.xtTheta.*ftn.*cos(alpha.*phi)+1i*ant1.xPhiTheta.*fpn.*sin(alpha.*phi));
                ant1.Jphi = ant1.Jphi+2*(1i*ant1.xtPhi.*ftn.*sin(alpha.*phi)+ant1.xPhiPhi.*fpn.*cos(alpha.*phi));
            end
            ant1.Jthe(1,1) = 0;
            ant1.Jthe(end,1) = 0;
            ant1.Jphi(1,1) = 0;
            ant1.Jphi(end,1) = 0;
            
            area = emission(obj, ant1, area, alpha, k, w, phiS);
        end
        
        function [ant, area] = mombasis(obj, ant, area, alpha, k, w, thetaI, phi, phiS, mu)
            DM = eye(ant.Segments);
            DMO = circshift(DM,1,2);
            for i=1:length(ant.T)-1
                iSegments = [i, i+1];
                for j=1:length(ant.T)-1
                    jSegments = [j, j+1];
                    if i==j
                    R = @(phimark) sqrt((ant.Coord(iSegments,3)/4).^2 ... 
                        +2*ant.Coord(iSegments,2).^2.*(1-cos(phimark)));
                    else
                    R = @(phimark) sqrt((ant.Coord(iSegments,1)-ant.Coord(jSegments,1)).^2 ...
                        +(ant.Coord(iSegments,2)-ant.Coord(jSegments,2)).^2 ...
                        +2.*ant.Coord(iSegments,2).*ant.Coord(jSegments,2).*(1-cos(phimark)));
                    end
                    
                    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func3 = @(phimark) sin(phimark).*sin(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    
                    %Should possibly integrate to 2*pi
                    G1 = ant.Coord(i,3).*ant.Coord(i+1,3).*integral(Func1, 0, pi, 'ArrayValued', true);
                    G2 = ant.Coord(i,3).*ant.Coord(i+1,3).*integral(Func2, 0, pi, 'ArrayValued', true);
                    G3 = ant.Coord(i,3).*ant.Coord(i+1,3).*integral(Func3, 0, pi, 'ArrayValued', true);
                
                    %Ztt
                    ant.Z(i,j) = sum(((ant.T(iSegments).*ant.T(jSegments)).*(sin(ant.gamma(iSegments)).*sin(ant.gamma(jSegments)).*G2 ...
                        +cos(ant.gamma(iSegments)).*cos(ant.gamma(jSegments)).*G1) ...
                        -1./k.^2.*ant.TD(iSegments).*ant.TD(jSegments).*G1));
                    %Ztphi
                    ant.Z(i+ant.Segments,j) = sum(1i.*(sin(ant.gamma(iSegments)).*ant.T(iSegments).*ant.T(jSegments).*G3 ...
                        +(1./k^2).*(alpha./ant.Coord(jSegments,2)).*ant.TD(iSegments).*ant.T(jSegments).*G1));
                    %Zphit
                    ant.Z(i,j+ant.Segments) = sum(1i.*ant.T(iSegments).*ant.T(jSegments).*sin(ant.gamma(jSegments)).*G3 ... 
                        +(1./k^2).*(alpha./ant.Coord(iSegments,2)).*ant.T(iSegments).*ant.TD(jSegments).*G1);
                    %Zphiphi
                    ant.Z(i+ant.Segments,j+ant.Segments) = sum(-(ant.T(iSegments).*ant.T(jSegments).*... 
                        (G2-1/k.^2.*alpha.^2./(ant.Coord(iSegments,2).*ant.Coord(jSegments,2)).*G1)));
                end
                
                J0 = besselj(alpha-1, k*ant.Coord(i,2)*sin(thetaI));
                J1 = besselj(alpha, k*ant.Coord(i,2)*sin(thetaI));
                J2 = besselj(alpha+1, k*ant.Coord(i,2)*sin(thetaI));
                %% planewave b equations
                if ant.E0(i) ~= 0
                ant.btTheta(i) = -ant.E0(i)*1i/(w*mu)*pi*1i^(alpha)*(ant.T(i))*ant.Coord(i,3)...
                *exp(1i*k*ant.Coord(i,1)*cos(thetaI))*(cos(thetaI)...
                *sin(ant.gamma(i))*1i*(J2-J0)-2*sin(thetaI)*cos(ant.gamma(i))*J1);
    
                ant.bPhiTheta(i) = ant.E0(i)*1i/(w*mu)*pi*1i^(alpha)*(ant.T(i))*ant.Coord(i,3)...
                *exp(1i*k*ant.Coord(i,1)*cos(thetaI))*(cos(thetaI)...
                *(J2+J0));
    
                ant.btPhi(i) = -ant.E0(i)*1i/(w*mu)*pi*1i^(alpha)*(ant.T(i))*ant.Coord(i,3)...
                *exp(1i*k*ant.Coord(i,1)*sin(thetaI))*(sin(ant.gamma(i))...
                *(J2-J0));
    
                ant.bPhiPhi(i) = -ant.E0(i)*1i/(w*mu)*pi*1i^(alpha+1)*(ant.T(i))*ant.Coord(i,3)...
                *exp(1i*k*ant.Coord(i,1)*sin(thetaI))*(J2-J0);
            
                else
                %% general b equations, missing hat vectors in integral
                ant.btTheta(i) = -2*pi*1i/(w*mu)*(ant.T(i))*ant.Coord(i,3);
                ant.bPhiTheta(i) = -2*pi*1i/(w*mu)*(ant.T(i))*ant.Coord(i,3);
                ant.btPhi(i) = -2*pi*1i/(w*mu)*(ant.T(i))*ant.Coord(i,3);
                ant.bPhiPhi(i) = -2*pi*1i/(w*mu)*(ant.T(i))*ant.Coord(i,3);
                end
            end
            
            ant.invZ = ant.Z^(-1);
            bThe = [ant.btTheta, ant.bPhiTheta];
            xTheta = ant.invZ*bThe.';
            bPhi = [ant.btPhi, ant.bPhiPhi];
            xPhi = ant.invZ*bPhi.';
            
            ant.xtTheta = xTheta(1:ant.Segments);
            ant.xPhiTheta = xTheta(ant.Segments+1:2*ant.Segments);    
            ant.xtPhi = xPhi(1:ant.Segments);
            ant.xPhiPhi = xPhi(ant.Segments+1:2*ant.Segments);
            
            ftn = ant.tHat(i,:).*(ant.T).*exp(1i.*alpha.*phi)./ant.Coord(:,2);%T, alpha, n. Expansions function
            fpn = [0 1 0].*(ant.T1+ant.T2).*exp(1i.*alpha.*phi)./ant.Coord(:,2);%Phi, alpha, n. Expansions function
            
            if alpha == 0
                ant.Jthe=ant.xtTheta.*ftn;
                ant.Jphi=ant.xtPhi.*fpn;
            else
                ant.Jthe = ant.Jthe+2*(ant.xtTheta.*ftn.*cos(alpha.*phi)+1i*ant.xPhiTheta.*fpn.*sin(alpha.*phi));
                ant.Jphi = ant.Jphi+2*(1i*ant.xtPhi.*ftn.*sin(alpha.*phi)+ant.xPhiPhi.*fpn.*cos(alpha.*phi));
            end
            ant.Jthe(1,1) = 0;
            ant.Jthe(end,1) = 0;
            ant.Jphi(1,1) = 0;
            ant.Jphi(end,1) = 0;
            
            area = emission(obj, ant, area, alpha, k, w, phiS);
        end
       
        
        function area = emission(obj, ant, area, alpha, k, w, phiS)
            rz = (area.z-ant.Coord(:,1));
            for i=1:ant.Segments
                rx = (area.x+ant.Coord(:,2)+area.SingularityProtection);
                r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                B = -(1i*w*area.mu0)/(2*pi)*(exp(-1i*k*r)./r);
                if alpha == 0
                    area.Ethethe = B/2 .* ant.xtTheta(i) .* ant.btTheta(i) + area.Ethethe;
                    area.Ephiphi = B/2 .* ant.xtPhi(i) .* ant.btPhi(i) + area.Ephiphi;
                else    
                    area.Ethethe = B*(ant.xtTheta(i)*ant.btTheta(i)...
                        +ant.xPhiTheta(i)*ant.bPhiTheta(i))*cos(alpha*phiS)...
                        +area.Ethethe;
                
                    area.Ephithe = 1i*B*(ant.xtTheta(i)*ant.btPhi(i)...
                        +ant.xPhiTheta(i)*ant.bPhiPhi(i))*sin(alpha*phiS)...
                        +area.Ephithe;
                
                    area.Ethephi = 1i*B*(ant.xtPhi(i)*ant.btTheta(i)...
                        +ant.xPhiPhi(i)*ant.bPhiTheta(i))*sin(alpha*phiS)...
                        +area.Ethephi;
                
                    area.Ephiphi = B*(ant.xtPhi(i)*ant.btPhi(i)...
                        +ant.xPhiPhi(i)*ant.bPhiPhi(i))*cos(alpha*phiS)...
                        +area.Ephiphi;
                end
        
                rx = (area.x-ant.Coord(:,2)-area.SingularityProtection);
                r = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                B = -(1i*w*area.mu0)/(2*pi)*(exp(-1i*k*r)./r);
                if alpha == 0
                    area.Ethethe = B/2 .* ant.xtTheta(i) .* ant.btTheta(i) + area.Ethethe;
                    area.Ephiphi = B/2 .* ant.xtPhi(i) .* ant.btPhi(i) + area.Ephiphi;
                else
                    area.Ethethe = B*(ant.xtTheta(i)*ant.btTheta(i)...
                        +ant.xPhiTheta(i)*ant.bPhiTheta(i))*cos(alpha*phiS)...
                        +area.Ethethe;
                
                    area.Ephithe = 1i*B*(ant.xtTheta(i)*ant.btPhi(i)...
                        +ant.xPhiTheta(i)*ant.bPhiPhi(i))*sin(alpha*phiS)...
                        +area.Ephithe;
                
                    area.Ethephi = 1i*B*(ant.xtPhi(i)*ant.btTheta(i)...
                        +ant.xPhiPhi(i)*ant.bPhiTheta(i))*sin(alpha*phiS)...
                        +area.Ethephi;
                
                    area.Ephiphi = B*(ant.xtPhi(i)*ant.btPhi(i)...
                        +ant.xPhiPhi(i)*ant.bPhiPhi(i))*cos(alpha*phiS)...
                        +area.Ephiphi;
                end
            end
        end
        
    end
end

