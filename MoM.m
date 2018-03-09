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
            %mom2on1 Calculates the influence of antenna 2 on antenna 1
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
            
            ftn = sqrt(ant.tHat(1,:).^2+ant.tHat(3,:).^2).*(ant1.T1+ant1.T2).*exp(1i.*alpha.*phi)./ant1.Coord(:,2);%T, alpha, n. Expansions function
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
            for i=1:length(ant.T1)
                iSegments = [i, i+1, i, i+1];
                for j=1:length(ant.T1)
                    jSegments = [j, j+1, j+1, j];
                    for h=1:4
                        if iSegments(1,h) == jSegments(1,h)
                            R{1,h} = @(phimark) sqrt((ant.CoordTest(iSegments(1,h),3)/4).^2 ...
                            +2*ant.CoordTest(iSegments(1,h),2).^2.*(1-cos(phimark)));
                        else
                            R{1,h} = @(phimark) sqrt((ant.CoordTest(iSegments(1,h),1)-ant.CoordTest(jSegments(1,h),1)).^2 ...
                            +(ant.CoordTest(iSegments(1,h),2)-ant.CoordTest(jSegments(1,h),2)).^2 ...
                            +2.*ant.CoordTest(iSegments(1,h),2).*ant.CoordTest(jSegments(1,h),2).*(1-cos(phimark)));
                        end
              
                    Func1{1,h} = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R{1,h}(phimark))./R{1,h}(phimark);
                    Func2{1,h} = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R{1,h}(phimark))./R{1,h}(phimark);
                    
                    %Should possibly integrate to 2*pi
                    G1{1,h} = ant.CoordTest(iSegments(1,h),3).*ant.CoordTest(jSegments(1,h),3).*integral(Func1{1,h}, 0, pi, 'ArrayValued', true);
                    G2{1,h} = ant.CoordTest(iSegments(1,h),3).*ant.CoordTest(jSegments(1,h),3).*integral(Func2{1,h}, 0, pi, 'ArrayValued', true);
                    end
                    
                    %Ztt
                    ant.Z(i,j) =...
                        ant.T1(i).*ant.T1(j).*...
                        (sin(ant.gammaTest(i)).*sin(ant.gammaTest(j)).*G2{1,1} ...
                        +cos(ant.gammaTest(i)).*cos(ant.gammaTest(j)).*G1{1,1}) ...
                        -1./k.^2.*ant.T1D(i).*ant.T1D(j).*G1{1,1}(1)...
                        +...
                        ant.T1(i).*ant.T2(j).*...
                        (sin(ant.gammaTest(i)).*sin(ant.gammaTest(j+1)).*G2{1,3} ...
                        +cos(ant.gammaTest(i)).*cos(ant.gammaTest(j+1)).*G1{1,3}) ...
                        -1./k.^2.*ant.T1D(i).*ant.T2D(j).*G1{1,3}...
                        +...
                        ant.T2(i).*ant.T2(j).*...
                        (sin(ant.gammaTest(i+1)).*sin(ant.gammaTest(j+1)).*G2{1,2} ...
                        +cos(ant.gammaTest(i+1)).*cos(ant.gammaTest(j+1)).*G1{1,2}) ...
                        -1./k.^2.*ant.T2D(i).*ant.T2D(j).*G1{1,2}(1)...
                        +...
                        ant.T2(i).*ant.T1(j).*...
                        (sin(ant.gammaTest(i+1)).*sin(ant.gammaTest(j)).*G2{1,4} ...
                        +cos(ant.gammaTest(i+1)).*cos(ant.gammaTest(j)).*G1{1,4}) ...
                        -1./k.^2.*ant.T2D(i).*ant.T1D(j).*G1{1,4};
                end
                
                J0 = besselj(alpha-1, k*ant.CoordTest([i, i+1],2)*sin(thetaI));
                J1 = besselj(alpha, k*ant.CoordTest([i, i+1],2)*sin(thetaI));
                J2 = besselj(alpha+1, k*ant.CoordTest([i, i+1],2)*sin(thetaI));
                %% planewave b equations
                if ant.E0(i) ~= 0
                    ant.btTheta(i) = -ant.E0(i)*1i/(w*mu)*pi*1i^(alpha)*...
                    (ant.T1(i).*ant.CoordTest(i,3)...
                    .*exp(1i.*k.*ant.CoordTest(i,1).*cos(thetaI)).*(cos(thetaI)...
                    .*sin(ant.gammaTest(i)).*1i.*(J2(1)-J0(1))-2.*sin(thetaI).*cos(ant.gammaTest(i)).*J1(1))...
                    +ant.T2(i).*ant.CoordTest(i+1,3)...
                    .*exp(1i.*k.*ant.CoordTest(i+1,1).*cos(thetaI)).*(cos(thetaI)...
                    .*sin(ant.gammaTest(i+1)).*1i.*(J2(2)-J0(2))-2.*sin(thetaI).*cos(ant.gammaTest(i+1)).*J1(2)));
                else
                %% general b equations, missing hat vectors in integral
                    ant.btTheta(i) = -2*pi*1i^(alpha)*1i/(w*mu).*(ant.T1(i).*ant.CoordTest(i,3)+ant.T2(i).*ant.CoordTest(i+1,3));

                end
            end
            
            ant.btTheta(1,1)=0;
            ant.btTheta(end,1)=0;
            
            ant.invZ = ant.Z^(-1);
            
            ant.xtTheta = ant.invZ*ant.btTheta.';
            
            tHatLen = sqrt(ant.tHatTest(1:end,1).^2+ant.tHatTest(1:end,3).^2);
            ftn = tHatLen(1:end-1).*ant.T1(:)./ant.CoordTest(1:end-1,2)...
                +tHatLen(2:end).*ant.T2(:)./ant.CoordTest(2:end,2);%T, alpha, n. Expansions function
            
            if alpha == 0
                ant.Jthe = ant.xtTheta.*ftn;
            else
                ant.Jthe = ant.Jthe+2*ant.xtTheta.*ftn.*cos(alpha.*phi);
            end
            
            ant.Jthe(1,1) = 0;
            ant.Jthe(end,1) = 0;
            area = emission(obj, ant, area, alpha, k, w, phiS);
        end
       
        function area = emissionNew(obj, ant, area, alpha, k, w, phiS)
                rz = (area.z-ant.CoordTest(:,1));
                rx = (area.x+ant.CoordTest(:,2)+area.SingularityProtection);
                r = sqrt((rz).^2+(rx).^2);
                B = -(1i*w*area.mu0)/(2*pi)*(exp(-1i*k*r)./r);
                if alpha == 0
                    area.Ethethe = B(1:end-1,:)/2 .* ant.xtTheta .* ant.btTheta.'...
                        +B(2:end,:)/2 .* ant.xtTheta .* ant.btTheta.'...
                        + area.Ethethe;
                else    
                    area.Ethethe = B(1:end-1,:).*ant.xtTheta.*ant.btTheta.'...
                                +B(2:end,:) .* ant.xtTheta .* ant.btTheta.'...
                                .*cos(alpha*phiS)+area.Ethethe;
                end
        
                rx = (area.x-ant.CoordTest(:,2)-area.SingularityProtection);
                r = sqrt((rz).^2+(rx).^2);
                B = -(1i*w*area.mu0)/(2*pi)*(exp(-1i*k*r)./r);
                if alpha == 0
                    area.Ethethe = B(1:end-1,:)/2 .* ant.xtTheta .* ant.btTheta.'...
                        +B(2:end,:)/2 .* ant.xtTheta .* ant.btTheta.'...
                        + area.Ethethe;
                else    
                    area.Ethethe = B(1:end-1,:).*ant.xtTheta.*ant.btTheta.'...
                                +B(2:end,:) .* ant.xtTheta .* ant.btTheta.'...
                                .*cos(alpha*phiS)+area.Ethethe;
                end
        end
        
        function area = emission(obj, ant, area, alpha, k, w, phiS)
            rz = (area.z-ant.CoordTest(:,1));
            for i=1:length(ant.T1)
                rx = (area.x+ant.CoordTest(:,2)+area.SingularityProtection);
                r1 = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                r2 = sqrt((rz(i+1,:).').^2+(rx(i+1,:)).^2);
                
                B1 = -((1i.*w.*area.mu0)./(2.*pi.*r1)).*(exp(-1i.*k.*r1));
                B2 = -((1i.*w.*area.mu0)./(2.*pi.*r2)).*(exp(-1i.*k.*r2));
                if alpha == 0
                    area.Ethethe = B1/2 .* ant.xtTheta(i) .* ant.btTheta(i) ...
                                   +B2/2 .* ant.xtTheta(i) .* ant.btTheta(i) + area.Ethethe;
                else    
                    area.Ethethe = B1.*ant.xtTheta(i).*ant.btTheta(i)...
                        +B2.*ant.xtTheta(i).*ant.btTheta(i)...
                        +area.Ethethe;
                end
        
                rx = (area.x-ant.CoordTest(:,2)-area.SingularityProtection);
                
                r1 = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                r2 = sqrt((rz(i+1,:).').^2+(rx(i+1,:)).^2);
                
                B1 = -((1i.*w.*area.mu0)./(2.*pi.*r1)).*(exp(-1i.*k.*r1));
                B2 = -((1i.*w.*area.mu0)./(2.*pi.*r2)).*(exp(-1i.*k.*r2));
                if alpha == 0
                    area.Ethethe = B1/2 .* ant.xtTheta(i) .* ant.btTheta(i) ...
                                   +B2/2 .* ant.xtTheta(i) .* ant.btTheta(i) + area.Ethethe;
                else
                    area.Ethethe = B1.*ant.xtTheta(i).*ant.btTheta(i)...
                        +B2.*ant.xtTheta(i).*ant.btTheta(i)...
                        +area.Ethethe;
                end
            end
        end
    end
end

