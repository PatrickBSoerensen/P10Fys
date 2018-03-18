classdef MoM
    %MOM holds methods for computing current density for an antenna, and
    %the resulting current density given other antennas influence.
    %Additionally it holds methods for computing the radiated E-field
    %   Detailed explanation goes here
    
    properties
        G1
        G2
    end
    
    methods
        function obj = MoM(ant)
            %MOM Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function [obj, ant, area] = mombasis(obj, ant, area, alpha, k, w, thetaI, phi, phiS, mu)
if alpha == 0
                iSegments = 1:length(ant.T1);
                jSegments = 1:length(ant.T2);
            
                G1 = cell([length(ant.T1) length(ant.T1)]);
                G2 = cell([length(ant.T1) length(ant.T1)]);
                
                for i = 1:length(ant.T1)
                    for h = 1:length(ant.T2)
                        if iSegments(i) == jSegments(h)
                            R = @(phimark) sqrt((ant.CoordTest(iSegments(i),3)/4).^2 ...
                            +2*ant.CoordTest(iSegments(i),2).^2.*(1-cos(phimark)));
                        else
                            R = @(phimark) sqrt((ant.CoordTest(iSegments(i),1)-ant.CoordTest(jSegments(h),1)).^2 ...
                            +(ant.CoordTest(iSegments(i),2)-ant.CoordTest(jSegments(h),2)).^2 ...
                            +2.*ant.CoordTest(iSegments(i),2).*ant.CoordTest(jSegments(h),2).*(1-cos(phimark)));
                        end
                        
                    Func1 = @(phimark, alpha) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func2 = @(phimark, alpha) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    
                    %Should possibly integrate to 2*pi
                    obj.G1{i,h} = @(alpha) ant.CoordTest(iSegments(i),3).*ant.CoordTest(jSegments(h),3).*integral(@(phimark)Func1(phimark,alpha), 0, pi);%, 'ArrayValued', true);
                    obj.G2{i,h} = @(alpha) ant.CoordTest(iSegments(i),3).*ant.CoordTest(jSegments(h),3).*integral(@(phimark)Func2(phimark,alpha), 0, pi);%, 'ArrayValued', true);
                    end
                end
           clear Func1; clear Func2;
end
                for i = 1:length(ant.T1)
                    for h = 1:length(ant.T2)
                    G1{i,h} = obj.G1{i,h}(alpha);
                    G2{i,h} = obj.G2{i,h}(alpha);
                    end
                end
                G1 = cell2mat(G1);
                G2 = cell2mat(G2);
                
            Z1 = ant.T1*ant.T1.'.*...
                (sin(ant.gammaTest)*sin(ant.gammaTest.').*G2  ...
                +cos(ant.gammaTest)*cos(ant.gammaTest.').*G1)...
                -(1./k.^2)*ant.T1D*ant.T1D.'.*G1;
            
            Z2 = ant.T1*ant.T2.'.*...
                (sin(ant.gammaTest)*sin(ant.gammaTest.').*G2...
                +cos(ant.gammaTest)*cos(ant.gammaTest.').*G1)...
                -(1./k.^2)*ant.T1D*ant.T2D.'.*G1;    
            
            Z3 = ant.T2*ant.T2.'.*...
                (sin(ant.gammaTest)*sin(ant.gammaTest.').*G2 ...
                +cos(ant.gammaTest)*cos(ant.gammaTest.').*G1) ...
                -(1./k.^2)*ant.T2D*ant.T2D.'.*G1;
            
            Z4 = ant.T2*ant.T1.'.*...
                (sin(ant.gammaTest)*sin(ant.gammaTest.').*G2 ...
                +cos(ant.gammaTest)*cos(ant.gammaTest.').*G1) ...
                -(1./k.^2)*ant.T2D*ant.T1D.'.*G1;

            ant.Z = Z1+Z2+Z3+Z4;
            
            J0 = besselj(alpha-1, k*ant.CoordTest(:,2)*sin(thetaI));
            J1 = besselj(alpha, k*ant.CoordTest(:,2)*sin(thetaI));
            J2 = besselj(alpha+1, k*ant.CoordTest(:,2)*sin(thetaI));
            
            %% planewave b equations
%             maybe ant.E0.* front faktor
%             b1= pi.*1i.^(alpha).*...
%                 (ant.T1.*ant.CoordTest(:,3)...
%                 .*exp(1i.*k.*ant.CoordTest(:,1).*cos(thetaI)).*...
%                 (cos(thetaI).*sin(ant.gammaTest).*1i.*(J2-J0)...
%                 -2.*sin(thetaI).*cos(ant.gammaTest).*J1));
%             
%             b2= pi.*1i.^(alpha).*...
%                 (ant.T2.*ant.CoordTest(:,3)...
%                 .*exp(1i.*k.*ant.CoordTest(:,1).*cos(thetaI)).*...
%                 (cos(thetaI).*sin(ant.gammaTest).*1i.*(J2-J0)...
%                 -2.*sin(thetaI).*cos(ant.gammaTest).*J1));
%             
%             ant.btTheta = b1+b2;
%%
%%General b expression for arbritary wave expression. Specified as tHat dot thetaHat direction          
           E=@(phi) ant.E0.*exp(1i.*k.*ant.CoordTest(:,1).*cos(thetaI)).*exp(1i.*k.*sin(thetaI).*ant.CoordTest(:,3).*cos(phi));           
           ant.btTheta=(ant.T1+ant.T2).*ant.CoordTest(:,3).*...
               integral(@(phi)exp(-1i.*alpha.*phi).*E(phi).*cos(thetaI).*sin(ant.gammaTest).*cos(phi)-sin(thetaI).*cos(ant.gammaTest),0,2*pi, 'ArrayValued', true);

            invZ = ant.Z^(-1);
            
            ant.xtTheta = invZ*ant.btTheta;
            
            
            ant.btTheta(1) = 0;
            ant.btTheta(end) = 0;
            
            ant.xtTheta(1) = 0;
            ant.xtTheta(end) = 0;
            
            tHatLen = sqrt(ant.tHatTest(:,1).^2+ant.tHatTest(:,3).^2);
            ftn = tHatLen.*ant.T1./ant.CoordTest(:,2)...
                +tHatLen.*ant.T2./ant.CoordTest(:,2);%T, alpha, n. Expansions function
            
            if alpha == 0
                ant.Jthe = ant.xtTheta.*ftn;
            else
                ant.Jthe = ant.Jthe+2*ant.xtTheta.*ftn.*cos(alpha.*phi);
            end
            ant.Jthe(1) = 0;
            ant.Jthe(end) = 0;
            area = emission(obj, ant, area, alpha, k, w, phiS);
        end

        function area = emission(obj, ant, area, alpha, k, w, phiS)
            rz = (area.z-ant.CoordTest(:,1));
            for i=1:length(ant.T1)
                rx = (area.x+ant.CoordTest(:,2)+area.SingularityProtection);
                r1 = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                
                B1 = -((1i.*w.*area.mu0)./(2.*pi.*r1)).*(exp(-1i.*k.*r1));
                if alpha == 0
                    area.Ethethe = area.Ethethe + B1/2 .* ant.xtTheta(i) .* ant.btTheta(i);
                else    
                    area.Ethethe = area.Ethethe + B1.*ant.xtTheta(i).*ant.btTheta(i);
                end
        
                rx = (area.x-ant.CoordTest(:,2)-area.SingularityProtection);
                
                r1 = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
                
                B1 = -((1i.*w.*area.mu0)./(2.*pi.*r1)).*(exp(-1i.*k.*r1));
                if alpha == 0
                    area.Ethethe = area.Ethethe + B1/2 .* ant.xtTheta(i) .* ant.btTheta(i);
                else
                    area.Ethethe = area.Ethethe + B1.*ant.xtTheta(i).*ant.btTheta(i);
                end
            end
        end
    end
end

