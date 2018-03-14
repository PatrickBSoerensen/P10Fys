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
%              [obj.G1, obj.G2] = DistanceSetup(obj, ant);
        end
        
        function [G1, G2] = DistanceSetup(obj, ant)
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
              
                    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    
                    %Should possibly integrate to 2*pi
                    G1{i,h} = ant.CoordTest(iSegments(i),3).*ant.CoordTest(jSegments(h),3).*integral(Func1, 0, pi);%, 'ArrayValued', true);
                    G2{i,h} = ant.CoordTest(iSegments(i),3).*ant.CoordTest(jSegments(h),3).*integral(Func2, 0, pi);%, 'ArrayValued', true);
                    end
                end
                
            obj.G1 = triu(cell2mat(G1));
            obj.G2 = triu(cell2mat(G2));
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
              
                    Func1 = @(phimark) cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    Func2 = @(phimark) cos(phimark).*cos(alpha.*phimark).*exp(-1i.*k.*R(phimark))./R(phimark);
                    
                    %Should possibly integrate to 2*pi
                    G1{i,h} = ant.CoordTest(iSegments(i),3).*ant.CoordTest(jSegments(h),3).*integral(Func1, 0, pi);%, 'ArrayValued', true);
                    G2{i,h} = ant.CoordTest(iSegments(i),3).*ant.CoordTest(jSegments(h),3).*integral(Func2, 0, pi);%, 'ArrayValued', true);
                    end
                end
                
            obj.G1 = triu(cell2mat(G1));
            obj.G2 = triu(cell2mat(G2));
            clear G1; clear G2; clear Func1; clear Func2;
            end
            
%             %Ztt
%             ant.Z = ...
%                 ant.T1(:)*ant.T1(:).'.*...
%                 (sin(ant.gammaTest(1:end-1))*sin(ant.gammaTest(1:end-1).').*obj.G2(:,:) ...%1
%                 +cos(ant.gammaTest(1:end-1))*cos(ant.gammaTest(1:end-1).').*obj.G1(:,:)) ...
%                 -1./k.^2.*ant.T1D(:)*ant.T1D(:).'.*obj.G1(:,:)...
%                 +...
%                 ant.T1(:)*ant.T2(:).'.*...
%                 (sin(ant.gammaTest(1:end-1))*sin(ant.gammaTest(2:end).').*obj.G2(:,:) ...%3
%                 +cos(ant.gammaTest(1:end-1))*cos(ant.gammaTest(2:end).').*obj.G1(:,:)) ...
%                 -1./k.^2.*ant.T1D(:)*ant.T2D(:).'.*obj.G1(:,:)...
%                 +...
%                 ant.T2(:)*ant.T2(:).'.*...
%                 (sin(ant.gammaTest(2:end))*sin(ant.gammaTest(2:end).').*obj.G2(:,:) ...%2
%                 +cos(ant.gammaTest(2:end))*cos(ant.gammaTest(2:end).').*obj.G1(:,:)) ...
%                 -1./k.^2.*ant.T2D(:)*ant.T2D(:).'.*obj.G1(:,:)...
%                 +...
%                 ant.T2(:)*ant.T1(:).'.*...
%                 (sin(ant.gammaTest(2:end))*sin(ant.gammaTest(1:end-1).').*obj.G2(:,:) ...%4
%                 +cos(ant.gammaTest(2:end))*cos(ant.gammaTest(1:end-1).').*obj.G1(:,:))...
%                 -1./k.^2.*ant.T2D(:)*ant.T1D(:).'.*obj.G1(:,:);

            ant.Z = ...
                ant.T1(:)*ant.T1(:).'.*...
                (sin(ant.gamma(2:end-1))*sin(ant.gamma(2:end-1).').*obj.G2(:,:) ...%1
                +cos(ant.gamma(2:end-1))*cos(ant.gamma(2:end-1).').*obj.G1(:,:)) ...
                -1./k.^2.*ant.T1D(:)*ant.T1D(:).'.*obj.G1(:,:)...
                +...
                ant.T1(:)*ant.T2(:).'.*...
                (sin(ant.gamma(2:end-1))*sin(ant.gamma(2:end-1).').*obj.G2(:,:) ...%1
                +cos(ant.gamma(2:end-1))*cos(ant.gamma(2:end-1).').*obj.G1(:,:)) ...
                -1./k.^2.*ant.T1D(:)*ant.T2D(:).'.*obj.G1(:,:)...
                +...
                ant.T2(:)*ant.T2(:).'.*...
                (sin(ant.gamma(2:end-1))*sin(ant.gamma(2:end-1).').*obj.G2(:,:) ...%1
                +cos(ant.gamma(2:end-1))*cos(ant.gamma(2:end-1).').*obj.G1(:,:)) ...
                -1./k.^2.*ant.T2D(:)*ant.T2D(:).'.*obj.G1(:,:)...
                +...
                ant.T2(:)*ant.T1(:).'.*...
                (sin(ant.gamma(2:end-1))*sin(ant.gamma(2:end-1).').*obj.G2(:,:) ...%1
                +cos(ant.gamma(2:end-1))*cos(ant.gamma(2:end-1).').*obj.G1(:,:)) ...
                -1./k.^2.*ant.T2D(:)*ant.T1D(:).'.*obj.G1(:,:);   
            
            ant.Z = ant.Z+triu(ant.Z,1).';
            
            J0 = besselj(alpha-1, k*ant.Coord(2:end-1,2)*sin(thetaI));
            J1 = besselj(alpha, k*ant.Coord(2:end-1,2)*sin(thetaI));
            J2 = besselj(alpha+1, k*ant.Coord(2:end-1,2)*sin(thetaI));
                
            %% planewave b equations
%             ant.btTheta = -ant.E0(1:end-1).*1i./(w.*mu).*pi.*1i.^(alpha).*...
%                 (ant.T1.*ant.CoordTest(1:end-1,3)...
%                 .*exp(1i.*k.*ant.CoordTest(1:end-1,1).*cos(thetaI)).*(cos(thetaI)...
%                 .*sin(ant.gammaTest(1:end-1)).*1i.*(J2(1:end-1)-J0(1:end-1))...
%                 -2.*sin(thetaI).*cos(ant.gammaTest(1:end-1)).*J1(1:end-1)))...
%                 -ant.E0(2:end).*1i./(w.*mu).*pi.*1i.^(alpha).*...
%                 (ant.T2.*ant.CoordTest(2:end,3)...
%                 .*exp(1i.*k.*ant.CoordTest(2:end,1).*cos(thetaI)).*(cos(thetaI)...
%                 .*sin(ant.gammaTest(2:end)).*1i.*(J2(2:end)-J0(2:end))...
%                 -2.*sin(thetaI).*cos(ant.gammaTest(2:end)).*J1(2:end)));

            ant.btTheta = -ant.E0(2:end-1).*1i./(w.*mu).*pi.*1i.^(alpha).*...
                (ant.T1.*ant.T2.*ant.Coord(2:end-1,3)...
                .*exp(1i.*k.*ant.Coord(2:end-1,1).*cos(thetaI)).*(cos(thetaI)...
                .*sin(ant.gamma(2:end-1)).*1i.*(J2-J0)...
                -2.*sin(thetaI).*cos(ant.gamma(2:end-1)).*J1));
            
            ant.btTheta(1)=0;
            ant.btTheta(end)=0;
            
            invZ = ant.Z^(-1);
            
            ant.xtTheta =invZ*ant.btTheta;
            
            ant.xtTheta(1)=0;
            ant.xtTheta(end)=0;
            
            tHatLen = sqrt(ant.tHat(2:end-1,1).^2+ant.tHat(2:end-1,3).^2);
            ftn = tHatLen.*ant.T1.*ant.T2./ant.Coord(2:end-1,2);%...
%                 +tHatLen(2:end).*ant.T2./ant.CoordTest(2:end,2);%T, alpha, n. Expansions function
            
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
            rz = (area.z-ant.Coord(2:end-1,1)-area.SingularityProtection);
            for i=1:length(ant.T1)
                rx = (area.x+ant.Coord(2:end-1,2)+area.SingularityProtection);
                r1 = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
%                 r2 = sqrt((rz(i+1,:).').^2+(rx(i+1,:)).^2);
                
                B1 = -((1i.*w.*area.mu0)./(2.*pi.*r1)).*(exp(-1i.*k.*r1));
%                 B2 = -((1i.*w.*area.mu0)./(2.*pi.*r2)).*(exp(-1i.*k.*r2));
                if alpha == 0
                    area.Ethethe = B1/2 .* ant.xtTheta(i) .* ant.btTheta(i);% + area.Ethethe;%.* ant.CoordTest(i,3)...
%                                    +B2/2 .* ant.xtTheta(i) .* ant.btTheta(i) .* ant.CoordTest(i+1,3)...
%                                    + area.Ethethe;
                else    
                    area.Ethethe = B1.*ant.xtTheta(i).*ant.btTheta(i)+ area.Ethethe;% .* ant.CoordTest(i,3)...
%                         +B2.*ant.xtTheta(i).*ant.btTheta(i) .* ant.CoordTest(i+1,3)...
%                         +area.Ethethe;
                end
        
                rx = (area.x-ant.Coord(2:end-1,2)-area.SingularityProtection);
                
                r1 = sqrt((rz(i,:).').^2+(rx(i,:)).^2);
%                 r2 = sqrt((rz(i+1,:).').^2+(rx(i+1,:)).^2);
                
                B1 = -((1i.*w.*area.mu0)./(2.*pi.*r1)).*(exp(-1i.*k.*r1));
%                 B2 = -((1i.*w.*area.mu0)./(2.*pi.*r2)).*(exp(-1i.*k.*r2));
                if alpha == 0
                    area.Ethethe = B1/2 .* ant.xtTheta(i) .* ant.btTheta(i);% + area.Ethethe;%.* ant.CoordTest(i,3) ...
%                                    +B2/2 .* ant.xtTheta(i) .* ant.btTheta(i) .* ant.CoordTest(i+1,3)...
%                                    + area.Ethethe;
                else
                    area.Ethethe = B1.*ant.xtTheta(i).*ant.btTheta(i) + area.Ethethe;%.* ant.CoordTest(i,3)...
%                         +B2.*ant.xtTheta(i).*ant.btTheta(i) .* ant.CoordTest(i+1,3)...
%                         +area.Ethethe;
                end
            end
        end
    end
end

