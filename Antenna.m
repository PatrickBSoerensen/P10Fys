classdef Antenna
    properties
        Length;
        PointsLine;
        PointsCircle;
        Points;
        Segments;
        Radii;
        Centre;
        Lin;
        LinTest;
        CircBot;
        CircBotTest;
        CircTop;
        CircTopTest;
        Coord;
        CoordTest;
        zHat;
        tHat;
        gamma;
        tHatTest;
        zHatTest;
        gammaTest;
        Jthe;
        Jphi;
        Z;
        invZ;
        btTheta;
        btPhi;
        bPhiTheta;
        bPhiPhi;
        xtTheta;
        xPhiTheta;
        xtPhi;
        xPhiPhi;
        T1;
        T2;
        T1D;
        T2D;
        E0;
    end
    
    methods
        function ant = Antenna(length, pointsline, pointscircle, radii, centre, generator)
            %Setting up parameters
            ant.Length = length;
            ant.PointsLine = pointsline;
            ant.PointsCircle = pointscircle;
            ant.Radii = radii;
            ant.Centre = centre;
            ant.Points = 2*pointscircle + pointsline-2;
            ant.Segments = 2*pointscircle + pointsline-3;
            %Splitting cylinder part of antenna
            cylinderlen = length-2.*radii;
            Lin = linspace(-cylinderlen./2, cylinderlen./2, pointsline);
            ant.Lin = Lin(2:pointsline-1);
            seglen = cylinderlen/(pointsline-1);
            ant.LinTest = linspace(-cylinderlen./2+seglen./2, cylinderlen./2-seglen./2, pointsline-1);
            %Splitting lower circ part
            ant.CircBot = CreateCirc(ant, 1, 0);
            ant.CircBotTest = CreateCirc(ant, 1, 1);
            %Splitting upper circ part
            ant.CircTop = CreateCirc(ant, 0, 0);            
            ant.CircTopTest = CreateCirc(ant, 0, 1);
            %Creating full coordinate matrix
            ant.Coord = CreateCoord(ant, 0);
            ant.CoordTest = CreateCoord(ant, 1);
            [ant.tHat, ant.zHat, ant.gamma] = UnitVecs(ant);
            [ant.tHatTest, ant.zHatTest, ant.gammaTest] = UnitVecsTest(ant);
            %Testing functions triangle
            [ant.T1, ant.T2, ant.T1D, ant.T2D] = TriangleBasis(ant);
            %Current density
            ant.Jthe = (1:ant.Segments).';
            ant.Jphi = (1:ant.Segments).';
%           ant.Z = zeros(2*ant.Segments,2*ant.Segments);
            ant.Z = zeros(ant.Segments,ant.Segments);
            ant.Jthe = (1:ant.Segments-1).';
            ant.Jphi = (1:ant.Segments-1).';
            ant.Z = zeros(ant.Segments-1,ant.Segments-1);
            ant.invZ = ant.Z;
            ant.btTheta = (1:ant.Segments);
            ant.xPhiTheta = (ant.Segments+1:2*ant.Segments);    
            ant.T1 = sqrt((ant.CoordTest(:,1)-ant.Coord(1:ant.Segments-1,1)).^2 ... 
            ant.btTheta = (1:ant.Segments-1);
            ant.btPhi = (1:ant.Segments-1);
            ant.bPhiTheta = (1:ant.Segments-1);
            ant.bPhiPhi = (1:ant.Segments-1);
            ant.xtTheta = (1:ant.Segments-1);
            ant.xPhiTheta = (ant.Segments:2*ant.Segments-1);    
            ant.xtPhi = (1:ant.Segments-1);
            ant.xPhiPhi = (ant.Segments:2*ant.Segments-1);
            %Field limiter
            if generator
                ant.E0 = FieldSetup(ant, ant.Length/10);
            else
                ant.E0(1:ant.Segments) = 0;
            end
            
        end
        
        function [T1, T2, T1D, T2D] = TriangleBasis(ant)
%           Rising basis function
            T1 = sqrt((ant.CoordTest(:,1)-ant.Coord(1:ant.Segments,1)).^2 ... 
                +(ant.CoordTest(:,2)-ant.Coord(1:ant.Segments,2)).^2)...
                ./sqrt((ant.Coord(1:ant.Segments,1)-ant.Coord(2:ant.Segments+1,1)).^2 ...
                +(ant.Coord(1:ant.Segments,2)-ant.Coord(2:ant.Segments+1,2)).^2);
            T1 = T1(1:end-1);
%           Falling basis function
            T2 = sqrt((ant.Coord(2:ant.Segments+1,1)-ant.CoordTest(:,1)).^2 ... 
                +(ant.Coord(2:ant.Segments+1,2)-ant.CoordTest(:,2)).^2) ...
                ./sqrt((ant.Coord(1:ant.Segments,1)-ant.Coord(2:ant.Segments+1,1)).^2 ...
                +(ant.Coord(1:ant.Segments,2)-ant.Coord(2:ant.Segments+1,2)).^2);
            T2 = T2(2:end);
            %Rising basis derivative
            T1D = 1./sqrt((ant.Coord(1:ant.Segments,1)-ant.Coord(2:ant.Segments+1,1)).^2 ...
                +(ant.Coord(1:ant.Segments,2)-ant.Coord(2:ant.Segments+1,2)).^2);
            T1D = T1D(1:end-1);
            %Falling basis derivative
            T2D = -1./sqrt((ant.Coord(1:ant.Segments,1)-ant.Coord(2:ant.Segments+1,1)).^2 ...
                +(ant.Coord(1:ant.Segments,2)-ant.Coord(2:ant.Segments+1,2)).^2);
            T2D = T2D(2:end);
            
        end
        
        function [tHat, zHat, gamma] = UnitVecs(ant)
            %tHat proper       
            tHat(1:ant.PointsCircle,1) ... 
            = -ant.Radii.*sin(linspace(-pi/2, 0, ant.PointsCircle));
            tHat(ant.PointsCircle+1:ant.PointsCircle+ant.PointsLine-2,1) ... 
            = 0;
            tHat(ant.PointsCircle+ant.PointsLine-1: ... 
            2*ant.PointsCircle+ant.PointsLine-2,1)...
            = -ant.Radii.*sin(linspace(0, pi/2, ant.PointsCircle));%x coord
        
            tHat(:,2) = 0;%y coord
            tHat(:,3) = 1;%z coord
            tHat(1:ant.PointsCircle,3) = ... 
            ant.Radii.*cos(linspace(-pi/2, 0, ant.PointsCircle));%z coord
            tHat(ant.PointsCircle+ant.PointsLine-1: ... 
            2*ant.PointsCircle+ant.PointsLine-2,3)...
            = ant.Radii.*cos(linspace(0, pi/2, ant.PointsCircle));%z coord
    
            tHat = tHat./sqrt(tHat(:,1).^2+tHat(:,2).^2+tHat(:,3).^2);%normalizing to unit
            %zHat proper
            zHat = tHat; %For dimensions
            zHat(:,1) = 0; %x coord
            zHat(:,2) = 0; %y coord
            zHat(:,3) = 1; %z coord
            %gamma proper
            gamma = acos(dot(tHat,zHat,2));    
        end
        
        function [tHat, zHat, gamma] = UnitVecsTest(ant)
            arclength = pi/(2*ant.PointsCircle);
            distbot = linspace(-pi/2+arclength/2, 0-arclength/2, ant.PointsCircle-1);
            disttop = linspace(0+arclength/2, pi/2-arclength/2, ant.PointsCircle-1);
            
            %tHat test       
            tHat(1:ant.PointsCircle-1,1) ... 
            = -ant.Radii.*sin(distbot);
            tHat(ant.PointsCircle:ant.PointsCircle+ant.PointsLine-2,1) ... 
            = 0;
            tHat(ant.PointsCircle+ant.PointsLine-1: ... 
            2*ant.PointsCircle+ant.PointsLine-3,1)...
            = -ant.Radii.*sin(disttop);%x coord
        
            tHat(:,2) = 0;%y coord
            tHat(:,3) = 1;%z coord
            tHat(1:ant.PointsCircle-1,3) = ... 
            ant.Radii.*cos(distbot);%z coord
            tHat(ant.PointsCircle+ant.PointsLine-1: ... 
            2*ant.PointsCircle+ant.PointsLine-3,3)...
            = ant.Radii.*cos(disttop);%z coord
    
            tHat = tHat./sqrt(tHat(:,1).^2+tHat(:,2).^2+tHat(:,3).^2);%normalizing to unit
            %zHat test
            zHat = tHat; %For dimensions
            zHat(:,1) = 0; %x coord
            zHat(:,2) = 0; %y coord
            zHat(:,3) = 1; %z coord
            %gamma test
            gamma = acos(dot(tHat,zHat,2));    
        end
        
        function circ = CreateCirc(ant, bot, test)
            arclength = pi/(2*ant.PointsCircle);
            circ = [];
            if test
                if bot
                    dist = linspace(-pi/2+arclength/2, 0-arclength/2, ant.PointsCircle-1);
                    circ(:,2) = (-ant.Length/2)+ant.Radii+ant.Radii*sin(dist);
                else
                    dist = linspace(0+arclength/2, pi/2-arclength/2, ant.PointsCircle-1);
                    circ(:,2) = (ant.Length/2)-ant.Radii+ant.Radii*sin(dist);
                end
                circ(:,1) = ant.Radii*cos(dist);
            else
                if bot
                    dist = linspace(-pi/2, 0, ant.PointsCircle);
                    circ(:,2) = (-ant.Length/2)+ant.Radii+ant.Radii*sin(dist);
                else
                    dist = linspace(0, pi/2, ant.PointsCircle);
                    circ(:,2) = (ant.Length/2)-ant.Radii+ant.Radii*sin(dist);
                end
                circ(:,1) = ant.Radii*cos(dist);
            end
        end
        
        function coord = CreateCoord(ant, test)
           coord = [];
           if test
               segcirc = ant.PointsCircle-1;
               circbot = ant.CircBotTest;
               circtop = ant.CircTopTest;
               lin = ant.LinTest;
           else
               segcirc = ant.PointsCircle;
               circbot = ant.CircBot;
               circtop = ant.CircTop;
               lin = ant.Lin;
           end
           coord(1:segcirc,1) = circbot(:,2)+ant.Centre(1);
           coord(1:segcirc,2) = circbot(:,1)+ant.Centre(2);
           coord(1:segcirc,3) = ant.Radii*pi/(2*ant.PointsCircle);
           
           coord(segcirc+1:segcirc+ant.PointsLine-2+test,1) = lin+ant.Centre(1);
           coord(segcirc+1:segcirc+ant.PointsLine-2+test,2) = ant.Radii+ant.Centre(2);
           coord(segcirc+test:segcirc+ant.PointsLine-2+test,3) = (ant.Length-2.*ant.Radii)/(ant.PointsLine-1);
          
           coord(segcirc+ant.PointsLine-1+test:2*segcirc+ant.PointsLine-2+test,1) = circtop(:,2)+ant.Centre(1);
           coord(segcirc+ant.PointsLine-1+test:2*segcirc+ant.PointsLine-2+test,2) = circtop(:,1)+ant.Centre(2);
           coord(segcirc+ant.PointsLine-1+test:2*segcirc+ant.PointsLine-2+test,3) = ant.Radii*pi/(2*ant.PointsCircle);
        end
        
        function E0 = FieldSetup(ant, lim)
            E0 = (1:ant.Segments);
            E0(:) = 0;
            
            upper1 = 0 < ant.CoordTest(:, 1);
            upper2 = ant.CoordTest(:, 1) <= lim/2;
            upper = logical(upper1.*upper2);
            
            mid = ant.CoordTest(:,1) == 0;
            
            lower1 = 0 > ant.CoordTest(:, 1);
            lower2 = ant.CoordTest(:, 1) >= -lim/2;
            lower = logical(lower1.*lower2);
            
            amount = round(ant.Segments/20);
            
            if sum(mid)
                E0(mid) = 1;
            end
            
            E0(upper) = (ant.CoordTest(round(ant.Segments/2+amount),1)-ant.CoordTest(upper, 1))...
            ./(ant.CoordTest(round(ant.Points/2+amount),1));
            E0(lower) = (ant.CoordTest(round(ant.Points/2-amount),1)-ant.CoordTest(lower, 1))...
            ./(ant.CoordTest(round(ant.Points/2-amount),1));
        
        end
    end
end