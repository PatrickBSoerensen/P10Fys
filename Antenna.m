classdef Antenna
    properties
        Length;
        PointsLine;
        PointsCircle;
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
    end
    
    methods
        function ant = Antenna(length, pointsline, pointscircle, radii, centre)
            %Setting up parameters
            ant.Length = length;
            ant.PointsLine = pointsline;
            ant.PointsCircle = pointscircle;
            ant.Radii = radii;
            ant.Centre = centre;
            ant.Segments = 2*pointscircle + pointsline-2;
            %Splitting cylinder part of antenna
            cylinderlen = length-2.*radii;
            Lin = linspace(-cylinderlen./2, cylinderlen./2, pointsline);
            ant.Lin = Lin(2:pointsline-1);
            seglen = cylinderlen/pointsline;
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
        end
        
        function [tHat, zHat, gamma] = UnitVecs(ant)
        %Unit vectors should be part of Antenna class
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

        zHat = tHat; %For dimensions
        zHat(:,1) = 0; %x coord
        zHat(:,2) = 0; %y coord
        zHat(:,3) = 1; %z coord
        
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
           coord(1:segcirc,1) = circbot(:,2);
           coord(1:segcirc,2) = circbot(:,1);%Should be expanded with centre
           coord(1:segcirc,3) = pi/(2*ant.PointsCircle)*ant.Radii;
           
           coord(segcirc+1:segcirc+ant.PointsLine-2+test,1) = lin;
           coord(segcirc+1:segcirc+ant.PointsLine-2+test,2) = ant.Radii;%Should be expanded with centre
           coord(segcirc+test:segcirc+ant.PointsLine-2+test,3) = (ant.Length-2.*ant.Radii)/ant.PointsLine;
          
           coord(segcirc+ant.PointsLine-1+test:2*segcirc+ant.PointsLine-2+test,1) = circtop(:,2);
           coord(segcirc+ant.PointsLine-1+test:2*segcirc+ant.PointsLine-2+test,2) = circtop(:,1);%Should be expanded with centre
           coord(segcirc+ant.PointsLine-1+test:2*segcirc+ant.PointsLine-2+test,3) = pi/(2*ant.PointsCircle)*ant.Radii;%Should be expanded with centre
        end
        
        function green = Green(ant, k)
           coord = CreateCoord(ant);
           
           f0=@(z) (0.5-z).*(1.0-z)/(0.5-0)/(1.0-0);
           f1=@(z) (0-z).*(1.0-z)/(0-0.5)/(1.0-0.5);
           f2=@(z) (0-z).*(0.5-z)/(0-1)/(0.5-1);
           %For lin
           for i = 1:2*ant.SegmentsCircle+ant.SegmentsLine
%                     dist = sqrt((coord(i,1) - coord(:,1)).^2 + (coord(i,2) - coord(:,2)).^2 ...
%                     +2*coord(i,1)*coord(:,2)*(1-cos(phi-phimark)));
           
               if i < length(ant.CircBot)
               
%                   r = @(z) sqrt((z-ant.Lin(i)).^2+().^2);
               end
               if length(ant.CircBot) < i < length(ant.CircBot)+length(ant.Lin) 
                green = @(z) (exp(1i.*k.*r(z))./4.*pi.*r(z)).*(1+(1i./r(z).*k)-1./(r(z).*k).^2 ...
                -(ant.Lin(i)-z).^2./r(z).^2.*(1+3.*1i./(k.*r(z))-3./(k.*r(z)).^2));
               end
               if length(ant.CircBot)+length(ant.Lin) < i < length(ant.CircBot)+length(ant.Lin)+length(ant.CircTop)
               end
           end
           %For lin mirror
           %For circ
           
           %For circ mirror
           %Add integrals
           
        end
    end
end