classdef Antenna
    properties
        Length;
        SegmentsLine;
        SegmentsCircle;
        Radii;
        Centre;
        Lin;
        LinTest;
        CircBot;
        CircBotTest;
        CircTop;
        CircTopTest;
        Coord;
    end
    
    methods
        function ant = Antenna(length, segmentsline, segmentscircle, radii, centre)
            %Setting up parameters
            ant.Length = length;
            ant.SegmentsLine = segmentsline;
            ant.SegmentsCircle = segmentscircle;
            ant.Radii = radii;
            ant.Centre = centre;
            %Splitting cylinder part of antenna
            cylinderlen = length-2.*radii;
            Lin = linspace(-cylinderlen./2, cylinderlen./2, segmentsline);
            ant.Lin = Lin(2:segmentsline-1);
            ant.LinTest = linspace(-cylinderlen./2+segmentsline./2, cylinderlen./2-segmentsline./2, segmentsline-1);
            %Splitting lower circ part
            arclength = pi/(2*segmentscircle);
            dist = linspace(-pi/2, 0, segmentscircle);
            circbot = [];
            circbot(:,1) = radii*cos(dist);
            circbot(:,2) = (-length/2)+radii+radii*sin(dist);
            ant.CircBot = circbot;
            dist = linspace(-pi/2+arclength, 0-arclength, segmentscircle-1);
            circbot = [];
            circbot(:,1) = radii*cos(dist);
            circbot(:,2) = (-length/2)+radii+radii*sin(dist);
            ant.CircBotTest = circbot;
            %Splitting upper circ part
            dist = linspace(0, pi/2, segmentscircle);
            circtop = [];
            circtop(:,1) = radii*cos(dist);
            circtop(:,2) = (length/2)-radii+radii*sin(dist);
            ant.CircTop = circtop;            
            dist = linspace(-pi/2+arclength, 0-arclength, segmentscircle-1);
            circtop = [];
            circtop(:,1) = radii*cos(dist);
            circtop(:,2) = (length/2)-radii+radii*sin(dist);
            ant.CircTopTest = circtop;
            ant.Coord = CreateCoord(ant);
        end
        
        function coord = CreateCoord(ant)
           coord = [];
           coord(1:ant.SegmentsCircle,1) = ant.CircBot(:,2);
           coord(1:ant.SegmentsCircle,2) = ant.CircBot(:,1);%Should be expanded with centre
           coord(1:ant.SegmentsCircle,3) = pi/(2*ant.SegmentsCircle);
           
           a=ant.SegmentsCircle+1;
           b=ant.SegmentsCircle+ant.SegmentsLine-1;
           
           coord(ant.SegmentsCircle+1:ant.SegmentsCircle+ant.SegmentsLine-2,1) = ant.Lin(1,:);
           coord(ant.SegmentsCircle+1:ant.SegmentsCircle+ant.SegmentsLine-2,2) = ant.Radii;%Should be expanded with centre
           coord(ant.SegmentsCircle+1:ant.SegmentsCircle+ant.SegmentsLine-2,3) = (ant.Length-2.*ant.Radii)/ant.SegmentsLine;
           
           coord(ant.SegmentsCircle+ant.SegmentsLine-1:2*ant.SegmentsCircle+ant.SegmentsLine-2,1) = ant.CircTop(:,2);
           coord(ant.SegmentsCircle+ant.SegmentsLine-1:2*ant.SegmentsCircle+ant.SegmentsLine-2,2) = ant.CircTop(:,1);%Should be expanded with centre
           coord(ant.SegmentsCircle+ant.SegmentsLine-1:2*ant.SegmentsCircle+ant.SegmentsLine-2,3) = pi/(2*ant.SegmentsCircle);%Should be expanded with centre
           
           
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