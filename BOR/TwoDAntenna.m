classdef TwoDAntenna < Antenna
    %TWODANTENNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gdc
        ns
        sf
        dl
        bt
        Solver
        
    end
    
    methods
        function obj = TwoDAntenna(length, pointsline, pointscircle, radii, centre, generator)
            %TWODANTENNA Construct an instance of this class
            %   Detailed explanation goes here
            obj@Antenna(length, pointsline, pointscircle, radii, centre, generator);
            [obj.gdc, obj.ns, obj.sf, obj.dl, obj.bt] = CreateGeometry(obj);
            obj.Solver = TwoDTriangulation(obj);
        end
        
        function [gdc, ns, sf, dl, bt] = CreateGeometry(ant)
            rect = [3; 4; -ant.Length/2+ant.Radii; ant.Length/2-ant.Radii; ...
                ant.Length/2-ant.Radii; -ant.Length/2+ant.Radii; ...
                -ant.Radii; -ant.Radii; ant.Radii; ant.Radii];
            UpperCirc = [1; ant.Length/2-ant.Radii; ant.Centre(2); ant.Radii];
            UpperCirc = [UpperCirc;zeros(length(rect) - length(UpperCirc),1)];
            LowerCirc = [1;  -ant.Length/2+ant.Radii; ant.Centre(2); ant.Radii];
            LowerCirc = [LowerCirc;zeros(length(rect) - length(LowerCirc),1)];
            gdc = [rect, UpperCirc, LowerCirc];
            ns = char('R', 'UC', 'LC');
            ns = ns';
%           sf = '((UC-R)+(LC-R))'; Only end semi circles
            sf = 'R+UC+LC';
            [dl, bt] = decsg(gdc,sf, ns);
            [dl, bt] = csgdel(dl,bt);
        end
        
        function solver = TwoDTriangulation(ant)
            solver = createpde;
            geometryFromEdges(solver,ant.dl);
            %0.004 er højeste for at få trekanter i enderne
            generateMesh(solver, 'Hmin', 0.000001, 'GeometricOrder', 'quadratic');
        end
    end
end

