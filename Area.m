classdef Area
    %AREA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SingularityProtection;
        x;
        z;
        Ethethe;
        Ephithe;
        Ethephi; 
        Ephiphi;
        mu0;
        rx;
        rz;
        r;
    end
    
    methods
        function obj = Area(SingularityProtection, xsteps, zsteps, xmin, xmax, zmin, zmax, mu0)
            %AREA Construct an instance of this class
            %   Detailed explanation goes here
            if SingularityProtection
                obj.SingularityProtection = (2*xmax).^2/xsteps.^2;
            else
                obj.SingularityProtection = 0;
            end
            
            obj.x = linspace(xmin, xmax, xsteps);
            obj.z = linspace(zmin, zmax, zsteps);
            
            obj.Ethethe = 0;
            obj.Ephithe = 0;
            obj.Ethephi = 0; 
            obj.Ephiphi = 0;
            
            obj.mu0 = mu0;
        end
    end
end

