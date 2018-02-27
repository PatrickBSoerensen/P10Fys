classdef area
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
    end
    
    methods
        function obj = area(SingularityProtection, xsteps, zsteps, xmin, xmax, zmin, zmax, mu0)
            %AREA Construct an instance of this class
            %   Detailed explanation goes here
            obj.SingularityProtection = SingularityProtection;
            obj.x = linspace(xmin, xmax, xsteps);
            obj.z = linspace(zmin, zmax, zsteps);
            
            obj.Ethethe = 0;
            obj.Ephithe = 0;
            obj.Ethephi = 0; 
            obj.Ephiphi = 0;
            
            obj.mu0 = mu0;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

