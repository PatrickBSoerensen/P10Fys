model = createpde;
gd = importGeometry(model,'ShortAntMesh.stl');

%             geometryFromEdges(model,gd);
            %0.004 er h�jeste for at f� trekanter i enderne
            generateMesh(model, 'Hmin', 0.00000001,'Hmax',5, 'GeometricOrder', 'quadratic');
            figure(1)
            pdeplot3D(model, 'FaceAlpha', 0)