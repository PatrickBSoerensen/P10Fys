% model = createpde;
% gd = importGeometry(model,'ShortAntMesh.stl');


% stl2 = stlread('AntBinMesh2556.stl');
%             geometryFromEdges(model,gd);
            %0.004 er højeste for at få trekanter i enderne
% generateMesh(model, 'Hmin', 0.001,'Hmax', 3, 'GeometricOrder', 'quadratic');
%             figure(1)
% pdeplot3D(model, 'FaceAlpha', 0)
%Load the binary stl file
stl = stlread('AntBinMesh.stl');
%% faces and unique vertices
faces=zeros(size(stl.faces));
for i=1:length(stl.vertices)
    a = stl.vertices==stl.vertices(i,:);
    a = sum(a,2);
    b = a==3;
    k=find(b);
    for j=1:length(k)
        faces(stl.faces(:,:) == k(j)) = k(1);
    end
    if all(faces)
        break
    end
end
UV=stl.vertices(faces);
%% Gibson connectivity list
for i=1:length(UV)
    a = faces == i;
end
%% Coupling matrix
% N=length(UV);
% gmat=zeros(N,N);
% for i=1:N
%     xi=p(1,i); yi=p(2,i);
%     for j=1:i
%         xj=p(1,j); yj=p(2,j);
%         r_ij=sqrt((xi-xj)^2+(yi-yj)^2);
%         gmat(i,j)=1i/4*besselh(0,1,k0*n1*r_ij)*k0^2*(eps2-eps1);
%     end
% end
