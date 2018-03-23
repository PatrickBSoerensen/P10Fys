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
UV = zeros(size(stl.vertices));
faces=stl.faces;
count = 1;
for i=1:length(stl.vertices)
    a = stl.vertices==stl.vertices(i,:);
    a = sum(a,2);
    b = a==3;
    k=find(b);
    UV(k(1),:) = stl.vertices(k(1),:);
    for j=2:length(k)
        c=stl.faces == k(j);
        faces(c) = k(1);
    end
end
for i=1:length(UV)
    UV(i,4) = i;
end
k = find(~any(UV(:,1:3),2));
UV( ~any(UV(:,1:3),2), : ) = []; 
for i=1:length(UV)
    faces(faces==UV(i,4)) = i;
end
UV(:,4) = []; 

%% Gibson connectivity list
ConnectCell = cell(length(UV),2);
for i=1:length(UV)
    a = faces==i;
    b = sum(a,2);
    b = find(b);
    a(~any(a,2),:) = logical(1-a(~any(a,2),:));
    a = logical(1-a);
    
    Connected = unique(faces(a));
    Connected = Connected(Connected>i);
    ConnectCell{i,1} = Connected;
    
    ConnectCell{i,2} = unique(b);
end

EdgeNumber = 1;

for i=1:length(UV)
    temp = faces(ConnectCell{i,2},:);
    nodes = ConnectCell{i,1};
    temp = sort(temp,2);
    for j=1:length(nodes)
        a = temp(:,2:end) == nodes(j);
        b = sum(a,2);
        b = find(b);
        NodesOfInterrest = temp(b,:);
        c = NodesOfInterrest ~= i;
        d = logical(NodesOfInterrest ~= nodes(j));
        c = logical(c.*d);
        
        Triangles = temp(b,:);
        EdgePoints = [i,nodes(j)];
        NotEdgePoints = Triangles(c);
        
        L = abs(sqrt(sum((UV(EdgePoints(1),:)-UV(EdgePoints(2),:)).^2)));
        
        LforAP = abs(sqrt(sum((UV(EdgePoints(1),:)-UV(NotEdgePoints(1),:)).^2)));
        LforAM = abs(sqrt(sum((UV(EdgePoints(1),:)-UV(NotEdgePoints(2),:)).^2)));
        
        AP = (L*LforAP)/2;
        AM = (L*LforAM)/2;
        
        %Central point for a cluster
        EdgeList(EdgeNumber, 1) = i;
        %Connection point
        EdgeList(EdgeNumber, 2) = nodes(j);
        %Plus triangle point
        EdgeList(EdgeNumber, 3) = NotEdgePoints(1);
        %Minus triangle point
        EdgeList(EdgeNumber, 4) = NotEdgePoints(2);
        
        %Plus
        Basis{EdgeNumber,1} = @(r) (NotEdgePoints(1,:)-r) * L/(2*AP);
        %Minus
        Basis{EdgeNumber,2} = @(r) (NotEdgePoints(1,:)-r) * L/(2*AM);
        %Plus
        BasisDeriv(EdgeNumber,1) = -L/AP;
        %Minus
        BasisDeriv(EdgeNumber,2) = L/AM;
        
        EdgeNumber = EdgeNumber + 1;
    end
end
%% MoM?
    Z = zeros(EdgeNumber, EdgeNumber);
    b = 1:EdgeNumber;
    g = @(r) exp(1i*k*r)/(4*pi*r);
    
    for i=1:EdgeNumber
        
        Z(i,j) = integral()+integral();
        
    end