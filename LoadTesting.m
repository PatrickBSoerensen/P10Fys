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
    for n=2:length(k)
        c=stl.faces == k(n);
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
r = @(x,y,z) sqrt(x.^2+y.^2+z.^2);

for i=1:length(UV)
    temp = faces(ConnectCell{i,2},:);
    nodes = ConnectCell{i,1};
    temp = sort(temp,2);
    for n=1:length(nodes)
        a = temp(:,2:end) == nodes(n);
        b = sum(a,2);
        b = find(b);
        NodesOfInterrest = temp(b,:);
        c = NodesOfInterrest ~= i;
        d = logical(NodesOfInterrest ~= nodes(n));
        c = logical(c.*d);
        
        Triangles = temp(b,:);
        EdgePoints = [i,nodes(n)];
        NotEdgePoints = Triangles(c);
        
        L = abs(sqrt(sum((UV(EdgePoints(1),:)-UV(EdgePoints(2),:)).^2)));
        
        LforAP = abs(sqrt(sum((UV(EdgePoints(1),:)-UV(NotEdgePoints(1),:)).^2)));
        LforAM = abs(sqrt(sum((UV(EdgePoints(1),:)-UV(NotEdgePoints(2),:)).^2)));
        
        AP = (L*LforAP)/2;
        AM = (L*LforAM)/2;
        
        %Central point for a cluster
        EdgeList(EdgeNumber, 1) = i;
        %Connection point
        EdgeList(EdgeNumber, 2) = nodes(n);
        %Plus triangle point
        EdgeList(EdgeNumber, 3) = NotEdgePoints(1);
        %Minus triangle point
        EdgeList(EdgeNumber, 4) = NotEdgePoints(2);
        
        %Plus
        Basis{EdgeNumber,1} = @(r) (NotEdgePoints(1,:)-r) * L/(2*AP);
        %Minus
        Basis{EdgeNumber,2} = @(r) (r-NotEdgePoints(1,:)) * L/(2*AM);
        %Plus
        BasisDeriv(EdgeNumber,1) = -L/AP;
        %Minus
        BasisDeriv(EdgeNumber,2) = L/AM;
        
        EdgeNumber = EdgeNumber + 1;
    end
end
%% MoM?
EdgeNumber = EdgeNumber-1;
    Z = zeros(EdgeNumber, EdgeNumber);
    b = 1:EdgeNumber;
    k = 4;
    g = @(r) exp(1i.*k.*r)./(4.*pi.*r);
    
    for i=1:2
    for m=1:EdgeNumber
        mLim(1,:) = UV(EdgeList(m,1),:);
        mLim(2,:) = UV(EdgeList(m,2),:);
        mLim(3,:) = UV(EdgeList(m,3),:);
        mUpper = max(mLim);
        mLower = min(mLim);
        basem = Basis{m,i};
        basederm = BasisDeriv(m,i);
        
        for n=1:EdgeNumber
            nLim(1,:) = UV(EdgeList(n,1),:);
            nLim(2,:) = UV(EdgeList(n,2),:);
            nLim(3,:) = UV(EdgeList(n,4),:);
            nUpper = max(nLim);
            nLower = min(nLim);
            
            basen = Basis{n,i};
            basedern = BasisDeriv(n,i);
            
            a = integral3(@(x,y,z) basem(r(x,y,z)).*g(r(x,y,z)), mLower(1,1), mUpper(1,1), mLower(1,2), mUpper(1,2), mLower(1,3), mUpper(1,3));
                     
            b = integral3(@(x,y,z) basen(r(x,y,z)), nLower(1,1), nUpper(1,1), nLower(1,2), nUpper(1,2), nLower(1,3), nUpper(1,3));
                     
            c = -1/k^2.*integral3(@(x,y,z) basederm*g(r(x,y,z)), mLower(1,1), mUpper(1,1), mLower(1,2), mUpper(1,2), mLower(1,3), mUpper(1,3));
                     
            d = basedern;
        
            Z(m,n) = a.*b-c.*d + Z(m,n);
            
        end
    end
    end