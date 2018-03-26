%% konstants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
%%
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
    d=find(b);
    UV(d(1),:) = stl.vertices(d(1),:);
    for n=2:length(d)
        c=stl.faces == d(n);
        faces(c) = d(1);
    end
end
for i=1:length(UV)
    UV(i,4) = i;
end
d = find(~any(UV(:,1:3),2));
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
r = @(v1,v2,v3, alpha, beta) (1-alpha-beta)*v1 + alpha*v2 + beta*v3;

    for i=1:length(faces)
        C(i,1) = sum(UV(faces(i,:),1))/3;
        C(i,2) = sum(UV(faces(i,:),2))/3;
        C(i,3) = sum(UV(faces(i,:),3))/3;
        
    
    L1 = UV(faces(i,1),:)-C(i,:);
    L2 = UV(faces(i,3),:)-UV(faces(i,1),:);
        A(i,1) = sqrt(sum((cross(L1,L2)).^2,2))/2;
        
        
    L1 = UV(faces(i,1),:)-C(i,:);
    L2 = UV(faces(i,2),:)-UV(faces(i,1),:);
        A(i,2) = sqrt(sum((cross(L1,L2)).^2,2))/2;
        
    L1 = UV(faces(i,2),:)-C(i,:);
    L2 = UV(faces(i,3),:)-UV(faces(i,2),:);
        A(i,3) = sqrt(sum((cross(L1,L2)).^2,2))/2;
    end
    
    L1 = UV(faces(:,2),:)-UV(faces(:,1),:);
    L2 = UV(faces(:,3),:)-UV(faces(:,1),:);
        
    Atot = sqrt(sum((cross(L1,L2)).^2,2))/2;
    A = A./Atot;
    
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
        
        L = UV(EdgePoints(1),:)-UV(EdgePoints(2),:);
        
        LforAP = UV(EdgePoints(1),:)-UV(NotEdgePoints(1),:);
        LforAM = UV(EdgePoints(1),:)-UV(NotEdgePoints(2),:);
        
        AP =  sqrt(sum((cross(L,LforAP)).^2,2))./2;
        AM =  sqrt(sum((cross(L,LforAM)).^2,2))./2;
        
        L = sqrt(sum(L.^2));
        
        %Central point for a cluster
        EdgeList(EdgeNumber, 1) = i;
        %Connection point
        EdgeList(EdgeNumber, 2) = nodes(n);
        %Plus triangle point
        EdgeList(EdgeNumber, 3) = NotEdgePoints(1);
        %Minus triangle point
        EdgeList(EdgeNumber, 4) = NotEdgePoints(2);
        
        %Plus
        Basis{EdgeNumber,1} = @(r) (UV(NotEdgePoints(1),:) - r);
        BasisLA(EdgeNumber,1) =  L./(2*AP);
        %Minus
        Basis{EdgeNumber,2} = @(r) (r - UV(NotEdgePoints(2),:));
        BasisLA(EdgeNumber,2) = L./(2*AM);
        %Plus
        BasisDeriv(EdgeNumber,1) = -L./AP;
        %Minus
        BasisDeriv(EdgeNumber,2) = L./AM;
        
        EdgeNumber = EdgeNumber + 1;
    end
end
%% MoM?
    loader = 0;
    Z = BasisLA(:,1)*BasisLA(:,2).';
    b = 1:EdgeNumber-1;
    g = @(r) exp(1i.*k.*r)./(4.*pi.*r);
    
if ~loader
    for i=1:2
        for m=1:length(EdgeList)
            mv1 = UV(EdgeList(m,1),:);
            mv2 = UV(EdgeList(m,2),:);
            mv3 = UV(EdgeList(m,3),:);
            
            basem = Basis{m,i};
            basederm = BasisDeriv(m,i);
            
            test1 = @(alpha, beta) basem(r(mv1,mv2,mv3, alpha, beta)).*g(r(mv1,mv2,mv3, alpha, beta));
            betamax = @(alpha) alpha - 1;
            
            z1 = integral2(test1, 0, 1, 0, betamax);
            z3 = 1/k^2.*integral2(basederm*g(r(mv1,mv2,mv3, alpha, beta)), mLower(1,1), 0, 1, 0, betamax);
            for n=1:length(EdgeList)
                nv1 = UV(EdgeList(n,1),:);
                nv2 = UV(EdgeList(n,2),:);
                nv3 = UV(EdgeList(n,4),:);
                
                basen = Basis{n,i};
                basedern = BasisDeriv(n,i);
                     
                z2 = integral2(@(alpha, beta) basen(r(nv1,nv2,nv3, alpha, beta)), nLower(1,1), 0, 1, 0, betamax);
                    
                z4 = basedern;
        
                Z(m,n) = z1.*z2-z3.*z4 + Z(m,n);
            
            end
            b(m) = -1i/(w*mu0)*integral3(@(x,y,z) basem(r(x,y,z)), mLower(1,1), mUpper(1,1), mLower(1,2), mUpper(1,2), mLower(1,3), mUpper(1,3));
        end
    end
else
    load('IWaited8HoursForThis');
end
            
%% MoM time Pimped
for f=1:length(faces)
    for basis=1:2
        
    end
end
%% current  
        
%             for i=1:2
%                 for m=1:EdgeNumber
%                     
%         mLim(1,:) = UV(EdgeList(m,1),:);
%         mLim(2,:) = UV(EdgeList(m,2),:);
%         mLim(3,:) = UV(EdgeList(m,3),:);
%         mUpper = max(mLim);
%         mLower = min(mLim);
%                     basem = Basis{m,i};
%                     for n = 1:EdgeNumber
%        
%             nLim(1,:) = UV(EdgeList(n,1),:);
%             nLim(2,:) = UV(EdgeList(n,2),:);
%             nLim(3,:) = UV(EdgeList(n,4),:);
%             nUpper = max(nLim);
%             nLower = min(nLim);
%             basen = Basis{n,i};
%             ftn = integral3(@(x,y,z) basem(r(x,y,z)), mLower(1,1), mUpper(1,1), mLower(1,2), mUpper(1,2), mLower(1,3), mUpper(1,3))./CoordTest(:,2)...
%                 +integral3(@(x,y,z) basen(r(x,y,z)), mLower(1,1), mUpper(1,1), mLower(1,2), mUpper(1,2), mLower(1,3), mUpper(1,3))./CoordTest(:,2);
%             
%             if alpha == 0
%                 Jthe = xtTheta.*ftn;
%             else
%                 Jthe = Jthe+2*ant.xtTheta.*ftn.*cos(alpha.*phi);
%             end
%             ant.Jthe(1) = 0;
%             ant.Jthe(end) = 0;           
%                     end
%                 end
%             end
    %% Emission
%     r = linspace(0,10,100);
%     E = zeros(100,100);
%              for i=1:length(r)
%                 
%                 g = exp(-1i.*k.*r)./(4.*pi.*r);
%                 G = g.*((1+1i./(k*r)-1./((k*r).^2)) - ...
%                             ((rx(i,:)).^2)./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
%                 E(i,:) = E(i,:)+G.*ant.Jthe(i);
%              end