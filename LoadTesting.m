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

BasisNumberOuter = 1;

C=zeros(size(faces));
A=zeros(size(faces));
    for i=1:length(faces)
        C(i,1) = sum(UV(faces(i,:),1))/3;
        C(i,2) = sum(UV(faces(i,:),2))/3;
        C(i,3) = sum(UV(faces(i,:),3))/3;
        
    
    L1 = UV(faces(i,1),:)-C(i,:);
    L2 = UV(faces(i,3),:)-UV(faces(i,1),:);
        A(i,1) = sqrt(sum((cross(L1,L2)).^2,2))/2;
        
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
        EdgeList(BasisNumberOuter, 1) = i;
        %Connection point
        EdgeList(BasisNumberOuter, 2) = nodes(n);
        %Plus triangle point
        EdgeList(BasisNumberOuter, 3) = NotEdgePoints(1);
        %Minus triangle point
        EdgeList(BasisNumberOuter, 4) = NotEdgePoints(2);
        
        %Plus
        Basis{BasisNumberOuter,1} = @(r) (UV(NotEdgePoints(1),:) - r);
        BasisLA(BasisNumberOuter,1) =  L./(2*AP);
        BasisLA(BasisNumberOuter,2) =  L;
        %Minus
        Basis{BasisNumberOuter,2} = @(r) (r - UV(NotEdgePoints(2),:));
        BasisLA(BasisNumberOuter,3) = L./(2*AM);
        BasisLA(BasisNumberOuter,4) = L;
        %Plus
        BasisDeriv(BasisNumberOuter,1) = -L./AP;
        %Minus
        BasisDeriv(BasisNumberOuter,2) = L./AM;
        
        BasisNumberOuter = BasisNumberOuter + 1;
    end
end
    %% MoM faces loop
    BasisNumberOuter = BasisNumberOuter-1;
    loader = 0;
    Z = zeros(BasisNumberOuter,BasisNumberOuter);
    b = 1:BasisNumberOuter-1;
    g = @(p,q) exp(1i.*k.*sqrt(sum((p-q).^2)))./(4.*pi.*sqrt(sum((p-q).^2)));
    wp = 1/3;
    Ei=1;
    
for f=1:length(faces)
    OuterBasis(1,:) = faces(f,1:2);
    OuterBasis(2,:) = faces(f,2:3);
    OuterBasis(3,1) = faces(f,1);
    OuterBasis(3,2) = faces(f,3);
    
    for h=1:length(faces)
        InnerBasis(1,:) = faces(h,1:2);
        InnerBasis(2,:) = faces(h,2:3);
        InnerBasis(3,1) = faces(h,1);
        InnerBasis(3,2) = faces(h,3);
        
        for i=1:3
            BasisNumberOuter = sum(EdgeList(:,1:2)==OuterBasis(i,:),2);
            BasisNumberOuter = find(BasisNumberOuter==2);
            
            BasisNumberInner = sum(EdgeList(:,1:2)==InnerBasis(i,:),2);
            BasisNumberInner = find(BasisNumberInner==2);
            v1 = UV(faces(f,1),:);
            v2 = UV(faces(f,2),:);
            v3 = UV(faces(f,3),:);
            
            u1 = UV(faces(h,1),:);
            u2 = UV(faces(h,2),:);
            u3 = UV(faces(h,3),:);
            
            [I111, I112, I11] = NearTriangleZ(v1, v2, v3);
            
            Z(BasisNumberOuter, BasisNumberInner) = I111+I112+I11 + Z(BasisNumberOuter, BasisNumberInner);
            
            [I111, I112, I11] = NearTriangleZ(u1, u2, u3);
            
            Z(BasisNumberOuter, BasisNumberInner) = I111+I112+I11 + Z(BasisNumberOuter, BasisNumberInner);
            
            if BasisNumberOuter
                edges(f,i) = BasisNumberOuter;
            end
        end
    end
end

%% MoM edge loop
if loader
    load('IWaited8HoursForThis');
else
        for m=1:length(EdgeList)
            mv1 = UV(EdgeList(m,1),:);
            mv2 = UV(EdgeList(m,2),:);
            mv3 = UV(EdgeList(m,3),:);
            mv4 = UV(EdgeList(m,4),:);
            mv(1,:) = mv1;
            mv(2,:) = mv2;
            mv(3,:) = mv3;
            mv(4,:) = mv4;
            BaseMP = Basis{m,1};
            BaseMM = Basis{m,2};
            BaseDMP = BasisDeriv(m,1);
            BaseDMM = BasisDeriv(m,2);
            
            z1M = @(q) wp*(1/4*BaseMP(mv1)+1/k^2)*g(mv1,q);
            z2M = @(q) wp*(1/4*BaseMP(mv2)+1/k^2)*g(mv2,q);
            z3M = @(q) wp*(1/4*BaseMP(mv3)+1/k^2)*g(mv3,q);
            z4M = @(q) wp*(1/4*BaseMM(mv4)-1/k^2)*g(mv4,q);
            
            for n=1:length(EdgeList)
                nv1 = UV(EdgeList(n,1),:);
                nv2 = UV(EdgeList(n,2),:);
                nv3 = UV(EdgeList(n,3),:);
                nv4 = UV(EdgeList(n,4),:);
                nv(1,:) = nv1;
                nv(2,:) = nv2;
                nv(3,:) = nv3;
                nv(4,:) = nv4;
                BaseNP = Basis{n,1};
                BaseNM = Basis{n,2};
                BaseDNP = BasisDeriv(n,1);
                BaseDNM = BasisDeriv(n,2);
                
                for sing=1:9
                    Near(sing) = sum(sum(nv == mv(sing)));
                    NearOverlap{sing} = nv == mv(sing);
                end
                    Nearsum = sum(sum(Near));
                
               if Nearsum>=4
                   [I111, I112, I11] = NearTriangleZ(mv1, nv2, nv3);
                   
                   Z(m,n) = I111.*I112.*I11 + Z(m,n);
                   
                   [I111, I112, I11] = NearTriangleZ(nv1, mv2, mv3);
                   
                   Z(m,n) = I111.*I112.*I11 + Z(m,n); 
               else
                   
                z1N = @(q) wp*(1/4*BaseNP(nv1)+1/k^2)*g(nv1,q);
                z2N = @(q) wp*(1/4*BaseNP(nv2)+1/k^2)*g(nv2,q);
                z3N = @(q) wp*(1/4*BaseNP(nv3)+1/k^2)*g(nv3,q);
                z4N = @(q) wp*(1/4*BaseNM(nv4)-1/k^2)*g(nv4,q);
                
                for ns = 1:4
                    for ms = 1:4
                    Z(m,n) = (BasisLA(m,2)*BasisLA(n,4))/(4*pi)...
                    * dot(z1M(nv(ns,:)),z1N(mv(ms,:)))...
                    + dot(z1M(nv(ns,:)),z2N(mv(ms,:)))...
                    + dot(z1M(nv(ns,:)),z3N(mv(ms,:)))...
                    + dot(z1M(nv(ns,:)),z4N(mv(ms,:)))...
                    + dot(z2M(nv(ns,:)),z1N(mv(ms,:)))...
                    + dot(z2M(nv(ns,:)),z2N(mv(ms,:)))...
                    + dot(z2M(nv(ns,:)),z3N(mv(ms,:)))...
                    + dot(z2M(nv(ns,:)),z4N(mv(ms,:)))...
                    + dot(z3M(nv(ns,:)),z1N(mv(ms,:)))...
                    + dot(z3M(nv(ns,:)),z2N(mv(ms,:)))...
                    + dot(z3M(nv(ns,:)),z3N(mv(ms,:)))...
                    + dot(z3M(nv(ns,:)),z4N(mv(ms,:)))...
                    + dot(z4M(nv(ns,:)),z1N(mv(ms,:)))...
                    + dot(z4M(nv(ns,:)),z2N(mv(ms,:)))...
                    + dot(z4M(nv(ns,:)),z3N(mv(ms,:)))...
                    + dot(z4M(nv(ns,:)),z4N(mv(ms,:)))...
                    + Z(m,n);
                    end
                end
                end
            end
             b(m) = sum(BasisLA(m,2)/2*wp*(BaseMP(mv1)+BaseMP(mv2)+BaseMP(mv3))*Ei);
        end
end
             zinv = inv(Z);
            xtesst = b*zinv;
            
%             Z2 = zeros(EdgeNumber,EdgeNumber);
%% Alternate Edge loop
for m=1:length(EdgeList)
    
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

function [I111, I112, I11] = NearTriangleZ(v1, v2, v3)
                   a = dot((v1-v3),(v1-v3));
                   b = dot((v1-v3),(v1-v2));
                   c = dot((v1-v2),(v1-v2));
                   I111 = log((b+sqrt(a).*sqrt(c))./(b-c-sqrt(c).*sqrt(a-2*b+c)))./(40*sqrt(c))...
                          + log((-b+c+sqrt(c).*sqrt(a-2*b+c)./(-b+sqrt(a).*sqrt(c))))./(40*sqrt(c))...
                          + (sqrt(a).*sqrt(a-2*b+c)-sqrt(c).*sqrt(a-2*b+c))./(120*(a-2*b+c).^(3/2))...
                          + ((2*a-5*b+3*c).*log((a-b+sqrt(a).*sqrt(a-2*b+c).*(c-b+sqrt(c).*sqrt(a-2*b+c)))) ...
                          ./(b-a+sqrt(a).*sqrt(a-2*b+c).*(b-c+sqrt(c).*sqrt(a-2*b+c))))./(120*(a-2*b+c).^(3/2))...
                          +(-sqrt(a).*sqrt(c)+sqrt(a).*sqrt(a-2*b+c))./(120*a.^(3/2))...
                          + ((2*a+b).*log(((b+sqrt(a).*sqrt(c)).*(a-b+sqrt(a).*sqrt(a-2*b+c)))...
                          ./((-b+sqrt(a).*sqrt(c)).*(-a+b+sqrt(a).*sqrt(a-2*b+c)))))./(120*a.^(3/2));
                      
                   I112 = log((b+sqrt(a).*sqrt(c))./(b-c+sqrt(c).*sqrt(a-2*b+c)))./(120*sqrt(c))...
                          + log((a-b+sqrt(a).*sqrt(a-2*b+c))./(-b+sqrt(a).*sqrt(c)))./(120*sqrt(a))...
                          + (-sqrt(a).*sqrt(a-2*b+c)+sqrt(c).*sqrt(a-2*b+c))./(120*(a-2*b+c).^(3/2))...
                          + (2*a-3*b+c).*log((a-b+sqrt(a).*sqrt(a-2*b+c))./(b-c+sqrt(c).*sqrt(a-2*b+c)))./(120*(a-2*b+c).^(3/2))...
                          + (sqrt(a).*sqrt(a-2*b+c)-sqrt(c).*sqrt(a-2*b+c))./(120*(a-2*b+c).^(3/2))...
                          + (a-3*b+2*c).*log((-b+c+sqrt(c).*sqrt(a-2*b+c))./(-a+b+sqrt(a).*sqrt(a-2*b+c)))./(120*(a-2*b+c).^(3/2))...
                          + (-3*sqrt(a).*sqrt(c)+3*sqrt(c).*sqrt(a-2*b+c))./(120*(c).^(3/2))...
                          + (3*b+2*c).*log((-b+c+sqrt(c).*sqrt(a-2*b+c))./(-b+sqrt(a).*sqrt(c)))./(120*c.^(3/2))...
                          + (-3*sqrt(a).*sqrt(c)+3*sqrt(a).*sqrt(a-2*b+c))./(120*(c).^(3/2))...
                          + (2*a+3*b).*log((b+sqrt(a).*sqrt(c))./(-a+b+sqrt(a).*sqrt(a-2*b+c)))./(120*a.^(3/2));
                   
                      I11  = (-log((-b+sqrt(a).*sqrt(c))./(a-b+sqrt(a).*sqrt(a-2*b+c))))./(24*sqrt(a))...
                        + (log((b+sqrt(a).*sqrt(c))./(b-c+sqrt(c).*sqrt(a-2*b+c))))./(24*sqrt(c))...
                        + (-sqrt(a).*sqrt(c)+sqrt(a).*sqrt(a-2*b+c))./(24*a.^(3/2))...
                        + (a+b).*log((b+sqrt(a).*sqrt(c)./(-a+b+sqrt(a).*sqrt(a-2*b+c))))./(24*a.^(3/2))...
                        + (log((a-b+sqrt(a).*sqrt(a-2*b+c))./(b-c+sqrt(c).*sqrt(a-2*b+c))))./(24*sqrt(a-2*b+c))...
                        - (log((b+sqrt(a).*sqrt(c))./(-b+c+sqrt(c).*sqrt(a-2*b+c))))./(12*sqrt(c)) ...
                        + (sqrt(a).*sqrt(a-2*b+c)-sqrt(c).*sqrt(a-2*b+c))./(24*(a-2*b+c).^(3/2))...
                        + ((a-3*b+2*c).*log((-b+c+sqrt(c).*sqrt(a-2*b+c))./(-a+b+sqrt(a).*sqrt(a-2*b+c))))./(24*(a-2*b+c).^(3/2));
end

function [I2] = NearTriangleB(v1, v2, v3)
                   a = dot((v1-v3),(v1-v3));
                   b = dot((v1-v3),(v1-v2));
                   c = dot((v1-v2),(v1-v2));
                   
                   I2   = (log(((a-b+sqrt(a).*sqrt(a-2*b+c)).*(b+sqrt(a).*sqrt(c))./((-b+sqrt(a).*sqrt(c)).*(-a+b+sqrt(a).*sqrt(a-2*b+c))))))./(6*sqrt(a))...
                        + (log(((b+sqrt(a).*sqrt(c)).*(-b+c+sqrt(c).*sqrt(a-2*b+c))./((b-c+sqrt(c).*sqrt(a-2*b+c)).*(-b+sqrt(a).*sqrt(a-2*b+c))))))./(6*sqrt(c))...
                        + (log((a-b+sqrt(a).*sqrt(a-2*b+c)).*(-b+c+sqrt(c).*sqrt(a-2*b+c))./((b-c+sqrt(c).*sqrt(a-2*b+c)).*(-a+b+sqrt(a).*sqrt(a-2*b+c)))))./(6*sqrt(a-2*b+c));
end