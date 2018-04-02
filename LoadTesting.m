%% konstants
eps0=8.854187817*10^-12; %F/m
mu0=4*pi*10^-7; %N/A^2
c=1/sqrt(eps0*mu0); %m/s

f=146.5*10^6; 
lambda=c/f;
w=2*pi*f;
k=w/c;
%% load STL file into matlab
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
%% Calculating areas for Simplex Coords
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
%% Basis Function setup
BasisNumberOuter = 1;

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

for f=1:length(faces)
    
    v1 = UV(faces(f,1),:);
    v2 = UV(faces(f,2),:);
    v3 = UV(faces(f,3),:);
    
    [I111(f), I112(f), I11(f)] = NearTriangleZ(v1, v2, v3);
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
    clear BasisNumberOuter
    OuterBasis(1,:) = faces(f,1:2);
    OuterBasis(2,:) = faces(f,2:3);
    OuterBasis(3,1) = faces(f,1);
    OuterBasis(3,2) = faces(f,3);
    
    v1 = UV(faces(f,1),:);
    v2 = UV(faces(f,2),:);
    v3 = UV(faces(f,3),:);
    v(1,:) = v1;
    v(2,:) = v2;
    v(3,:) = v3;
    
    for i=1:3
        if isempty(find(sum(EdgeList(:,1:2)==OuterBasis(i,:),2)==2,1))
            value = find(sum(fliplr(EdgeList(:,1:2))==OuterBasis(i,:),2)==2);
            BasisNumberOuter(i) = find(sum(fliplr(EdgeList(:,1:2))==OuterBasis(i,:),2)==2);
        else
            value = find(sum(EdgeList(:,1:2)==OuterBasis(i,:),2)==2);
            BasisNumberOuter(i) = value;
        end
    v(3+i,:) = UV(EdgeList(value,4),:);
    end
    
    for h=1:length(faces)
        clear BasisNumberInner
        InnerBasis(1,:) = faces(h,1:2);
        InnerBasis(2,:) = faces(h,2:3);
        InnerBasis(3,1) = faces(h,1);
        InnerBasis(3,2) = faces(h,3);
       
        u1 = UV(faces(h,1),:);
        u2 = UV(faces(h,2),:);
        u3 = UV(faces(h,3),:);
        
        u(1,:) = u1;
        u(2,:) = u2;
        u(3,:) = u3;
     
        for i=1:3
            if isempty(find(sum(EdgeList(:,1:2)==InnerBasis(i,:),2)==2,1))
                value = find(sum(fliplr(EdgeList(:,1:2))==InnerBasis(i,:),2)==2);
                BasisNumberInner(i) =  value;
            else
                value = find(sum(EdgeList(:,1:2)==InnerBasis(i,:),2)==2);
                BasisNumberInner(i) = value;
            end
        u(3+i,:) = UV(EdgeList(value,4),:);
        end
             for sing=1:6
                 if ~isempty(find(sum(u(1:6,:) == v(sing,:),2)==3,1))
                     a = find(sum(u(1:6,:) == v(sing,:),2)==3)
                    Near(sing) = a(1);
                 else
                     Near(sing) = 0;
                 end
             end
             Nearsum = sum(sum(Near));
             clear Near
            if Nearsum > 0
                Z(BasisNumberOuter, BasisNumberInner) = I111(f)+I112(f)+I11(f) + Z(BasisNumberOuter, BasisNumberInner);
            
                Z(BasisNumberOuter, BasisNumberInner) = I111(h)+I112(h)+I11(h) + Z(BasisNumberOuter, BasisNumberInner);
            else
                for i = 1:3
                    for j = 1:3
                BaseMP = Basis{BasisNumberOuter(i),1};
                BaseMM = Basis{BasisNumberOuter(i),2};
                BaseDMP = BasisDeriv(BasisNumberOuter(i),1);
                BaseDMM = BasisDeriv(BasisNumberOuter(i),2);
            
                BaseNP = Basis{BasisNumberInner(j),1};
                BaseNM = Basis{BasisNumberInner(j),2};
                BaseDNP = BasisDeriv(BasisNumberInner(j),1);
                BaseDNM = BasisDeriv(BasisNumberInner(j),2);
                
                %Determine if + og - is correctly placed
                z1M = @(q) wp*(1/4*BaseMP(UV(EdgeList(BasisNumberOuter(i),1),:))+1/k^2)*g(UV(EdgeList(BasisNumberOuter(i),1),:),q);
                z2M = @(q) wp*(1/4*BaseMP(UV(EdgeList(BasisNumberOuter(i),2),:))+1/k^2)*g(UV(EdgeList(BasisNumberOuter(i),2),:),q);
                z3M = @(q) wp*(1/4*BaseMP(UV(EdgeList(BasisNumberOuter(i),3),:))+1/k^2)*g(UV(EdgeList(BasisNumberOuter(i),3),:),q);
                z4M = @(q) wp*(1/4*BaseMM(UV(EdgeList(BasisNumberOuter(i),4),:))-1/k^2)*g(UV(EdgeList(BasisNumberOuter(i),4),:),q);
                outer(1,:) =  UV(EdgeList(BasisNumberOuter(i),1),:);
                outer(2,:) =  UV(EdgeList(BasisNumberOuter(i),2),:);
                outer(3,:) =  UV(EdgeList(BasisNumberOuter(i),3),:);
                outer(4,:) =  UV(EdgeList(BasisNumberOuter(i),4),:);
                
                z1N = @(q) wp*(1/4*BaseNP(UV(EdgeList(BasisNumberInner(j),1),:))+1/k^2)*g(UV(EdgeList(BasisNumberInner(j),1),:),q);
                z2N = @(q) wp*(1/4*BaseNP(UV(EdgeList(BasisNumberInner(j),2),:))+1/k^2)*g(UV(EdgeList(BasisNumberInner(j),2),:),q);
                z3N = @(q) wp*(1/4*BaseNP(UV(EdgeList(BasisNumberInner(j),3),:))+1/k^2)*g(UV(EdgeList(BasisNumberInner(j),3),:),q);
                z4N = @(q) wp*(1/4*BaseNM(UV(EdgeList(BasisNumberInner(j),4),:))-1/k^2)*g(UV(EdgeList(BasisNumberInner(j),4),:),q);
                inner(1,:) =  UV(EdgeList(BasisNumberInner(j),1),:);
                inner(2,:) =  UV(EdgeList(BasisNumberInner(j),2),:);
                inner(3,:) =  UV(EdgeList(BasisNumberInner(j),3),:);
                inner(4,:) =  UV(EdgeList(BasisNumberInner(j),4),:);
                
                %Clean up loop and functions
                Z(BasisNumberOuter(i), BasisNumberInner(j)) = ...
                        (BasisLA(BasisNumberOuter(i),2)*BasisLA(BasisNumberInner(j),4))/(4*pi)...
                        * dot(z1M(inner(1,:)),z1N(outer(1,:)))...
                        + dot(z1M(inner(2,:)),z2N(outer(2,:)))...
                        + dot(z1M(inner(3,:)),z3N(outer(3,:)))...
                        + dot(z1M(inner(4,:)),z4N(outer(4,:)))...
                        + dot(z2M(inner(1,:)),z1N(outer(1,:)))...
                        + dot(z2M(inner(2,:)),z2N(outer(2,:)))...
                        + dot(z2M(inner(3,:)),z3N(outer(3,:)))...
                        + dot(z2M(inner(4,:)),z4N(outer(4,:)))...
                        + dot(z3M(inner(1,:)),z1N(outer(1,:)))...
                        + dot(z3M(inner(2,:)),z2N(outer(2,:)))...
                        + dot(z3M(inner(3,:)),z3N(outer(3,:)))...
                        + dot(z3M(inner(4,:)),z4N(outer(4,:)))...
                        + dot(z4M(inner(1,:)),z1N(outer(1,:)))...
                        + dot(z4M(inner(2,:)),z2N(outer(2,:)))...
                        + dot(z4M(inner(3,:)),z3N(outer(3,:)))...
                        + dot(z4M(inner(4,:)),z4N(outer(4,:)))...
                        + Z(BasisNumberOuter(i), BasisNumberInner(j));
                    end
                end
            end
    end
end
hej=55;

b(m) = sum(BasisLA(m,2)/2*wp*(BaseMP(mv1)+BaseMP(mv2)+BaseMP(mv3))*Ei);

zinv = inv(Z);
xtesst = b/Z;
            
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