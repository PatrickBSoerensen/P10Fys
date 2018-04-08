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
%% Basis Function simplex setup
BasisNumberOuter = 1;
r = @(alpha, beta, v1, v2, v3) (1-alpha-beta)*v1+alpha*v2+beta*v3;
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
        
        %Plus
        BasisSimp{BasisNumberOuter,1} = @(r) (UV(NotEdgePoints(1),1) - r);
        BasisSimp{BasisNumberOuter,2} = @(r) (UV(NotEdgePoints(1),2) - r);
        BasisSimp{BasisNumberOuter,3} = @(r) (UV(NotEdgePoints(1),3) - r);
        %Minus
        BasisSimp{BasisNumberOuter,4} = @(r) (r - UV(NotEdgePoints(2),1));
        BasisSimp{BasisNumberOuter,5} = @(r) (r - UV(NotEdgePoints(2),2));
        BasisSimp{BasisNumberOuter,6} = @(r) (r - UV(NotEdgePoints(2),3));
  
        BasisNumberOuter = BasisNumberOuter + 1;
    end
end

%% pre analytic calculations
for f=1:length(faces)
    v1 = UV(faces(f,1),:);
    v2 = UV(faces(f,2),:);
    v3 = UV(faces(f,3),:);
    
    I2(f) = NearTriangleFaces(v1, v2, v3);
end

for f=1:length(faces)
    v1 = UV(faces(f,1),:);
    v2 = UV(faces(f,2),:);
    v3 = UV(faces(f,3),:);
    for permut=1:3
    [I111(f,permut), I112(f,permut), I11(f,permut)] = NearTriangleZ(v1, v2, v3, permut);
    end
end


    %% MoM faces loop
    NearAnalyticFixed = 1;
    BasisNumberOuter = BasisNumberOuter-1;
    loader = 0;
    Z = zeros(BasisNumberOuter,BasisNumberOuter);
%     b = 1:BasisNumberOuter-1;
    g = @(p,q) exp(1i.*k.*sqrt(sum((p-q).^2)))./(sqrt(sum((p-q).^2)));
    
    gg = @(r1,r2,r3,rm1,rm2,rm3)...
        (exp(1i.*k.*sqrt((r1-rm1+0.0001).^2+(r2-rm2+0.0001).^2+(r3-rm3+0.0001).^2))./(4.*pi.*sqrt((r1-rm1+0.0001).^2+(r2-rm2+0.0001).^2+(r3-rm3+0.0001).^2))...
        -1/sqrt((r1-rm1+0.0001).^2+(r2-rm2+0.0001).^2+(r3-rm3+0.0001).^2))...
        +1/sqrt((r1-rm1+0.0001).^2+(r2-rm2+0.0001).^2+(r3-rm3+0.0001).^2);
    Ei = 1;
    faces = sort(faces,2);
    
for f=1:length(faces)
    clear BasisNumberOuter
    OuterBasis(1,:) = faces(f,1:2);
    OuterBasis(2,:) = faces(f,2:3);
    OuterBasis(3,1) = faces(f,1);
    OuterBasis(3,2) = faces(f,3);
    
    v(1,:) = UV(faces(f,1),:);
    v(2,:) = UV(faces(f,2),:);
    v(3,:) = UV(faces(f,3),:);
     
    for i=1:3
        if isempty(find(sum(EdgeList(:,1:2)==OuterBasis(i,:),2)==2,1))
            value = find(sum(fliplr(EdgeList(:,1:2))==OuterBasis(i,:),2)==2);
            BasisNumberOuter(i) = value;
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
       
        u(1,:) = UV(faces(h,1),:);
        u(2,:) = UV(faces(h,2),:);
        u(3,:) = UV(faces(h,3),:);
        
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
        
        for i = 1:3
                BaseMP = Basis{BasisNumberOuter(i),1};
                BaseMM = Basis{BasisNumberOuter(i),2};
                BaseDMP = BasisDeriv(BasisNumberOuter(i),1);
                BaseDMM = BasisDeriv(BasisNumberOuter(i),2);
                
                outer(1,:) =  UV(EdgeList(BasisNumberOuter(i),1),:);
                outer(2,:) =  UV(EdgeList(BasisNumberOuter(i),2),:);
                outer(3,:) =  UV(EdgeList(BasisNumberOuter(i),3),:);
                outer(4,:) =  UV(EdgeList(BasisNumberOuter(i),4),:);
            for j = 1:3
                BaseNP = Basis{BasisNumberInner(j),1};
                BaseNM = Basis{BasisNumberInner(j),2};
                BaseDNP = BasisDeriv(BasisNumberInner(j),1);
                BaseDNM = BasisDeriv(BasisNumberInner(j),2);
                
                nearOuter = EdgeList(BasisNumberOuter(i),:);
                nearInner = EdgeList(BasisNumberInner(j),:);
                
                inner(1,:) =  UV(EdgeList(BasisNumberInner(j),1),:);
                inner(2,:) =  UV(EdgeList(BasisNumberInner(j),2),:);
                inner(3,:) =  UV(EdgeList(BasisNumberInner(j),3),:);
                inner(4,:) =  UV(EdgeList(BasisNumberInner(j),4),:);
                for sing=1:4
                    if ~isempty(find(nearInner == nearOuter(sing),1))
                        %Overlap of faces
                        a = find(nearInner == nearOuter(sing),1);
                        Near(sing) = a(1);
                    else
                        Near(sing) = 0;
                    end
                end
             Nearsum = sum(Near);
             clear Near
            
            if Nearsum > 0
                %Responsible for 11580 matrix elements
                if NearAnalyticFixed
                    Z(BasisNumberOuter, BasisNumberInner) = I111(f)+I112(f)+I11(f) + Z(BasisNumberOuter, BasisNumberInner);
                    Z(BasisNumberOuter, BasisNumberInner) = I111(h)+I112(h)+I11(h) + Z(BasisNumberOuter, BasisNumberInner);
                else
                    BaseMP1 = BasisSimp{BasisNumberOuter(i),1};
                    BaseMP2 = BasisSimp{BasisNumberOuter(i),2};
                    BaseMP3 = BasisSimp{BasisNumberOuter(i),3};
                    BaseMM1 = BasisSimp{BasisNumberOuter(i),4};
                    BaseMM2 = BasisSimp{BasisNumberOuter(i),5};
                    BaseMM3 = BasisSimp{BasisNumberOuter(i),6};
                    
                    BaseNP1 = BasisSimp{BasisNumberInner(j),1};
                    BaseNP2 = BasisSimp{BasisNumberInner(j),2};
                    BaseNP3 = BasisSimp{BasisNumberInner(j),3};
                    
                    BaseNM1 = BasisSimp{BasisNumberInner(j),4};
                    BaseNM2 = BasisSimp{BasisNumberInner(j),5};
                    BaseNM3 = BasisSimp{BasisNumberInner(j),6};
                    
                    r1 = @(alpha1, beta1) (1-alpha1-beta1)*outer(1,1) + alpha1*outer(2,1) + beta1*outer(3,1);
                    r2 = @(alpha1, beta1) (1-alpha1-beta1)*outer(1,2) + alpha1*outer(2,2) + beta1*outer(3,2);
                    r3 = @(alpha1, beta1) (1-alpha1-beta1)*outer(1,3) + alpha1*outer(2,3) + beta1*outer(3,3);
                    
                    rmark1 = @(alpha2, beta2) (1-alpha2-beta2)*inner(1,1) + alpha2*inner(2,1) + beta2*inner(3,1);
                    rmark2 = @(alpha2, beta2) (1-alpha2-beta2)*inner(1,2) + alpha2*inner(2,2) + beta2*inner(3,2);
                    rmark3 = @(alpha2, beta2) (1-alpha2-beta2)*inner(1,3) + alpha2*inner(2,3) + beta2*inner(3,3);
                    
                    gg = @(r1,r2,r3,rmark1,rmark2,rmark3) 1./sqrt((r1-rmark1+0.0001).^2+(r2-rmark2+0.0001).^2+(r3-rmark3+0.0001).^2);
                    
                    betamax = @(alpha) alpha - 1;
                    testInner = @(alpha1,beta1) integral2(@(alpha2,beta2) BaseNP1(r(alpha2, beta2, inner(1,1),inner(2,1),inner(3,1)))...
                        .*gg(r1(alpha1,beta1),r2(alpha1,beta1),r3(alpha1,beta1),rmark1(alpha2,beta2),rmark2(alpha2,beta2),rmark2(alpha2,beta2)), 0, 1, 0, betamax);
                    test3 = @(alpha1,beta1) -integral2(@(alpha2,beta2) gg(r1(alpha1,beta1),r2(alpha1,beta1),r3(alpha1,beta1),rmark1(alpha2,beta2),rmark2(alpha2,beta2),rmark2(alpha2,beta2))...
                        , 0, 1, 0, betamax);
                    
                    test = integral2(@(alpha1,beta1) BaseMP1(r(alpha1, beta1, outer(1,1),outer(2,1),outer(3,1))).*testInner(alpha1,beta1), 0, 1, 0, betamax);
                    
                    test4 = integral2(@(alpha1,beta1) test3(alpha2,beta2)/k^2, 0, 1, 0, betamax);
                    
                    Z(BasisNumberOuer, BasisNumberInner) = (BasisLA(BasisNumberOuter(i),1)*BasisLA(BasisNumberInner(j),3))...
                        *(integral2(@(alpha1, beta1)...
                        (1/4*dot(BaseMP(r(alpha1, beta1)),...
                        BaseNP(rmark(alpha2, beta2)))-1/k^2)...
                        .*gg(r(alpha1, beta1)), 0, 1, 0, betamax));...
%                         +integral2(@(alpha, beta)integral2(@(alpha, beta)(1/4*dot(BaseMP(r),BaseNM(rmark))+1/k^2)*gg(r), 0, 1, 0, betamax), 0, 1, 0, betamax)...
%                         +integral2(@(alpha, beta)integral2(@(alpha, beta)(1/4*dot(BaseMM(r),BaseNP(rmark))+1/k^2)*gg(r), 0, 1, 0, betamax), 0, 1, 0, betamax)...
%                         +integral2(@(alpha, beta)integral2(@(alpha, beta)(1/4*dot(BaseMM(r),BaseNM(rmark))-1/k^2)*gg(r), 0, 1, 0, betamax), 0, 1, 0, betamax));
                end
            else
                %Responsible for 145236 matrix elements
                
                FaceEdgeOuter(1:2) = EdgeList(BasisNumberOuter(i),1:2);
                FaceEdgeOuterP = [FaceEdgeOuter  EdgeList(BasisNumberOuter(i),3)];
                FaceEdgeOuterP = sort(FaceEdgeOuterP);
                FaceEdgeOuterM = [FaceEdgeOuter  EdgeList(BasisNumberOuter(i),4)];
                FaceEdgeOuterM = sort(FaceEdgeOuterM);
                PlusOuter = find(sum(faces==FaceEdgeOuterP,2)==3);
                MinusOuter = find(sum(faces==FaceEdgeOuterM,2)==3);
                
                z1Mp = A(PlusOuter,:).*(1/4*BaseMP(UV(EdgeList(BasisNumberOuter(i),1),:)));
                z2Mp = A(PlusOuter,:).*(1/4*BaseMP(UV(EdgeList(BasisNumberOuter(i),2),:)));
                z1Mm = A(MinusOuter,:).*(1/4*BaseMM(UV(EdgeList(BasisNumberOuter(i),1),:)));
                z2Mm = A(MinusOuter,:).*(1/4*BaseMM(UV(EdgeList(BasisNumberOuter(i),2),:)));
                z3M = A(PlusOuter,:).*(1/4*BaseMP(UV(EdgeList(BasisNumberOuter(i),3),:)));
                z4M = A(MinusOuter,:).*(1/4*BaseMM(UV(EdgeList(BasisNumberOuter(i),4),:)));
                
                FaceEdgeInner(1:2) = EdgeList(BasisNumberInner(j),1:2);
                FaceEdgeInnerP = [FaceEdgeInner  EdgeList(BasisNumberInner(j),3)];
                FaceEdgeInnerP = sort(FaceEdgeInnerP); 
                FaceEdgeInnerM = [FaceEdgeInner  EdgeList(BasisNumberInner(j),4)];
                FaceEdgeInnerM = sort(FaceEdgeInnerM);
                PlusInner = find(sum(faces==FaceEdgeInnerP,2)==3);
                MinusInner = find(sum(faces==FaceEdgeInnerM,2)==3);
                
                z1Np = A(PlusInner,:).*(1/4*BaseNP(UV(EdgeList(BasisNumberInner(j),1),:)));
                z2Np = A(PlusInner,:).*(1/4*BaseNP(UV(EdgeList(BasisNumberInner(j),2),:)));
                z1Nm = A(MinusInner,:).*(1/4*BaseNM(UV(EdgeList(BasisNumberInner(j),1),:)));
                z2Nm = A(MinusInner,:).*(1/4*BaseNM(UV(EdgeList(BasisNumberInner(j),2),:)));
                z3N = A(PlusInner,:).*(1/4*BaseNP(UV(EdgeList(BasisNumberInner(j),3),:)));
                z4N = A(MinusInner,:).*(1/4*BaseNM(UV(EdgeList(BasisNumberInner(j),4),:)));
                                
                g11 = g(outer(1,:),inner(1,:));
                g12 = g(outer(1,:),inner(2,:));
                g13 = g(outer(1,:),inner(3,:));
                g14 = g(outer(1,:),inner(4,:));
                g22 = g(outer(2,:),inner(2,:));
                g23 = g(outer(2,:),inner(3,:));
                g24 = g(outer(2,:),inner(4,:));
                g33 = g(outer(3,:),inner(3,:));
                g34 = g(outer(3,:),inner(4,:));
                g44 = g(outer(4,:),inner(4,:));
                
                %Clean up loop and functions
                Z(BasisNumberOuter(i), BasisNumberInner(j)) = ...
                        (BasisLA(BasisNumberOuter(i),2)*BasisLA(BasisNumberInner(j),4))/(4*pi)...
                        *((dot(z1Mp, z1Np)-1/k^2)*g11 ...
                        + (dot(z1Mp, z2Np)-1/k^2)*g12 ...
                        + (dot(z1Mp, z3N)-1/k^2) *g13 ...
                        + (dot(z1Mp, z4N)+1/k^2) *g14...
                        + (dot(z1Mp, z1Nm)+1/k^2)*g11...
                        + (dot(z1Mp, z2Nm)+1/k^2)*g12...
                        ...
                        + (dot(z2Mp, z1Np)-1/k^2)*g12 ...
                        + (dot(z2Mp, z2Np)-1/k^2)*g22 ...
                        + (dot(z2Mp, z3N)-1/k^2)*g23 ...
                        + (dot(z2Mp, z4N)+1/k^2)*g24 ...
                        + (dot(z2Mp, z1Nm)+1/k^2)*g12 ...
                        + (dot(z2Mp, z2Nm)+1/k^2)*g22 ...
                        ...
                        + (dot(z3M, z1Np)-1/k^2)*g13 ...
                        + (dot(z3M, z2Np)-1/k^2)*g23 ...
                        + (dot(z3M, z3N)-1/k^2)*g33 ...
                        + (dot(z3M, z4N)+1/k^2)*g34 ...
                        + (dot(z3M, z1Nm)+1/k^2)*g13 ...
                        + (dot(z3M, z2Nm)+1/k^2)*g23 ...
                        ...
                        + (dot(z4M, z1Np)-1/k^2)*g14 ...
                        + (dot(z4M, z2Np)-1/k^2)*g24 ...
                        + (dot(z4M, z3N)-1/k^2) *g34...
                        + (dot(z4M, z4N)+1/k^2) *g44...
                        + (dot(z4M, z1Nm)+1/k^2) *g14...
                        + (dot(z4M, z2Nm)+1/k^2) *g24...
                        ...
                        + (dot(z1Mp, z1Np)-1/k^2)*g11 ...
                        + (dot(z1Mp, z2Np)-1/k^2)*g12 ...
                        + (dot(z1Mp, z3N)-1/k^2)*g13 ...
                        + (dot(z1Mp, z4N)+1/k^2)*g14 ...
                        + (dot(z1Mp, z1Nm)+1/k^2)*g11 ...
                        + (dot(z1Mp, z2Nm)+1/k^2)*g12 ...
                        ...
                        + (dot(z2Mm, z1Np)-1/k^2)*g12 ...
                        + (dot(z2Mm, z2Np)-1/k^2)*g22 ...
                        + (dot(z2Mm, z3N)-1/k^2) *g23...
                        + (dot(z2Mm, z4N)+1/k^2) *g24...
                        + (dot(z2Mm, z1Nm)+1/k^2) *g12...
                        + (dot(z2Mm, z2Nm)+1/k^2)) *g22...
                        ...
                        + Z(BasisNumberOuter(i), BasisNumberInner(j));
                    end
            end
        b(BasisNumberOuter(i)) = BasisLA(BasisNumberOuter(i),2);
        b1(BasisNumberOuter(i),:) = BaseMP(UV(EdgeList(BasisNumberOuter(i),1),:));
        b2(BasisNumberOuter(i),:) = BaseMM(UV(EdgeList(BasisNumberOuter(i),1),:));
        b3(BasisNumberOuter(i),:) = BaseMP(UV(EdgeList(BasisNumberOuter(i),2),:));
        b4(BasisNumberOuter(i),:) = BaseMM(UV(EdgeList(BasisNumberOuter(i),2),:));
        b5(BasisNumberOuter(i),:) = BaseMP(UV(EdgeList(BasisNumberOuter(i),3),:));
        b6(BasisNumberOuter(i),:) = BaseMM(UV(EdgeList(BasisNumberOuter(i),4),:));
        end
    end
end
%%
b = b(:,1)./2.*A(PlusOuter,:).*Ei.*(sqrt(b1(:,1).^2+b1(:,2).^2+b1(:,3).^2)...
    +sqrt(b2(:,1).^2+b2(:,2).^2+b2(:,3).^2)+sqrt(b3(:,1).^2+b3(:,2).^2+b3(:,3).^2)...
    +sqrt(b4(:,1).^2+b4(:,2).^2+b4(:,3).^2)+sqrt(b5(:,1).^2+b5(:,2).^2+b5(:,3).^2)...
    +sqrt(b6(:,1).^2+b6(:,2).^2+b6(:,3).^2));
% Z\b is a newer faster version of inv(Z)*b
xtesst = Z\b(:,1);
            
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

function [I111, I112, I11] = NearTriangleZ(v1, v2, v3, permut)
                   if permut == 1
                   a = dot((v1-v3),(v1-v3));
                   b = dot((v1-v3),(v1-v2));
                   c = dot((v1-v2),(v1-v2));
                   elseif permut == 2
                   a = dot((v2-v3),(v2-v3));
                   b = dot((v2-v3),(v2-v1));
                   c = dot((v2-v1),(v2-v1));
                   elseif permut == 3
                   a = dot((v3-v1),(v3-v1));
                   b = dot((v3-v2),(v3-v1));
                   c = dot((v3-v2),(v3-v2));
                   end    
                   
                   I111 = log((b+sqrt(a).*sqrt(c))./(b-c-sqrt(c).*sqrt(a-2*b+c)))./(40*sqrt(c))...
                          + log((-b+c+sqrt(c).*sqrt(a-2*b+c)./(-b+sqrt(a).*sqrt(c))))./(40*sqrt(c))...
                          + (sqrt(a).*sqrt(a-2*b+c)-sqrt(c).*sqrt(a-2*b+c))./(60*(a-2*b+c).^(3/2))...
                          + ((2*a-5*b+3*c).*log((a-b+sqrt(a).*sqrt(a-2*b+c).*(c-b+sqrt(c).*sqrt(a-2*b+c)))) ...
                          ./(b-a+sqrt(a).*sqrt(a-2*b+c).*(b-c+sqrt(c).*sqrt(a-2*b+c))))./(120*(a-2*b+c).^(3/2))...
                          +(-sqrt(a).*sqrt(c)+sqrt(a).*sqrt(a-2*b+c))./(60*a.^(3/2))...
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

function I2 = NearTriangleFaces(v1, v2, v3)
                   a = dot((v1-v3),(v1-v3));
                   b = dot((v1-v3),(v1-v2));
                   c = dot((v1-v2),(v1-v2));
                   
                   I2   = (log(((a-b+sqrt(a).*sqrt(a-2*b+c)).*(b+sqrt(a).*sqrt(c))./((-b+sqrt(a).*sqrt(c)).*(-a+b+sqrt(a).*sqrt(a-2*b+c))))))./(6*sqrt(a))...
                        + (log(((b+sqrt(a).*sqrt(c)).*(-b+c+sqrt(c).*sqrt(a-2*b+c))./((b-c+sqrt(c).*sqrt(a-2*b+c)).*(-b+sqrt(a).*sqrt(a-2*b+c))))))./(6*sqrt(c))...
                        + (log((a-b+sqrt(a).*sqrt(a-2*b+c)).*(-b+c+sqrt(c).*sqrt(a-2*b+c))./((b-c+sqrt(c).*sqrt(a-2*b+c)).*(-a+b+sqrt(a).*sqrt(a-2*b+c)))))./(6*sqrt(a-2*b+c));
end