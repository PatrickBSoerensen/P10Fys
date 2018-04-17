classdef ArbitraryAntenna
    %ARBRITARYANTENNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        function [p, t] = RemoveEqualPoints(stl)
            p = zeros(size(stl.vertices));
            t = stl.faces;
            
            for i=1:length(stl.vertices)
                a = stl.vertices==stl.vertices(i,:);
                a = sum(a,2);
                b = a==3;
                d=find(b);
                p(d(1),:) = stl.vertices(d(1),:);
                
                for n=2:length(d)
                    c=stl.faces == d(n);
                    t(c) = d(1);
                end
            end
            
            for i=1:length(p)
                p(i,4) = i;
            end
            
            p( ~any(p(:,1:3),2), : ) = []; 
            
            for i=1:length(p)
                t(t==p(i,4)) = i;
            end
            p(:,4) = [];
        end
        
        function [ConnectCell] = GibsonConnect(p, t)
            
            ConnectCell = cell(length(p),2);
            for i=1:length(p)
                a = t==i;
                b = sum(a,2);
                b = find(b);
                a(~any(a,2),:) = logical(1-a(~any(a,2),:));
                a = logical(1-a);
    
                Connected = unique(t(a));
                Connected = Connected(Connected>i);
                ConnectCell{i,1} = Connected;
    
                ConnectCell{i,2} = unique(b);
            end
        end
        
        function [Area, Atot, Center] = TriangleAreas(p, t)
            TotTri = length(t);
            Area = zeros(TotTri,3);
            Center = zeros(TotTri,3);
                        
            for m=1:TotTri
                Center(m,:) = sum(p(t(m,:),:))/3;
    
                L1 = p(t(m,1),:)-Center(m,:);
                L2 = p(t(m,3),:)-p(t(m,1),:);
                Area(m,1) = sqrt(sum((cross(L1,L2)).^2,2))/2;
        
                L2 = p(t(m,2),:)-p(t(m,1),:);
                Area(m,2) = sqrt(sum((cross(L1,L2)).^2,2))/2;
        
                L1 = p(t(m,2),:)-Center(m,:);
                L2 = p(t(m,3),:)-p(t(m,2),:);
                Area(m,3) = sqrt(sum((cross(L1,L2)).^2,2))/2;
            end
    
            L1 = p(t(:,2),:)-p(t(:,1),:);
            L2 = p(t(:,3),:)-p(t(:,1),:);
        
            Atot = sqrt(sum((cross(L1,L2)).^2,2))/2;
            Area = Area./Atot;
        end
        
        function [SubTri, Integral] = SubTriangles(p, t, Center)
            TotTri = length(t);
            Integral = 1:TotTri;
            
            for i=1:TotTri
                n1 = t(i,1);
                n2 = t(i,2);
                n3 = t(i,3); 
                M = Center(i,:);
                r1 = p(n1,:);
                r2 = p(n2,:);
                r3 = p(n3,:);
                r12=r2-r1;
                r23=r3-r2;
                r13=r3-r1;
                C1=r1+(1/3)*r12;
                C2=r1+(2/3)*r12;
                C3=r2+(1/3)*r23;
                C4=r2+(2/3)*r23;
                C5=r1+(1/3)*r13;
                C6=r1+(2/3)*r13;
                a1=1/3*(C1+C5+r1);
                a2=1/3*(C1+C2+M);
                a3=1/3*(C2+C3+r2);
                a4=1/3*(C2+C3+M);
                a5=1/3*(C3+C4+M);
                a6=1/3*(C1+C5+M);
                a7=1/3*(C5+C6+M);
                a8=1/3*(C4+C6+M);
                a9=1/3*(C4+C6+r3);
                SubTri(:,:,i)=...
                    [a1 a2 a3 a4 a5 a6 a7 a8 a9];
            end
        end
        
        function I2 = SelfTerm(p, t)
            TotTri = length(t);
            I2 = 1:TotTri;
            
            for i=1:TotTri
                n1 = t(i,1);
                n2 = t(i,2);
                n3 = t(i,3); 
                v1 = p(n1,:);
                v2 = p(n2,:);
                v3 = p(n3,:);
                a = dot((v1-v3),(v1-v3));
                b = dot((v1-v3),(v1-v2));
                c = dot((v1-v2),(v1-v2));
                d=a-2*b+c;
    
                N1=(a-b+sqrt(a)*sqrt(d))*(b+sqrt(a)*sqrt(c));
                D1=(-a+b+sqrt(a)*sqrt(d))*(-b+sqrt(a)*sqrt(c));
    
                N2=(-b+c+sqrt(c)*sqrt(d))*(b+sqrt(a)*sqrt(c));
                D2=(b-c+sqrt(c)*sqrt(d))*(-b+sqrt(a)*sqrt(c));
    
                N3=(a-b+sqrt(a)*sqrt(d))*(-b+c+sqrt(c)*sqrt(d));
                D3=(b-c+sqrt(c)*sqrt(d))*(-a+b+sqrt(a)*sqrt(d));
                    
                I2(i) = 1/6*(1/sqrt(a)*log(N1/D1) + 1/sqrt(c)*log(N2/D2) + 1/sqrt(d)*log(N3/D3));
                I2(i) = 4*I2(i);
            end
        end
        
        function [EdgeList, Basis, BasisLA, BasisDeriv, BasisNumber, BasisArea] = BasisFunc(p, t, ConnectCell)
        BasisNumber = 1;

        for i=1:length(p)
            temp = t(ConnectCell{i,2},:);
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
        
                L = p(EdgePoints(1),:)-p(EdgePoints(2),:);
        
                LforAP = p(EdgePoints(1),:)-p(NotEdgePoints(1),:);
                LforAM = p(EdgePoints(1),:)-p(NotEdgePoints(2),:);
            
                AP =  sqrt(sum((cross(L,LforAP)).^2,2))./2;
                AM =  sqrt(sum((cross(L,LforAM)).^2,2))./2;
    
                L = sqrt(sum(L.^2));
        
                %Central point of a cluster
                EdgeList(BasisNumber, 1) = i;
                %Connection point
                EdgeList(BasisNumber, 2) = nodes(n);
                %Plus triangle point
                EdgeList(BasisNumber, 3) = NotEdgePoints(1);
                %Minus triangle point
                EdgeList(BasisNumber, 4) = NotEdgePoints(2);
        
                %Plus
                Basis{BasisNumber,1} = @(r) (p(NotEdgePoints(1),:) - r);
                BasisLA(BasisNumber,1) =  L./(2*AP);
                BasisLA(BasisNumber,2) =  L;
                %Minus
                Basis{BasisNumber,2} = @(r) (r - p(NotEdgePoints(2),:));
                BasisLA(BasisNumber,3) = L./(2*AM);
                BasisLA(BasisNumber,4) = L;
                %Plus
                BasisDeriv(BasisNumber,1) = -L./AP;
                %Minus
                BasisDeriv(BasisNumber,2) = L./AM;
               
                BasisArea(BasisNumber,1) = AP;
                BasisArea(BasisNumber,2) = AM;
                
                BasisNumber = BasisNumber + 1; 
            end
        end
        BasisNumber = BasisNumber - 1;
        end
        
        function [RhoP, RhoM, RhoP_, RhoM_] = BasisEvalCenter(t, EdgeList, Basis, Center, SubTri)
            for m=1:length(t)
                Edge(1,:) = t(m,1:2);
                Edge(2,:) = t(m,2:3);
                Edge(3,1) = t(m,1);
                Edge(3,2) = t(m,3);
    
                for i=1:3
                    if isempty(find(sum(EdgeList(:,1:2)==Edge(i,:),2)==2,1))
                        value = find(sum(fliplr(EdgeList(:,1:2))==Edge(i,:),2)==2);
                        BasisNumber(i) = value;
                    else
                        value = find(sum(EdgeList(:,1:2)==Edge(i,:),2)==2);
                        BasisNumber(i) = value;
                    end
                end
                
                for i=1:3
                    T = EdgeList(BasisNumber(i),1:2);
                    BasisP = Basis{BasisNumber(i), 1};
                    TP = [T EdgeList(BasisNumber(i),3)];
                    TP = sort(TP);
                    t1 = t(:,1)==TP(1);
                    t2 = t(:,2)==TP(2);
                    t3 = t(:,3)==TP(3);
                    TP = find(t1+t2+t3==3);
                    RhoP(BasisNumber(i),:) = BasisP(Center(TP,:));
                    SubP = SubTri(1,:,TP);
                    SubP = reshape(SubP, [3,9]).';
                    RhoP_(:,:,BasisNumber(i)) = BasisP(SubP);
                    BasisM = Basis{BasisNumber(i), 2};
                    TM = [T EdgeList(BasisNumber(i),4)];
                    TM = sort(TM);
                    t1 = t(:,1)==TM(1);
                    t2 = t(:,2)==TM(2);
                    t3 = t(:,3)==TM(3);
                    TM = find(t1+t2+t3==3);
                    SubM = SubTri(1,:,TM);
                    SubM = reshape(SubM, [3,9]).';
                    RhoM_(:,:,BasisNumber(i)) = BasisM(SubM);
                    RhoM(BasisNumber(i),:) = BasisM(Center(TM,:));
                end        
            end
        end
        
        function [Z, b, J, a] = MoM(p, t, EdgeList, BasisNumber, BasisLA, A, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, x, y, z)
        Z = zeros(BasisNumber,BasisNumber)+1i*zeros(BasisNumber,BasisNumber);
        Ei = zeros(BasisNumber,3);
        Ei(:,1) = x; Ei(:,2) = y; Ei(:,3) = z;
        b = zeros(BasisNumber,1);
    
        for y=1:length(t)
            OuterEdge(1,:) = t(y,1:2);
            OuterEdge(2,:) = t(y,2:3);
            OuterEdge(3,1) = t(y,1);
            OuterEdge(3,2) = t(y,3);
     
            for i=1:3
                if isempty(find(sum(EdgeList(:,1:2)==OuterEdge(i,:),2)==2,1))
                    value = find(sum(fliplr(EdgeList(:,1:2))==OuterEdge(i,:),2)==2);
                    EdgeNumberOuter(i) = value;
                else
                    value = find(sum(EdgeList(:,1:2)==OuterEdge(i,:),2)==2);
                    EdgeNumberOuter(i) = value;
                end
            end
            
            for h=1:length(t)
                InnerEdge(1,:) = t(h,1:2);
                InnerEdge(2,:) = t(h,2:3);
                InnerEdge(3,1) = t(h,1);
                InnerEdge(3,2) = t(h,3);
        
                for i=1:3
                    if isempty(find(sum(EdgeList(:,1:2)==InnerEdge(i,:),2)==2,1))
                        value = find(sum(fliplr(EdgeList(:,1:2))==InnerEdge(i,:),2)==2);
                        EdgeNumberInner(i) =  value;
                    else
                        value = find(sum(EdgeList(:,1:2)==InnerEdge(i,:),2)==2);
                        EdgeNumberInner(i) = value;
                    end
                end
        
                for i = 1:3
                    outer(1,:) =  p(EdgeList(EdgeNumberOuter(i),1),:);
                    outer(2,:) =  p(EdgeList(EdgeNumberOuter(i),2),:);
                    outer(3,:) =  p(EdgeList(EdgeNumberOuter(i),3),:);
                    outer(4,:) =  p(EdgeList(EdgeNumberOuter(i),4),:);
                    
                        FaceEdgeOuter(1:2) = EdgeList(EdgeNumberOuter(i),1:2);
                        FaceEdgeOuterP = [FaceEdgeOuter EdgeList(EdgeNumberOuter(i),3)];
                        FaceEdgeOuterP = sort(FaceEdgeOuterP);
                        FaceEdgeOuterM = [FaceEdgeOuter EdgeList(EdgeNumberOuter(i),4)];
                        FaceEdgeOuterM = sort(FaceEdgeOuterM);
                        PlusOuter = find(sum(t==FaceEdgeOuterP,2)==3);
                        MinusOuter = find(sum(t==FaceEdgeOuterM,2)==3);
                
                        zMP = A(PlusOuter,1).*(1/4*RhoP(EdgeNumberOuter(i),:));
                        zMM = A(MinusOuter,2).*(1/4*RhoM(EdgeNumberOuter(i),:));
                
                    for j = 1:3
                        inner(1,:) =  p(EdgeList(EdgeNumberInner(j),1),:);
                        inner(2,:) =  p(EdgeList(EdgeNumberInner(j),2),:);
                        inner(3,:) =  p(EdgeList(EdgeNumberInner(j),3),:);
                        inner(4,:) =  p(EdgeList(EdgeNumberInner(j),4),:);
                            
                        FaceEdgeInner(1:2) = EdgeList(EdgeNumberInner(j),1:2);
                        FaceEdgeInnerP = [FaceEdgeInner EdgeList(EdgeNumberInner(j), 3)];
                        FaceEdgeInnerP = sort(FaceEdgeInnerP);
                        FaceEdgeInnerM = [FaceEdgeInner EdgeList(EdgeNumberInner(j), 4)];
                        FaceEdgeInnerM = sort(FaceEdgeInnerM);
                        PlusInner = find(sum(t==FaceEdgeInnerP,2)==3);
                        MinusInner = find(sum(t==FaceEdgeInnerM,2)==3);
                
                        zNP = A(PlusInner,1).*(1/4*RhoP(EdgeNumberInner(j),:));
                        zNM = A(MinusInner,2).*(1/4*RhoM(EdgeNumberInner(j),:));
                        
                        SPI = reshape(SubTri(:,:,PlusInner),[3,9]).';
                        SMI = reshape(SubTri(:,:,MinusInner),[3,9]).';
                        
                        g = I2(y)/4;
                        
                        pp = sqrt(sum((Center(PlusOuter,:)-SPI).^2,2));
                        mp = sqrt(sum((Center(MinusOuter,:)-SPI).^2,2));
                        pm = sqrt(sum((Center(PlusOuter,:)-SMI).^2,2));
                        mm = sqrt(sum((Center(MinusOuter,:)-SMI).^2,2));
                        
                        gPP = exp(1i.*k.*pp)./pp;
                        gMP = exp(1i.*k.*mp)./mp;
                        gPM = exp(1i.*k.*pm)./pm;
                        gMM = exp(1i.*k.*mm)./mm;
                        
                        if PlusOuter==PlusInner
                            Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(36*pi)...
                            *((dot(zMP, zNP)-1/k^2) * g)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        else
                            Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                            *sum((dot(zMP, zNP)-1/k^2) * gPP)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        end
                        if PlusOuter==MinusInner 
                           Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(36*pi)...
                            *((dot(zMP, zNM)-1/k^2) * g)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        else
                            Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                            *sum((dot(zMP, zNM)+1/k^2) * gPM )...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        end
                        if MinusOuter==PlusInner        
                        Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(36*pi)...
                            *((dot(zMM, zNP)-1/k^2) * g)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        else
                            Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                            *sum((dot(zMM, zNP)+1/k^2) * gMP)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        end
                        if MinusOuter==MinusInner
                           Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(36*pi)...
                            *((dot(zMM, zNM)-1/k^2) * g)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                        else
                            Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                            (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                            *sum((dot(zMM, zNM)-1/k^2) * gMM)...
                            + Z(EdgeNumberOuter(i), EdgeNumberInner(j));           
                        end
                    end
                    
                    b1 = reshape(RhoP_(:,:,EdgeNumberOuter(i)),[9,3]);
                    b2 = reshape(RhoM_(:,:,EdgeNumberOuter(i)),[9,3]);
                    b3 = BasisLA(EdgeNumberOuter(i),2);
                    
                    bCollect(EdgeNumberOuter(i),:) = sum((b1+b2).*b3./2);
                end
            end
        end
        b = dot(Ei,bCollect,2) + b;
        % Z\b is a newer faster version of inv(Z)*b
        a = Z\b;
        J = a.*(RhoP+RhoM);
        end
        
        function [Exy, Exz, Eyz, xrange, yrange, zrange] = EField(Center, w, k, mu, J,...
                xmin, xmax, zmin, zmax, ymin, ymax, steps, Area, LA, RhoP_, RhoM_, a, SubTri, t, EdgeList)
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            zrange = linspace(zmin, zmax, steps);
            
            Exy = zeros(steps,steps);
            Exz = zeros(steps,steps);
            Eyz = zeros(steps,steps);
            
            for j=1:3
                if j == 1
                    rx = (xrange-Center(:,1));
                    ry = (yrange-Center(:,2));
                    rz = (0-Center(:,3));
                elseif j==2
                    rx = (xrange-Center(:,1));
                    ry = (0-Center(:,2));
                    rz = (zrange-Center(:,3));
                else
                    rx = (0-Center(:,1));
                    ry = (yrange-Center(:,2));
                    rz = (zrange-Center(:,3));
                end
                
                I = eye(500, 500);
                B = (1i.*w.*mu);
                
             for i=1:length(Center)
                Edge(1,:) = t(i,1:2);
                Edge(2,:) = t(i,2:3);
                Edge(3,1) = t(i,1);
                Edge(3,2) = t(i,3);
                
                for h=1:3
                    if isempty(find(sum(EdgeList(:,1:2)==Edge(h,:),2)==2,1))
                        value = find(sum(fliplr(EdgeList(:,1:2))==Edge(h,:),2)==2);
                        EdgeNumber(h) = value;
                    else
                        value = find(sum(EdgeList(:,1:2)==Edge(h,:),2)==2);
                        EdgeNumber(h) = value;
                    end
                end
                
                for y=1:3
                        FaceEdge(1:2) = EdgeList(EdgeNumber(y),1:2);
                        FaceEdgeP = [FaceEdge EdgeList(EdgeNumber(y),3)];
                        FaceEdgeP = sort(FaceEdgeP);
                        FaceEdgeM = [FaceEdge EdgeList(EdgeNumber(y),4)];
                        FaceEdgeM = sort(FaceEdgeM);
                        PlusEdge = find(sum(t==FaceEdgeP,2)==3);
                        MinusEdge = find(sum(t==FaceEdgeM,2)==3);
                
                 if j == 1
                     %Plus
                    rr(:,1) = rx(PlusEdge,:).';
                    rr(:,2) = ry(PlusEdge,:).';
                    rr(:,3) = rz(PlusEdge,:).';
                    r = sqrt((rz(PlusEdge,:)).^2+(ry(PlusEdge,:).').^2+(rx(PlusEdge,:)).^2);
               
                    g = B.*exp(1i.*k.*r)./(4.*pi.*r);
%                   G = g.*(I.*(1+1i./(k*r)-1./((k*r).^2)) - ...
%                          (rr*rr.')./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
                    
                    Exy = Exy + g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                        .* sum(dot(reshape(RhoP_(:,:,PlusEdge),[9,3]), ...
                           exp(1i.*k.*reshape(SubTri(:,:,PlusEdge),[9,3])),2));
                    
                    %Minus
                    rr(:,1) = rx(MinusEdge,:).';
                    rr(:,2) = ry(MinusEdge,:).';
                    rr(:,3) = rz(MinusEdge,:).';
                    r = sqrt((rz(MinusEdge,:)).^2+(ry(MinusEdge,:).').^2+(rx(MinusEdge,:)).^2);
               
                    g = B.*exp(1i.*k.*r)./(4.*pi.*r);
                    
                    Exy = Exy + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                        .* sum(dot(reshape(RhoP_(:,:,MinusEdge),[9,3]), ...
                           exp(1i.*k.*reshape(SubTri(:,:,MinusEdge),[9,3])),2));                    
                elseif j==2
                    rr(:,1) = rx(PlusEdge,:).';
                    rr(:,2) = ry(PlusEdge,:).';
                    rr(:,3) = rz(PlusEdge,:).';
                    r = sqrt((rz(PlusEdge,:).').^2+(ry(PlusEdge,:)).^2+(rx(PlusEdge,:)).^2);
               
                    g = B.*exp(1i.*k.*r)./(4.*pi.*r);
%                   G = g.*(I.*(1+1i./(k*r)-1./((k*r).^2)) - ...
%                          (rr*rr.')./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
                    
                    Exz = Exz + g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                        .* sum(dot(reshape(RhoP_(:,:,PlusEdge),[9,3]), ...
                           exp(1i.*k.*reshape(SubTri(:,:,PlusEdge),[9,3])),2));
                    
                    %Minus
                    rr(:,1) = rx(MinusEdge,:).';
                    rr(:,2) = ry(MinusEdge,:).';
                    rr(:,3) = rz(MinusEdge,:).';
                    r = sqrt((rz(MinusEdge,:).').^2+(ry(MinusEdge,:)).^2+(rx(MinusEdge,:)).^2);
               
                    g = B.*exp(1i.*k.*r)./(4.*pi.*r);
                    
                    Exz = Exz + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                        .* sum(dot(reshape(RhoP_(:,:,MinusEdge),[9,3]), ...
                           exp(1i.*k.*reshape(SubTri(:,:,MinusEdge),[9,3])),2));    
                else
                    rr(:,1) = rx(PlusEdge,:).';
                    rr(:,2) = ry(PlusEdge,:).';
                    rr(:,3) = rz(PlusEdge,:).';
                    r = sqrt((rz(PlusEdge,:)).^2+(ry(PlusEdge,:).').^2+(rx(PlusEdge,:)).^2);
               
                    g = B.*exp(1i.*k.*r)./(4.*pi.*r);
%                   G = g.*(I.*(1+1i./(k*r)-1./((k*r).^2)) - ...
%                          (rr*rr.')./(r.^2).*(1+3i./(k*r)-3./((k*r).^2)));
                    
                    Eyz = Eyz + g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                        .* sum(dot(reshape(RhoP_(:,:,PlusEdge),[9,3]), ...
                           exp(1i.*k.*reshape(SubTri(:,:,PlusEdge),[9,3])),2));
                       
                    %Minus
                    rr(:,1) = rx(MinusEdge,:).';
                    rr(:,2) = ry(MinusEdge,:).';
                    rr(:,3) = rz(MinusEdge,:).';
                    r = sqrt((rz(MinusEdge,:)).^2+(ry(MinusEdge,:).').^2+(rx(MinusEdge,:)).^2);
               
                    g = B.*exp(1i.*k.*r)./(4.*pi.*r);
                    
                    Eyz = Eyz + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                        .* sum(dot(reshape(RhoP_(:,:,MinusEdge),[9,3]), ...
                           exp(1i.*k.*reshape(SubTri(:,:,MinusEdge),[9,3])),2));
                end
                end
             end
            end
        end
    end
end

