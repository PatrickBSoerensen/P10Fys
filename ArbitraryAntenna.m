classdef ArbitraryAntenna
    %ARBRITARYANTENNA MoM method for Arbitary antenna with flat triangle mesh
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        function [p, t] = RemoveEqualPoints(stl)
            %Loading p/t matrixes from STL object
            p = zeros(size(stl.vertices));
            t = stl.faces;
            
            %Looping through vertices
            for i=1:length(stl.vertices)
                % Comparing current vertices to all others
                a = stl.vertices==stl.vertices(i,:);
                % Summing over rows
                a = sum(a,2);
                % Finding rows where all points where the same
                b = a==3;
                % Finding index of points that are the same
                d=find(b);
                % Writing the first point into p matrix as a unqiue point
                p(d(1),:) = stl.vertices(d(1),:);
                % Loops through the points ekvivalent to the first unique
                % one and replaces it in the t matrix
                for n=2:length(d)
                    c = stl.faces == d(n);
                    t(c) = d(1);
                end
            end
            % Moving the row index into a temporary index in column four
            for i=1:length(p)
                p(i,4) = i;
            end
            % Removes any empty points i.e. duplicates
            p( ~any(p(:,1:3),2), : ) = []; 
            
            % Updates vertices in t matrix with new index after removing
            % duplicate points from p
            for i=1:length(p)
                t(t==p(i,4)) = i;
            end
            % Removing the temporary index from p
            p(:,4) = [];
            
            t = sort(t,2);
        end
        
        function [ConnectCell] = Connectivity(p, t)
            %Setting up cell size
            ConnectCell = cell(length(p),2);
            %Looping through points
            for i=1:length(p)
                % Finding point i in t
                a = t==i;
                % Summing rows
                b = sum(a,2);
                % Finding indexes where the current point appears
                b = find(b);
                % Inverting logical indexing i.e. a row where the point
                % was found on was 1,0,0 is turned into 0,1,1
                a(~any(a,2),:) = logical(1-a(~any(a,2),:));
                a = logical(1-a);
                % Inserting unique points, i.e. first row in t is 1,2,3.
                % The first point is 1, therefore the first row in a looks
                % like 0,1,1 Connected then gets 2,3 representing the edges
                % 1-2 and 1-3
                Connected = unique(t(a));
                Connected = Connected(Connected>i);
                ConnectCell{i,1} = Connected;
    
                ConnectCell{i,2} = unique(b);
                % Final ConnectCell has row indices equal to the starting
                % point of an edge.
                % Column 1 contains arrays of endpoints, so
                % ConnectCell{1,1} is an array of end points. For instance
                % ConnectCell{1,1}(1) would return 2 signifiying the edge
                % 1-2, this edge is not contained in ConnectCell{2,1}.
                % Column 2 contains the triangle faces that the point is
                % part of i.e. ConnectCell{1,2} is the array of faces that
                % point 1 is contained in.
            end
        end
        
        function [Atot, Center] = TriangleAreas(p, t)
            TotTri = length(t);
            Center = zeros(TotTri,3);
                        
            for m=1:TotTri
                % Finding center of triangle
                Center(m,:) = sum(p(t(m,:),:))/3;
            end
            % Area of Triangles
            L1 = p(t(:,2),:)-p(t(:,1),:);
            L2 = p(t(:,3),:)-p(t(:,1),:);
        
            Atot = sqrt(sum((cross(L1,L2)).^2,2))/2;
        end
        
        function [SubTri, SubTriArea] = SubTriangles(p, t, Center)
            % Method for creating subtriangles
            TotTri = length(t);
            
            for i=1:TotTri
                % Triangle points index
                n1 = t(i,1);
                n2 = t(i,2);
                n3 = t(i,3);
                % Center of triangle
                M = Center(i,:);
                % Coordinates to triangle points
                r1 = p(n1,:);
                r2 = p(n2,:);
                r3 = p(n3,:);
                % Edges coordinates
                r12=r2-r1;
                r23=r3-r2;
                r13=r3-r1;
                % Subtriangle points on edges
                C1=r1+(1/3)*r12;
                C2=r1+(2/3)*r12;
                C3=r2+(1/3)*r23;
                C4=r2+(2/3)*r23;
                C5=r1+(1/3)*r13;
                C6=r1+(2/3)*r13;
                % Subtriangles centers
                a1=1/3*(C1+C5+r1);
                a2=1/3*(C1+C2+M);
                a3=1/3*(C2+C3+r2);
                a4=1/3*(C2+C3+M);
                a5=1/3*(C3+C4+M);
                a6=1/3*(C1+C5+M);
                a7=1/3*(C5+C6+M);
                a8=1/3*(C4+C6+M);
                a9=1/3*(C4+C6+r3);
                % Saving subtriangles centerpoints
                SubTri(:,:,i)=...
                    [a1 a2 a3 a4 a5 a6 a7 a8 a9];
                % SubTriangles areas
                L1 = r1-C1;
                L2 = C5-C1;
                L3 = M-C1;
                L4 = C2-M;
                L5 = C2-C3;
                L6 = r2-C2;
                L7 = M-C3;
                L8 = M-C4;
                L9 = C6-M;
                L10 = M-C5;
                L11 = C6-r3;
                
                A1 = sqrt(sum((cross(L1,L2)).^2,2))/2;
                A2 = sqrt(sum((cross(L3,L4)).^2,2))/2;
                A3 = sqrt(sum((cross(L5,L6)).^2,2))/2;
                A4 = sqrt(sum((cross(L4,L5)).^2,2))/2;
                A5 = sqrt(sum((cross(L7,L8)).^2,2))/2;
                A6 = sqrt(sum((cross(L2,L3)).^2,2))/2;
                A7 = sqrt(sum((cross(L9,L10)).^2,2))/2;
                A8 = sqrt(sum((cross(L8,L9)).^2,2))/2;
                A9 = sqrt(sum((cross(L9,L11)).^2,2))/2;
                
                SubTriArea(:,:,i)=...
                    [A1 A2 A3 A4 A5 A6 A7 A8 A9];
            end
        end
        
        function I2 = SelfTerm(p, t)
            % Method for calculating self coupling terms
            TotTri = length(t);
            I2 = 1:TotTri;
            % Looping through triangles
            for i=1:TotTri
                % Point index
                n1 = t(i,1);
                n2 = t(i,2);
                n3 = t(i,3);
                % vertices coordinates
                v1 = p(n1,:);
                v2 = p(n2,:);
                v3 = p(n3,:);
                % Intermediate calculations
                a = dot((v1-v3),(v1-v3));
                b = dot((v1-v3),(v1-v2));
                c = dot((v1-v2),(v1-v2));
                d=a-2*b+c;
                
                % More intermediate calculations, expression from Eibert &
                % Hansen
                N1=(a-b+sqrt(a)*sqrt(d))*(b+sqrt(a)*sqrt(c));
                D1=(-a+b+sqrt(a)*sqrt(d))*(-b+sqrt(a)*sqrt(c));
    
                N2=(-b+c+sqrt(c)*sqrt(d))*(b+sqrt(a)*sqrt(c));
                D2=(b-c+sqrt(c)*sqrt(d))*(-b+sqrt(a)*sqrt(c));
    
                N3=(a-b+sqrt(a)*sqrt(d))*(-b+c+sqrt(c)*sqrt(d));
                D3=(b-c+sqrt(c)*sqrt(d))*(-a+b+sqrt(a)*sqrt(d));
                    
                I2(i) = 1/6*(1/sqrt(a)*log(N1/D1) + 1/sqrt(c)*log(N2/D2) + 1/sqrt(d)*log(N3/D3));
                I2(i) = I2(i);
            end
        end
        
        function [EdgeList, Basis, BasisLA, BasisDeriv, BasisNumber, BasisArea] = BasisFunc(p, t, ConnectCell)
            % Method for creating basis functions and saving other
            % information relevant to these
            BasisNumber = 1;

            for i=1:length(p)
                % temp holds a list of point indexes for the faces that the
                % i'th point is a part of.
                temp = t(ConnectCell{i,2},:);
                % nodes is the array of points that the i'th point shares
                % an edge with
                nodes = ConnectCell{i,1};
                
                for n=1:length(nodes)
                    % Finding triangles that match the n'th point.
                    % Note it only checks the second and third point as the
                    % first point is known to be the index from ConnectCell
                    a = temp(:,2:end) == nodes(n);
                    % Summing rows of a and finding the indexes in b where
                    % b is different from 0
                    b = sum(a,2);
                    b = find(b);
                    % Saving the Nodes on opposite sides of the current edge
                    NodesOfInterrest = temp(b,:);
                    % Logical indexing for finding the Node/vertex that
                    % isn't part of the edge
                    c = NodesOfInterrest ~= i;
                    d = logical(NodesOfInterrest ~= nodes(n));
                    c = logical(c.*d);
                    
                    % Triangles that the edge is part of
                    Triangles = temp(b,:);
                    % The points that span the edge
                    EdgePoints = [i,nodes(n)];
                    % The verteces that are part of the triangles but not
                    % part of the edge
                    NotEdgePoints = Triangles(c);
                    
                    % Vector of edge
                    L = p(EdgePoints(1),:)-p(EdgePoints(2),:);
                    
                    % vector of other side in triangle to calculate area
                    % Here the choice of what is the plus and minus
                    % triangle is decided
                    LforAP = p(EdgePoints(1),:)-p(NotEdgePoints(1),:);
                    LforAM = p(EdgePoints(1),:)-p(NotEdgePoints(2),:);
                    
                    % Area calculation for plus and minus triangle
                    % respectively
                    AP =  sqrt(sum((cross(L,LforAP)).^2,2))./2;
                    AM =  sqrt(sum((cross(L,LforAM)).^2,2))./2;
                    
                    % Calculating the length of the edge
                    L = sqrt(sum(L.^2));
        
                    % Start point of the edge
                    EdgeList(BasisNumber, 1) = i;
                    % End point of the edge
                    EdgeList(BasisNumber, 2) = nodes(n);
                    % Plus triangle point
                    EdgeList(BasisNumber, 3) = NotEdgePoints(1);
                    % Minus triangle point
                    EdgeList(BasisNumber, 4) = NotEdgePoints(2);
        
                    % Plus basis function
                    Basis{BasisNumber,1} = @(r) (p(NotEdgePoints(1),:) - r);
                    % Frontterm for basis function
                    BasisLA(BasisNumber,1) =  L./(2*AP);
                    BasisLA(BasisNumber,2) =  L;
                    % Minus basis function
                    Basis{BasisNumber,2} = @(r) (r - p(NotEdgePoints(2),:));
                    % Fronterm for basis function
                    BasisLA(BasisNumber,3) = L./(2*AM);
                    BasisLA(BasisNumber,4) = L;
                    % Plus differentiated basis function
                    BasisDeriv(BasisNumber,1) = -L./AP;
                    % Minus differentiated basis function
                    BasisDeriv(BasisNumber,2) = L./AM;
                    
                    % Saving plus and minus triangle area
                    BasisArea(BasisNumber,1) = AP;
                    BasisArea(BasisNumber,2) = AM;
                
                    BasisNumber = BasisNumber + 1; 
                end
            end
            BasisNumber = BasisNumber - 1;
        end
        
        function [BasisNumber] = EdgeNumbering(EdgeList, Edge)
            for i=1:3
                if isempty(find(sum(EdgeList(:,1:2)==Edge(i,:),2)==2,1))
                    value = find(sum(fliplr(EdgeList(:,1:2))==Edge(i,:),2)==2);
                    BasisNumber(i) = value;
                else
                    value = find(sum(EdgeList(:,1:2)==Edge(i,:),2)==2);
                    BasisNumber(i) = value;
                end
            end
        end
        
        function [RhoP, RhoM, RhoP_, RhoM_] = BasisEvalCenter(t, EdgeList, Basis, Center, SubTri)
            % Looping through triangles
            for m=1:length(t)
                % Three edges of the triangle
                Edge(1,:) = t(m,1:2);
                Edge(2,:) = t(m,2:3);
                Edge(3,1) = t(m,1);
                Edge(3,2) = t(m,3);
                
                % Identifying Edge numbering
                BasisNumber = ArbitraryAntenna.EdgeNumbering(EdgeList, Edge);
                
                %Looping through the 3 edges for a triangle
                for i=1:3
                    % Part of triangle from edge
                    T = EdgeList(BasisNumber(i),1:2);
                    % Getting the plus Basis function that correspond to the
                    % edge
                    BasisP = Basis{BasisNumber(i), 1};
                    % The plus triangle
                    TP = [T EdgeList(BasisNumber(i),3)];
                    TP = sort(TP);
                    % Finding the index of the plus triangle in t
                    t1 = t(:,1)==TP(1);
                    t2 = t(:,2)==TP(2);
                    t3 = t(:,3)==TP(3);
                    TP = find(t1+t2+t3==3);
                    % Evaluating the basis function in center of the
                    % triangle
                    RhoP(BasisNumber(i),:) = BasisP(Center(TP,:));
                    % Evaluating the basis function in centres of
                    % subtriangles
                    SubP = SubTri(1,:,TP);
                    SubP = reshape(SubP, [9,3]);
                    RhoP_(:,:,BasisNumber(i)) = BasisP(SubP);
                    % Getting the minus Basis function that correspond to the
                    % edge
                    BasisM = Basis{BasisNumber(i), 2};
                    % The minus triangle
                    TM = [T EdgeList(BasisNumber(i),4)];
                    % Finding the index of the minus triangle in t
                    TM = sort(TM);
                    t1 = t(:,1)==TM(1);
                    t2 = t(:,2)==TM(2);
                    t3 = t(:,3)==TM(3);
                    TM = find(t1+t2+t3==3);
                    % Evaluating the basis function in center of the
                    % triangle
                    RhoM(BasisNumber(i),:) = BasisM(Center(TM,:));
                    % Evaluating the basis function in centres of
                    % subtriangles
                    SubM = SubTri(1,:,TM);
                    SubM = reshape(SubM, [9,3]);
                    RhoM_(:,:,BasisNumber(i)) = BasisM(SubM);
                end        
            end
        end
        
        function [Z, b, J, a] = MoM(p, t, EdgeList, BasisNumber, BasisLA, A, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, x, y, z)
            % alocating space
            Z = zeros(BasisNumber,BasisNumber)+1i*zeros(BasisNumber,BasisNumber);
            Ei = zeros(size(t));
            Ei(:,1) = x.*exp(-1i*k.*(Center(:,2)));
            Ei(:,2) = y.*exp(-1i*k.*(Center(:,1)));
            Ei(:,3) = z.*exp(-1i*k.*(Center(:,1)));

            % Outer loop over triangles
            for y=1:length(t)
                % Outer triangle edges
                OuterEdge(1,:) = t(y,1:2);
                OuterEdge(2,:) = t(y,2:3);
                OuterEdge(3,1) = t(y,1);
                OuterEdge(3,2) = t(y,3);
                
                % Locating points that complete triangles for basis
                % functions
                EdgeNumberOuter = ArbitraryAntenna.EdgeNumbering(EdgeList, OuterEdge);
                
                % Inner triangle loop
                for h=1:length(t)
                    % Finding Inner triangle edges
                    InnerEdge(1,:) = t(h,1:2);
                    InnerEdge(2,:) = t(h,2:3);
                    InnerEdge(3,1) = t(h,1);
                    InnerEdge(3,2) = t(h,3);
        
                    % Locating points that complete triangles for basis
                    % functions
                    EdgeNumberInner = ArbitraryAntenna.EdgeNumbering(EdgeList, InnerEdge);
                    
                    %Loop over outer triangles basis functions
                    for i = 1:3
                        % outer(1,:) and outer(2,:) are edge points
                        % outer(3,:) are the plus triangle point
                        % outer(4,:) are the minus triangle point
                        outer(1,:) =  p(EdgeList(EdgeNumberOuter(i),1),:);
                        outer(2,:) =  p(EdgeList(EdgeNumberOuter(i),2),:);
                        outer(3,:) =  p(EdgeList(EdgeNumberOuter(i),3),:);
                        outer(4,:) =  p(EdgeList(EdgeNumberOuter(i),4),:);
                        
                        %Finding all the points of the plus and minus
                        %triangle
                        FaceEdgeOuter(1:2) = EdgeList(EdgeNumberOuter(i),1:2);
                        FaceEdgeOuterP = [FaceEdgeOuter EdgeList(EdgeNumberOuter(i),3)];
                        FaceEdgeOuterP = sort(FaceEdgeOuterP);
                        FaceEdgeOuterM = [FaceEdgeOuter EdgeList(EdgeNumberOuter(i),4)];
                        FaceEdgeOuterM = sort(FaceEdgeOuterM);
                        % Finding the index of plus and minus triangles for
                        % the current basis for the outer
                        PlusOuter = find(sum(t==FaceEdgeOuterP,2)==3);
                        MinusOuter = find(sum(t==FaceEdgeOuterM,2)==3);
                        
                        % Intermediate loading of plus and minus functions
                        % evaluated in center points
                        zMP = A(PlusOuter,1).*(RhoP(EdgeNumberOuter(i),:));
                        zMM = A(MinusOuter,2).*(RhoM(EdgeNumberOuter(i),:));
                            
                        % Subtriangles for the inner triangles basis
                        % functions
                        SPO = reshape(SubTri(:,:,PlusOuter),[9,3]);
                        SMO = reshape(SubTri(:,:,MinusOuter),[9,3]);
                            
                        %Inner triangle basis functions loop
                        for j = 1:3
                            % inner(1,:) and inner(2,:) are edge points
                            % inner(3,:) are the plus triangle point
                            % inner(4,:) are the minus triangle point
                            inner(1,:) =  p(EdgeList(EdgeNumberInner(j),1),:);
                            inner(2,:) =  p(EdgeList(EdgeNumberInner(j),2),:);
                            inner(3,:) =  p(EdgeList(EdgeNumberInner(j),3),:);
                            inner(4,:) =  p(EdgeList(EdgeNumberInner(j),4),:);
                           
                            %Finding all the points of the plus and minus
                            %triangle
                            FaceEdgeInner(1:2) = EdgeList(EdgeNumberInner(j),1:2);
                            FaceEdgeInnerP = [FaceEdgeInner EdgeList(EdgeNumberInner(j), 3)];
                            FaceEdgeInnerP = sort(FaceEdgeInnerP);
                            FaceEdgeInnerM = [FaceEdgeInner EdgeList(EdgeNumberInner(j), 4)];
                            FaceEdgeInnerM = sort(FaceEdgeInnerM);
                            PlusInner = find(sum(t==FaceEdgeInnerP,2)==3);
                            MinusInner = find(sum(t==FaceEdgeInnerM,2)==3);
                        
                            % Intermediate loading of plus and minus functions
                            % evaluated in center points
                            zNP = A(PlusInner,1).*(RhoP(EdgeNumberInner(j),:));
                            zNM = A(MinusInner,2).*(RhoM(EdgeNumberInner(j),:));
                            
                            % Subtriangles for the inner triangles basis
                            % functions
                            SPI = reshape(SubTri(:,:,PlusInner),[9,3]);
                            SMI = reshape(SubTri(:,:,MinusInner),[9,3]);
                                                        
                            % Intermediate calculations
                            ppo = sqrt(sum((Center(PlusOuter,:)-SPI).^2,2));
                            mpo = sqrt(sum((Center(MinusOuter,:)-SPI).^2,2));
                            pmo = sqrt(sum((Center(PlusOuter,:)-SMI).^2,2));
                            mmo = sqrt(sum((Center(MinusOuter,:)-SMI).^2,2));
                            
                            ppi = sqrt(sum((Center(PlusInner,:)-SPO).^2,2));
                            mpi = sqrt(sum((Center(MinusInner,:)-SPO).^2,2));
                            pmi = sqrt(sum((Center(PlusInner,:)-SMO).^2,2));
                            mmi = sqrt(sum((Center(MinusInner,:)-SMO).^2,2));
                            
                            % greens for different couplings
                            gPPo = exp(1i.*k.*ppo)./ppo;
                            gMPo = exp(1i.*k.*mpo)./mpo;
                            gPMo = exp(1i.*k.*pmo)./pmo;
                            gMMo = exp(1i.*k.*mmo)./mmo;
                        
                            gPPi = exp(1i.*k.*ppi)./ppi;
                            gMPi = exp(1i.*k.*mpi)./mpi;
                            gPMi = exp(1i.*k.*pmi)./pmi;
                            gMMi = exp(1i.*k.*mmi)./mmi;
                        
                            if PlusOuter==PlusInner
                                % greens for the self term
                                g = I2(PlusOuter);
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *((dot(zMP, zNP)/4-1/k^2) * g)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            else
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMP, zNP)/4-1/k^2) * gPPo)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMP, zNP)/4-1/k^2) * gPPi)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            end
                            if PlusOuter==MinusInner
                                % greens for the self term
                                g = I2(PlusOuter);
                                
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *((dot(zMP, zNM)/4-1/k^2) * g)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            else
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMP, zNM)/4+1/k^2) * gPMo )...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                                
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMP, zNM)/4+1/k^2) * gPMi )...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            end
                            if MinusOuter==PlusInner        
                                % greens for the self term
                                g = I2(MinusOuter);
                                
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *((dot(zMM, zNP)/4-1/k^2) * g)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            else
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMM, zNP)/4+1/k^2) * gMPo)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMM, zNP)/4+1/k^2) * gMPi)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            end
                            if MinusOuter==MinusInner
                                % greens for the self term
                                g = I2(MinusOuter);
                                
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *((dot(zMM, zNM)/4-1/k^2) * g)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            else
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMM, zNM)/4-1/k^2) * gMMo)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            
                                Z(EdgeNumberOuter(i), EdgeNumberInner(j)) = ...
                                (BasisLA(EdgeNumberOuter(i),2)*BasisLA(EdgeNumberInner(j),4))/(4*pi)...
                                *sum(1/9*(dot(zMM, zNM)/4-1/k^2) * gMMi)...
                                + Z(EdgeNumberOuter(i), EdgeNumberInner(j));
                            end
                        end
                        % Intermediate variables for b calculation
                        b1 = RhoP_(:,:,EdgeNumberOuter(i));
                        b2 = RhoM_(:,:,EdgeNumberOuter(i));
                        b3 = BasisLA(EdgeNumberOuter(i),2)/2;
                        
                        b(EdgeNumberOuter(i),:) = sum(sum(b1.*Ei(PlusOuter,:)/9,2)+sum(b2.*Ei(MinusOuter,:)/9,2).*b3);
                    end    
                end
            end
<<<<<<< HEAD
=======
            % Final b calculation
            b = dot(Ei,bCollect,2);
            Upper = triu(Z,1);
            Lower = tril(Z,-1);
            Z = Z + Upper.' + Lower.';
>>>>>>> 56fefb2e9f2103f02173acd0c324b4ca9d130450
            % Z\b is a newer faster version of inv(Z)*b
            a = Z\b;
            J = a.*(BasisLA(:,1).*RhoP+BasisLA(:,3).*RhoM);
        end
        
        function [Eyz, Exz, Exy, xrange, yrange] = EFieldXY(center, pJ, w, k, mu,...
                xmin, xmax, ymin, ymax, steps, LA, RhoP_, RhoM_, a, SubTri, t, EdgeList)
            
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            zrange = linspace(ymin, ymax, steps);
            
            Eyz = zeros(steps,steps);
            Exz = zeros(steps,steps);
            Exy = zeros(steps,steps);
            
            B = (-1i.*w.*mu);
            
            for i=1:steps
                for m=1:steps
                    for j=1:length(center)
                        Edge(1,:) = t(j,1:2);
                        Edge(2,:) = t(j,2:3);
                        Edge(3,1) = t(j,1);
                        Edge(3,2) = t(j,3);
                
                        EdgeNumber = ArbitraryAntenna.EdgeNumbering(EdgeList, Edge);
                
                        for y=1:3
                            FaceEdge(1:2) = EdgeList(EdgeNumber(y),1:2);
                            FaceEdgeP = [FaceEdge EdgeList(EdgeNumber(y),3)];
                            FaceEdgeP = sort(FaceEdgeP);
                            FaceEdgeM = [FaceEdge EdgeList(EdgeNumber(y),4)];
                            FaceEdgeM = sort(FaceEdgeM);
                            PlusEdge = find(sum(t==FaceEdgeP,2)==3);
                            MinusEdge = find(sum(t==FaceEdgeM,2)==3);
                        
                            rx = (0 - pJ(PlusEdge));
                            ry = (yrange(m) - pJ(PlusEdge));
                            rz = (zrange(i) - pJ(PlusEdge));
                            r = sqrt(rx^2+ry^2+rz^2);
                        
                            g = B.*exp(1i.*k.*r)./(r);
                            temp = sum(1/9 * RhoP_(:,:,EdgeNumber(y)).* ...
                            exp(-1i.*k.*sqrt(sum(reshape(SubTri(:,:,PlusEdge),[9,3]).^2,2))));
                                               
                            Eyz(i,m) = Eyz(i,m) + g * a(PlusEdge) * LA(PlusEdge,1)*temp(1);
                            
                            rx = (xrange(i) - pJ(PlusEdge));
                            ry = (0 - pJ(PlusEdge));
                            rz = (zrange(m) - pJ(PlusEdge));
                            r = sqrt(rx^2+ry^2+rz^2);
                        
                            g = B.*exp(1i.*k.*r)./(r);
                            Exz(i,m) = Exz(i,m) + g * a(PlusEdge) * LA(PlusEdge,1)*temp(2);
                            
                            rx = (xrange(i) - pJ(PlusEdge));
                            ry = (yrange(m) - pJ(PlusEdge));
                            rz = (0 - pJ(PlusEdge));
                            r = sqrt(rx^2+ry^2+rz^2);
                        
                            g = B.*exp(1i.*k.*r)./(r);
                            Exy(i,m) = Exy(i,m) + g * a(PlusEdge) * LA(PlusEdge,1)*temp(3);
                            
                            temp = sum(1/9 * RhoM_(:,:,EdgeNumber(y)).* ...
                            exp(-1i.*k.*sqrt(sum(reshape(SubTri(:,:,MinusEdge),[9,3]).^2,2))));
                            rx = (0 - pJ(MinusEdge));
                            ry = (yrange(m) - pJ(MinusEdge));
                            rz = (zrange(i) - pJ(MinusEdge));
                            r = sqrt(rx^2+ry^2+rz^2);
                        
                            g = B.*exp(1i.*k.*r)./(r);
                            Eyz(i,m) = Eyz(i,m) + g * a(MinusEdge) * LA(MinusEdge,1) * temp(1);
                            
                            rx = (xrange(i) - pJ(MinusEdge));
                            ry = (0 - pJ(MinusEdge));
                            rz = (zrange(m) - pJ(MinusEdge));
                            r = sqrt(rx^2+ry^2+rz^2);
                        
                            g = B.*exp(1i.*k.*r)./(r);
                            Exz(i,m) = Exz(i,m) + g * a(MinusEdge) * LA(MinusEdge,1)*temp(2);

                            rx = (xrange(i) - pJ(MinusEdge));
                            ry = (yrange(m) - pJ(MinusEdge));
                            rz = (0 - pJ(MinusEdge));
                            r = sqrt(rx^2+ry^2+rz^2);
                        
                            g = B.*exp(1i.*k.*r)./(r);
                            Exy(i,m) = Exy(i,m) + g * a(MinusEdge) * LA(MinusEdge,1)*temp(3);
                        end
                    end
                end
            end
            
        end
        
        function [Exy, Exz, Eyz, xrange, yrange, zrange, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = EField(Center, w, k, mu, J,...
                xmin, xmax, zmin, zmax, ymin, ymax, steps, Area, LA, RhoP_, RhoM_, a, SubTri, t, EdgeList)
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            zrange = linspace(zmin, zmax, steps);
            
            Exy = zeros(steps,steps);
            Exyx = Exy; Exyy = Exy; Exyz = Exy;
            Exz = zeros(steps,steps);
            Exzx = Exz; Exzy = Exz; Exzz = Exz;
            Eyz = zeros(steps,steps);
            Eyzx = Eyz; Eyzy = Eyz; Eyzz = Eyz;
                        
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
                
                B = (-1i.*w.*mu);
                
                for i=1:length(Center)
                    Edge(1,:) = t(i,1:2);
                    Edge(2,:) = t(i,2:3);
                    Edge(3,1) = t(i,1);
                    Edge(3,2) = t(i,3);
                
                    EdgeNumber = ArbitraryAntenna.EdgeNumbering(EdgeList, Edge);
                
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
               
                            g = B.*exp(1i.*k.*r)./(r);
                            temp = sum(1/9 * RhoP_(:,:,PlusEdge).* ...
                            exp(1i.*k.*sqrt(sum(reshape(SubTri(:,:,PlusEdge),[9,3]).^2,2))));
                        
                            Exyx = Exyx + g .* a(PlusEdge) .* LA(PlusEdge,1) * temp(1);
                            Exyy = Exyy + g .* a(PlusEdge) .* LA(PlusEdge,1) * temp(2);
                            Exyz = Exyz + g .* a(PlusEdge) .* LA(PlusEdge,1) * temp(3);          
                       
                            %Minus
                            rr(:,1) = rx(MinusEdge,:).';
                            rr(:,2) = ry(MinusEdge,:).';
                            rr(:,3) = rz(MinusEdge,:).';
                            r = sqrt((rz(MinusEdge,:)).^2+(ry(MinusEdge,:).').^2+(rx(MinusEdge,:)).^2);
               
                            g = B.*exp(1i.*k.*r)./(r);
                            temp = sum( 1/9 *RhoM_(:,:,MinusEdge) .* ...
                            exp(1i.*k.*sqrt(sum(reshape(SubTri(:,:,MinusEdge),[9,3]).^2,2))));
                        
                            Exyx = Exyx + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(1);
                            Exyy = Exyy + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(2);
                            Exyz = Exyz + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(3);
                       
                        elseif j==2
                            rr(:,1) = rx(PlusEdge,:).';
                            rr(:,2) = ry(PlusEdge,:).';
                            rr(:,3) = rz(PlusEdge,:).';
                            r = sqrt((rz(PlusEdge,:).').^2+(ry(PlusEdge,:)).^2+(rx(PlusEdge,:)).^2);
                
                            g = B.*exp(1i.*k.*r)./(r);
                             temp = sum(1/9 * RhoP_(:,:,PlusEdge).* ...
                            exp(1i.*k.*sqrt(sum(reshape(SubTri(:,:,PlusEdge),[9,3]).^2,2))));
                        
                            Exzx = Exzx +   g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                            .* temp(1);
                            Exzy = Exzy +   g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                            .* temp(2);
                            Exzz = Exzz +   g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                            .* temp(3);
                       
                        
                            %Minus
                            rr(:,1) = rx(MinusEdge,:).';
                            rr(:,2) = ry(MinusEdge,:).';
                            rr(:,3) = rz(MinusEdge,:).';
                            r = sqrt((rz(MinusEdge,:).').^2+(ry(MinusEdge,:)).^2+(rx(MinusEdge,:)).^2);
               
                            g = B.*exp(1i.*k.*r)./(r);
                            temp = sum( 1/9 * RhoM_(:,:,MinusEdge).* ...
                            exp(1i.*k.*sqrt(sum(reshape(SubTri(:,:,MinusEdge),[9,3]).^2,2))));
                             
                            Exzx = Exzx + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(1);
                        
                            Exzy = Exzy + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(2);
                        
                            Exzz = Exzz + g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(3);
                       
                        else
                            rr(:,1) = rx(PlusEdge,:).';
                            rr(:,2) = ry(PlusEdge,:).';
                            rr(:,3) = rz(PlusEdge,:).';
                            r = sqrt((rz(PlusEdge,:)).^2+(ry(PlusEdge,:).').^2+(rx(PlusEdge,:)).^2);
               
                            g = B.*exp(1i.*k.*r)./(r);
                            temp = sum(1/9 *RhoP_(:,:,PlusEdge).* ...
                            exp(1i.*k.*sqrt(sum(reshape(SubTri(:,:,PlusEdge),[9,3]).^2,2))));
                            
                            Eyzx = Eyzx +  g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                            .* temp(1);
                        Eyzy = Eyzy +  g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                            .* temp(2);
                        Eyzz = Eyzz +  g .* a(PlusEdge) .* LA(PlusEdge,1) ...
                            .* temp(3);
                        
                            %Minus
                            rr(:,1) = rx(MinusEdge,:).';
                            rr(:,2) = ry(MinusEdge,:).';
                            rr(:,3) = rz(MinusEdge,:).';
                            r = sqrt((rz(MinusEdge,:)).^2+(ry(MinusEdge,:).').^2+(rx(MinusEdge,:)).^2);
               
                            g = B.*exp(1i.*k.*r)./(r);
                            
                            temp = sum(1/9 * RhoM_(:,:,MinusEdge) .* ...
                            exp(1i.*k.*sqrt(sum(reshape(SubTri(:,:,MinusEdge),[9,3]).^2,2))));
                            
                            Eyzx = Eyzx +  g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(1);
                            Eyzy = Eyzy +  g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(2);
                            Eyzz = Eyzz +  g .* a(MinusEdge) .* LA(MinusEdge,1) ...
                            .* temp(3);
                       
                        end    
                    end    
                end
            end
        end
        
        function [Exy, Exz, Eyz, xrange, yrange, zrange] = EFieldAlt(pJ, w, k, mu,...
                xmin, xmax, zmin, zmax, ymin, ymax, steps, LA, RhoP_, RhoM_, a, SubTri)
            
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            zrange = linspace(zmin, zmax, steps);
            
            Exy = zeros(steps,steps);
            Exz = zeros(steps,steps);
            Eyz = zeros(steps,steps);
            
            rx = (xrange-pJ(:,1));
            ry = (yrange-pJ(:,2));
            rz = (zrange-pJ(:,3));
            
            rzz= (0-pJ(:,3));
            ryy = (0-pJ(:,2));
            rxx = (0-pJ(:,1));
               
            B = (-1i.*w.*mu);
            
            for i=1:steps
                rr(:,1) = rx(:,i).';
                rr(:,2) = ry(:,i).';
                rr(:,3) = rz(:,i).';
                
                Rx = sqrt((rz(:,i)).^2+(ry(:,i).').^2+(rxx(:,i)).^2);
                Ry = sqrt((rz(:,i).').^2+(ryy(:,i)).^2+(rx(:,i)).^2);
                Rz = sqrt((rzz(:,i)).^2+(ry(:,i).').^2+(rx(:,i)).^2);
                
                gx = B.*exp(1i.*k.*Rx)./(Rx);
                gy = B.*exp(1i.*k.*Ry)./(Ry);
                gz = B.*exp(1i.*k.*Rz)./(Rz);
                                                        
                Exy = Exy + 1/9 * gz .* a(:) .* LA(:,1) ...
                    .* sum(dot(RhoP_(:,:,:), ...
                    exp(1i.*k.*SubTri(:,:,:)),2));
                                 
                Exy = Exy + 1/9 * gz .* a(:) .* LA(:,1) ...
                    .* sum(dot(RhoM_(:,:,:), ...
                    exp(1i.*k.*SubTri(:,:,:)),2));
                                 
                Exz = Exz + 1/9 * gy .* a(:) .* LA(:,1) ...
                    .* sum(dot(RhoP_(:,:,:), ...
                    exp(1i.*k.*SubTri(:,:,:)),2));
                
                Exz = Exz + 1/9 * gy .* a(:) .* LA(:,1) ...
                    .* sum(dot(RhoM_(:,:,:), ...
                    exp(1i.*k.*SubTri(:,:,:)),2));
                
                Eyz = Eyz + 1/9 * gx .* a(:) .* LA(:,1) ...
                    .* sum(dot(RhoP_(:,:,:), ...
                    exp(1i.*k.*SubTri(:,:,:)),2));
                
                Eyz = Eyz + 1/9 * gx .* a(:) .* LA(:,1) ...
                    .* sum(dot(RhoM_(:,:,:), ...
                    exp(1i.*k.*SubTri(:,:,:)),2));
   
            end 
            
        end 
    end    
end

