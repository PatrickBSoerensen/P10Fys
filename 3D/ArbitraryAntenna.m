 classdef ArbitraryAntenna
    %ARBRITARYANTENNA MoM method for Arbitary antenna with flat triangle mesh
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        function [p, t] = RemoveDuplicatePoints(stl)
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
        
        function [p, t] = RemoveDuplicatePointsPT(vertices, faces)
            %Loading p/t matrixes from STL object
            p = zeros(size(vertices));
            t = faces;
            
            %Looping through vertices
            for i=1:length(vertices)
                % Comparing current vertices to all others
                a = vertices==vertices(i,:);
                % Summing over rows
                a = sum(a,2);
                % Finding rows where all points where the same
                b = a==3;
                % Finding index of points that are the same
                d=find(b);
                % Writing the first point into p matrix as a unqiue point
                p(d(1),:) = vertices(d(1),:);
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
            L1 = p(t(:,1),:)-p(t(:,2),:);
            L2 = p(t(:,1),:)-p(t(:,3),:);
            Atot = sqrt(sum((cross(L1,L2)).^2,2))/2;
        end
        
        function [Center, SubTriRet] = CenterLift(Center, SubTri, R, Lift)
            ForRe = size(SubTri);
            for k=1:length(Center)
                    CenterToMove=Center(k,:);
                    SubTriToMove=reshape(SubTri(:,:,k), 3, [])';
                    SubTriRet(:,:,k) = reshape(SubTri(:,:,k), 3, [])';
                if Lift
                if Center(k,2) >=0.05 
                    % part of spherical surface
                    CenterOfSphere = [0, 0.05, 0];
                    Normal = CenterToMove-CenterOfSphere;
                    n=sqrt(Normal(1)^2+Normal(2)^2+Normal(3)^2);
                    Normal = Normal/n;
                    Center(k,:) = CenterOfSphere+Normal*R;
                    
                    Normal = SubTriToMove-CenterOfSphere;
                    n=sqrt(Normal(:,1).^2+Normal(:,2).^2+Normal(:,3).^2);
                    Normal = Normal./n;
                    SubTriRet(:,:,k) = CenterOfSphere+Normal*R;
                elseif Center(k,2) <=-0.05
                    % part of spherical surface
                    CenterOfSphere = [0, -0.05, 0];
                    Normal = CenterToMove-CenterOfSphere;
                    n=sqrt(Normal(1)^2+Normal(2)^2+Normal(3)^2);
                    Normal = Normal/n;
                    Center(k,:) = CenterOfSphere+Normal*R;
                    
                    Normal = SubTriToMove-CenterOfSphere;
                    n=sqrt(Normal(:,1).^2+Normal(:,2).^2+Normal(:,3).^2);
                    Normal = Normal./n;
                    SubTriRet(:,:,k) = CenterOfSphere+Normal*R;
                else
                    % Part of Cylinder
                    CenterOfCylinder = [0, CenterToMove(:,2), 0];
                    Normal = CenterToMove-CenterOfCylinder;
                    n=sqrt(Normal(1)^2+Normal(2)^2+Normal(3)^2);
                    Normal = Normal/n;
                    Center(k,:) = CenterOfCylinder+Normal*R;
                    
                    clear CenterOfCylinder
                    CenterOfCylinder(:,2) = SubTriToMove(:,2);
                    CenterOfCylinder(:,1) = 0;
                    CenterOfCylinder(:,3) = 0;
                    Normal = SubTriToMove-CenterOfCylinder;
                    n=sqrt(Normal(:,1).^2+Normal(:,2).^2+Normal(:,3).^2);
                    Normal = Normal./n;
                    SubTriRet(:,:,k) = CenterOfCylinder+Normal*R;
                end
                end
            end
        end
        
        function [SubTri] = SubTriangles(p, t, Center, itt)
            % Method for creating subtriangles
            TotTri = length(t);
            
            for i=1:TotTri
                % Center of triangle
                M = Center(i,:);
                % Coordinates to triangle points
                r1 = p(t(i,1),:);
                r2 = p(t(i,2),:);
                r3 = p(t(i,3),:);
                % Edges coordinates
                r12=r2-r1;
                r23=r3-r2;
                r13=r3-r1; 
                % Points on triangle edges
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
                if itt
                    [aa1] = ArbitraryAntenna.SubTrianglesIttMet(r1,C1,C5, a1);
                    [aa2] = ArbitraryAntenna.SubTrianglesIttMet(C1,C2,M, a2);        
                    [aa3] = ArbitraryAntenna.SubTrianglesIttMet(C2,r2,C3, a3);
                    [aa4] = ArbitraryAntenna.SubTrianglesIttMet(C2,C3,M, a4);            
                    [aa5] = ArbitraryAntenna.SubTrianglesIttMet(C3,M,C4, a5);
                    [aa6] = ArbitraryAntenna.SubTrianglesIttMet(C1,M,C5, a6);
                    [aa7] = ArbitraryAntenna.SubTrianglesIttMet(C5,M,C6, a7);
                    [aa8] = ArbitraryAntenna.SubTrianglesIttMet(C4,C6,M, a8);
                    [aa9] = ArbitraryAntenna.SubTrianglesIttMet(C4,r3,C6, a9);
                    
                    SubTri(:,:,i)=...
                    [aa1 aa2 aa3 aa4 aa5 aa6 aa7 aa8 aa9];
                else
                    SubTri(:,:,i)=...
                    [a1 a2 a3 a4 a5 a6 a7 a8 a9];
                end
                
            end
        end
        
        function [SubTri] = SubTrianglesIttMet(r1,r2,r3, Center)       
                % Center of triangle
                M = Center;
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
                SubTri(:,:)=...
                    [a1 a2 a3 a4 a5 a6 a7 a8 a9];
        end
         
        function [EdgeList, Basis, BasisLA, BasisCoord] = BasisFunc(p, t, ConnectCell)
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
                    % Minus basis function
                    Basis{BasisNumber,2} = @(r) (r - p(NotEdgePoints(2),:));
                  
                    % Plus basis function
                    BasisCoord{BasisNumber,1} = @(r) (p(NotEdgePoints(1),1) - r);
                    BasisCoord{BasisNumber,2} = @(r) (p(NotEdgePoints(1),2) - r);
                    BasisCoord{BasisNumber,3} = @(r) (p(NotEdgePoints(1),3) - r);
                    % Minus basis function
                    BasisCoord{BasisNumber,4} = @(r) (r - p(NotEdgePoints(2),1));
                    BasisCoord{BasisNumber,5} = @(r) (r - p(NotEdgePoints(2),2));
                    BasisCoord{BasisNumber,6} = @(r) (r - p(NotEdgePoints(2),3));
                  
                    
                    % Fronterm for basis function
                    BasisLA(BasisNumber,1) =  L./(2*AP);
                    BasisLA(BasisNumber,2) =  L;
                    BasisLA(BasisNumber,3) = L./(2*AM);
                
                    BasisNumber = BasisNumber + 1; 
                end
            end
        end
        
        function [RhoP, RhoM, RhoP_, RhoM_] = BasisEvalCenter(t, EdgeList, Basis, Center, SubTri)
            % Looping through triangles
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            for m=1:length(t)
                Plus     =find(PlusTri-m==0);
                Minus    =find(MinusTri-m==0);
                %Looping through the 3 edges for a triangle
                for i=1:length(Plus)
                    TP = Plus(i);
                    BasisP = Basis{TP,1};
                    % Evaluating the basis function in center of the
                    % triangle
                    RhoP(TP,:) = BasisP(Center(m,:));
                    % Evaluating the basis function in centres of
                    % subtriangles
                    SubP = SubTri(:,:,m);
                    RhoP_(:,:,TP) = BasisP(SubP);
                end
                for i=1:length(Minus)
                    TM = Minus(i);
                    BasisM = Basis{TM,2};
                    % Evaluating the basis function in center of the
                    % triangle
                    RhoM(TM,:) = BasisM(Center(m,:));
                    % Evaluating the basis function in centres of
                    % subtriangles
                    SubM = SubTri(:,:,m);
                    RhoM_(:,:,TM) = BasisM(SubM);
                end        
            end
        end
             
        function [PlusTri, MinusTri] = PMTri(t, EdgeList)
            for y=1:length(EdgeList)
                %Finding the points of the plus and minus triangle
                Edge(1:2) = EdgeList(y,1:2);
                FaceEdgeP = [Edge EdgeList(y,3)];
                FaceEdgeP = sort(FaceEdgeP,2);
                FaceEdgeM = [Edge EdgeList(y,4)];
                FaceEdgeM = sort(FaceEdgeM,2);
                % Finding the index of plus and minus triangles for
                % the current basis for the outer
                PlusTri(y) = find(sum(abs(t-FaceEdgeP),2)==0);
                MinusTri(y) = find(sum(abs(t-FaceEdgeM),2)==0); 
            end
        end
        
        function [Ei] = PointSource(w, mu, k, Center, SubTri, sub, PointPos, p) 
            if sub
                dims = size(SubTri);
                SubTri = reshape(SubTri, 3, [], dims(3));
                SubTri = permute(SubTri, [2 1 3]);
                Rx = PointPos(1)-SubTri(:,1,:);
                Ry = PointPos(2)-SubTri(:,2,:);
                Rz = PointPos(3)-SubTri(:,3,:);
            else
                Rx = PointPos(1)-Center(:,1);
                Ry = PointPos(2)-Center(:,2);
                Rz = PointPos(3)-Center(:,3);
            end
            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                                                        
            g = exp(1i.*k.*r)./(4*pi*r);
            G1 = (1+1i./(k*r)-1./(k*r).^2);
            G2 = (1+3i./(k*r)-3./(k*r).^2);

            RR = Rx.*Rx;
            Gxx = (G1-(RR./r.^2).*G2).*g;
            RR = Ry.*Rx;            
            Gxy = (-(RR./r.^2).*G2).*g;
            RR = Rz.*Rx;
            Gxz = (-(RR./r.^2).*G2).*g;

            RR = Ry.*Rx;
            Gyx = (-(RR./r.^2).*G2).*g;
            RR = Ry.*Ry;
            Gyy = (G1-(RR./r.^2).*G2).*g;
            RR = Ry.*Rz;
            Gyz = (-(RR./r.^2).*G2).*g;
                            
            RR = Rz.*Rx;
            Gzx = (-(RR./r.^2).*G2).*g;
            RR = Rz.*Ry;
            Gzy = (-(RR./r.^2).*G2).*g;
            RR = Rz.*Rz;
            Gzz = (G1-(RR./r.^2).*G2).*g;
            
            if sub
                dims = size(g);
                Ei(:,1) = sum(w^2*mu.*(Gxx .* p(1) + Gxy .* p(2) + Gxz .* p(3)))/dims(1);
                Ei(:,2) = sum(w^2*mu.*(Gyx .* p(1) + Gyy .* p(2) + Gyz .* p(3)))/dims(1);
                Ei(:,3) = sum(w^2*mu.*(Gzx .* p(1) + Gzy .* p(2) + Gzz .* p(3)))/dims(1);
            else
                Ei(:,1) = w^2*mu.*(Gxx .* p(1) + Gxy .* p(2) + Gxz .* p(3));
                Ei(:,2) = w^2*mu.*(Gyx .* p(1) + Gyy .* p(2) + Gyz .* p(3));
                Ei(:,3) = w^2*mu.*(Gzx .* p(1) + Gzy .* p(2) + Gzz .* p(3));
                
                Ei(:,1) = w^2*mu.*(Gxx .* p(1) + Gxy .* p(2) + Gxz .* p(3));
                Ei(:,2) = w^2*mu.*(Gyx .* p(1) + Gyy .* p(2) + Gyz .* p(3));
                Ei(:,3) = w^2*mu.*(Gzx .* p(1) + Gzy .* p(2) + Gzz .* p(3));
            end
        end
        
        function [Ei, EdgeV] = VoltageFeed(t, p, Center, FeedPos, v, EdgeList, BasisLA, Yagi, OG)
            
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            Ei = zeros(size(Center));
            EdgeV = 1:length(EdgeList);
            EdgeV(:) = 0;
            
            FeedEdges = [];
            counter = 1;
            if Yagi
               
                PointsOnFeed = find(abs(p(1:OG,2))<=1e-18+FeedPos(2));
            else
                PointsOnFeed = find(abs(p(:,2))<=1e-18+FeedPos(2));
            end
            for i = 1:length(PointsOnFeed)
                    StartPoint = EdgeList(:,1)==PointsOnFeed(i);
                for j = i+1:length(PointsOnFeed)
                    EndPoint = EdgeList(:,2)==PointsOnFeed(j);
                    if sum(StartPoint.*EndPoint)
                        FeedEdges(counter) = find(StartPoint.*EndPoint);
                        counter = counter+1;
                    end
                end
            end
            
            EdgeV(FeedEdges)=BasisLA(FeedEdges,2)*v;
            if ~isempty(FeedEdges)
                Ei(PlusTri(FeedEdges),2) = v/length(FeedEdges);
                Ei(MinusTri(FeedEdges),2) = v/length(FeedEdges);
            end
        end
         
        function [Z,  a, b] = MoM(w, mu, t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Center, k, SubTri, Ei, eps0)
            
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
        
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            
            SubAmount = size(SubTri);
            Quad = SubAmount(1);
            
            for m=1:length(EdgeList)
                    mPdist = sqrt(sum((Center(PlusTri(m),:)-SubTri(:,:,:)).^2,2));
                    mMdist = sqrt(sum((Center(MinusTri(m),:)-SubTri(:,:,:)).^2,2));
                    rhomP = RhoP(m,:);
                    rhomM = RhoM(m,:);
                for n=1:length(EdgeList)
                    rhonP_ = RhoP_(:,:,n);
                    rhonM_ = RhoM_(:,:,n);
                    
                    gmPnP = exp(-1i*k*mPdist(:,:,PlusTri(n)))./mPdist(:,:,PlusTri(n));
                    gmMnP = exp(-1i*k*mMdist(:,:,PlusTri(n)))./mMdist(:,:,PlusTri(n));
                
                    gmPnM = exp(-1i*k*mPdist(:,:,MinusTri(n)))./mPdist(:,:,MinusTri(n));
                    gmMnM = exp(-1i*k*mMdist(:,:,MinusTri(n)))./mMdist(:,:,MinusTri(n));
                        
            AmnP = mu/(4*pi)*(BasisLA(n,2)*sum(rhonP_.*gmPnP/(2*Quad))+BasisLA(n,2)*sum(rhonM_.*gmPnM/(2*Quad)));
            AmnM = mu/(4*pi)*(BasisLA(n,2)*sum(rhonP_.*gmMnP/(2*Quad))+BasisLA(n,2)*sum(rhonM_.*gmMnM/(2*Quad)));
            
            PhiP = -1/(4*pi*1i*w*eps0)*(BasisLA(n,2)*sum(gmPnP)/Quad-BasisLA(n,2)*sum(gmPnM)/Quad);
            PhiM = -1/(4*pi*1i*w*eps0)*(BasisLA(n,2)*sum(gmMnP)/Quad-BasisLA(n,2)*sum(gmMnM)/Quad);
            
            Z(m,n) = BasisLA(m,2)*(1i*w*(dot(AmnP,rhomP)/2+dot(AmnM,rhomM)/2)+PhiM-PhiP);
                end
            end
            
            b = BasisLA(:,2).*(dot(Ei(PlusTri,:),RhoP,2)/2+dot(Ei(MinusTri,:),RhoM,2)/2);
            
            a = Z\b;
        end
    
        function [Z, a, b] = MoMVectorized(w, mu, t, p, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Basis, Center, k, SubTri, Ei,...
                Reflector, GIxx, GIxy, GIxz, GIyx, GIyy, GIyz, GIzx, GIzy, GIzz, eps0)
            % alocating space
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
            
            SubAmount = size(SubTri);
            Quad = SubAmount(1);
         
            EdgesTotal = length(EdgeList);
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);

%             [BasisAnalytic, DistAnalytic] = ArbitraryAntenna.SelfTermInt(t, p, k, Basis, PlusTri, MinusTri, EdgeList);
%             EdgeLevelBasisPlus = BasisAnalytic(PlusTri).';
%             EdgeLevelBasisMinus = BasisAnalytic(MinusTri).';
            
            for m=1:EdgesTotal
                mPdist = sqrt(sum((Center(PlusTri(m),:)-SubTri(:,:,:)).^2,2));
                mMdist = sqrt(sum((Center(MinusTri(m),:)-SubTri(:,:,:)).^2,2));
                rhomP = repmat(RhoP(m,:),length(EdgeList),1);
                rhomM = repmat(RhoM(m,:),length(EdgeList),1);
                
                SamenPmP = find(PlusTri - PlusTri(m) == 0);
                SamenMmP = find(MinusTri - PlusTri(m) == 0);
                SamenMmM = find(MinusTri - MinusTri(m) == 0);
                SamenPmM = find(PlusTri - MinusTri(m) == 0);
                
                gmPnP = exp(1i*k*mPdist(:,:,PlusTri))./mPdist(:,:,PlusTri);
                gmMnP = exp(1i*k*mMdist(:,:,PlusTri))./mMdist(:,:,PlusTri);
                
                gmPnM = exp(1i*k*mPdist(:,:,MinusTri))./mPdist(:,:,MinusTri);
                gmMnM = exp(1i*k*mMdist(:,:,MinusTri))./mMdist(:,:,MinusTri);
                
                Acnst = mu/(4*pi);
                PPA = permute(sum(RhoP_.*gmPnP/(2*Quad)),[3 2 1]);
                MPA = permute(sum(RhoM_.*gmPnM/(2*Quad)),[3 2 1]);
                PMA = permute(sum(RhoP_.*gmMnP/(2*Quad)),[3 2 1]);
                MMA = permute(sum(RhoM_.*gmMnM/(2*Quad)),[3 2 1]);
                
                AmnP = Acnst.*BasisLA(:,2).*(PPA+MPA);
                AmnM = Acnst.*BasisLA(:,2).*(PMA+MMA);
                
                Pcnst = -1/(4*pi*1i*w*eps0); 
                PPPhi = permute(sum(gmPnP),[3 2 1])/(2*Quad);
                PMPhi = permute(sum(gmPnM),[3 2 1])/(2*Quad);
                MPPhi = permute(sum(gmMnP),[3 2 1])/(2*Quad);
                MMPhi = permute(sum(gmMnM),[3 2 1])/(2*Quad);
                 
%                 PPPhi(SamenPmP) = DistAnalytic(PlusTri(m));
%                 PMPhi(SamenMmP) = DistAnalytic(PlusTri(m));
%                 MMPhi(SamenMmM) = DistAnalytic(MinusTri(m));
%                 MPPhi(SamenPmM) = DistAnalytic(MinusTri(m));
                 
                PhiP = Pcnst*BasisLA(:,2).*(PPPhi-PMPhi);
                PhiM = Pcnst*BasisLA(:,2).*(MPPhi-MMPhi);
                
                if Reflector
                    GIx = [GIxx(:,PlusTri(m),:) GIxy(:,PlusTri(m),:) GIxz(:,PlusTri(m),:)];
                    GIy = [GIyx(:,PlusTri(m),:) GIyy(:,PlusTri(m),:) GIyz(:,PlusTri(m),:)];
                    GIz = [GIzx(:,PlusTri(m),:) GIzy(:,PlusTri(m),:) GIzz(:,PlusTri(m),:)];
                    GImP = GIx + GIy + GIz;
                    GImP = [GIxx(:,PlusTri(m),:) GIyy(:,PlusTri(m),:) GIzz(:,PlusTri(m),:)];
                    GIx = [GIxx(:,MinusTri(m),:) GIxy(:,MinusTri(m),:) GIxz(:,MinusTri(m),:)];
                    GIy = [GIyx(:,MinusTri(m),:) GIyy(:,MinusTri(m),:) GIyz(:,MinusTri(m),:)];
                    GIz = [GIzx(:,MinusTri(m),:) GIzy(:,MinusTri(m),:) GIzz(:,MinusTri(m),:)];
                    GImM = GIx+GIy+GIz;
                    GImM = [GIxx(:,MinusTri(m),:) GIyy(:,MinusTri(m),:) GIzz(:,MinusTri(m),:)];
             
                    GImPnP = GImP(:,:,PlusTri);
                    GImMnP = GImM(:,:,MinusTri);
                    
                    GImPnM = GImP(:,:,PlusTri);
                    GImMnM = GImM(:,:,MinusTri);
                    
                    Z(m,:) = (1i*w*mu/(4*pi).*BasisLA(:,2).*BasisLA(m,2).*(dot(permute(sum(RhoP_.*GImPnP/(2*Quad)),[3 2 1])...
                        +permute(sum(RhoM_.*GImPnM/(2*Quad)),[3 2 1]),rhomP,2)+...
                        dot(permute(sum(RhoP_.*GImMnP/(2*Quad)),[3 2 1])...
                        +permute(sum(RhoM_.*GImMnM/(2*Quad)),[3 2 1]),rhomM,2))).';   
                end
                
            PlusDotProd = dot(AmnP,rhomP,2);
            MinusDotProd = dot(AmnM,rhomM,2);
%                 PlusDotProd(SamenPmP)  = Acnst.*BasisLA(SamenPmP,2).*EdgeLevelBasisPlus(SamenPmP);%PlusTri(m)
%                 PlusDotProd(SamenMmP)  = -Acnst.*BasisLA(SamenMmP,2).*EdgeLevelBasisPlus(SamenMmP);% ;PlusTri(m)
%                 MinusDotProd(SamenMmM) = -Acnst.*BasisLA(SamenMmM,2).*EdgeLevelBasisMinus(SamenMmM);%;MinusTri(m)
%                 MinusDotProd(SamenPmM) = Acnst.*BasisLA(SamenPmM,2).*EdgeLevelBasisMinus(SamenPmM);%;MinusTri(m)
             
%                 PlusDotProd(SamenPmP)  = -Acnst/BasisLA(SamenPmP,2).*BasisAnalytic(PlusTri(m))/BasisLA(m,2);
%                 PlusDotProd(SamenMmP)  = -Acnst/BasisLA(SamenMmP,2).*BasisAnalytic(PlusTri(m))/BasisLA(m,2);
%                 MinusDotProd(SamenMmM) = Acnst/BasisLA(SamenMmM,2).*BasisAnalytic(MinusTri(m))/BasisLA(m,2);
%                 MinusDotProd(SamenPmM) = Acnst/BasisLA(SamenPmM,2).*BasisAnalytic(MinusTri(m))/BasisLA(m,2);
                
                Z(m,:) = BasisLA(m,2).*(1i*w*(PlusDotProd/2+MinusDotProd/2)+PhiM-PhiP).'+Z(m,:);
            end
            b = BasisLA(:,2).*(dot(Ei(PlusTri,:),RhoP,2)/2+dot(Ei(MinusTri,:),RhoM,2)/2);

            %System solution
            a=Z\b;
        end
   
        function [BasisAnalytic, DistAnalytic] = SelfTermInt(t, p, k, Basis, PlusTri, MinusTri, EdgeList)
            BasisAnalytic = 1:length(t); BasisAnalytic(:)=0;
            DistAnalytic = 1:length(t); DistAnalytic(:)=0;
            for i=1:length(t)
                    v1 = p(t(i,1),:);
                    v2 = p(t(i,2),:);
                    v3 = p(t(i,3),:);
                    
                    a = dot((v3-v1),(v3-v1));
                    b = dot((v3-v1),(v3-v2));
                    c = dot((v3-v2),(v3-v2));
                    
                    l1 = sqrt(c);
                    l2 = sqrt(a);
                    l3 = sqrt(a-2*b+c);
                    
                    ln1 = log(((l1+l2)^2-l3^2)/(l2^2-(l3-l1)^2));
                    ln2 = log(((l2+l3)^2-l1^2)/(l3^2-(l1-l2)^2));
                    ln3 = log(((l3+l1)^2-l2^2)/(l1^2-(l2-l3)^2));
              
                    I2 = 1/(3*l1)*ln1 + 1/(3*l2)*ln2 + 1/(3*l3)*ln3;
                    
                    DistAnalytic(i) = 1i*k+I2;
            end
            for j=1:2
            for i=1:length(t)
                A = t(i,1);
                B = t(i,2);
                C = t(i,3);
                if j==1
                    v1 = p(A,:);
                    v2 = p(B,:);
                    v3 = p(C,:);
                elseif j ==2
                    v1 = p(A,:);
                    v2 = p(C,:);
                    v3 = p(B,:);
                end
                        a = dot((v3-v1),(v3-v1));
                        b = dot((v3-v1),(v3-v2));
                        c = dot((v3-v2),(v3-v2));
                    
                        l1 = sqrt(c);
                        l2 = sqrt(a);
                        l3 = sqrt(a-2*b+c);
                    
                        ln1 = log(((l1+l2)^2-l3^2)/(l2^2-(l3-l1)^2));
                        ln2 = log(((l2+l3)^2-l1^2)/(l3^2-(l1-l2)^2));
                        ln3 = log(((l3+l1)^2-l2^2)/(l1^2-(l2-l3)^2));
              
                        L1L1 = 1/(20*l1)*ln1 + (l1^2+5*l2^2-l3^2)/(120*l2^3)*ln2...
                        +(l1^2-l2^2+5*l3^2)/(120*l3^3)*ln3+(l3-l1)/(60*l2^2)+(l2-l1)/60*l3^2;
    
                        L1L2 = (3*l1^2+l2^2-l3^2)/(80*l1^2)*ln1+(l1^2+3*l2^2-l3^2)/(80*l2^3)*ln2...
                        +1/(40*l3)*ln3+(l3-l2)/(40*l1^2)+(l3-l1)/(40*l2^2);
                
                        LI = 1/(8*l2)*ln1 + (l1^2+5*l2^2-l3^2)/(48*l2^3)*ln2...
                        +(l1^2-l2^2+5*l3^2)/(48*l3^3)*ln3 + (l3-l1)/(24*l2^2)+(l2-l1)/(24*l3^2);
                for n = 1:3
                    vn = p(t(i,n),:);
                    for m=1:3
                        vm = p(t(i,m),:);
                        
%                         a11 = dot(v1,v1); a12 = dot(v1,v2); a13 = dot(v1,v3);
%                         a22 = dot(v2,v2); a23 = dot(v2,v3); a33 = dot(v3,v3);
%                         a1n = dot(v1,vn); a1m = dot(v1,vm); a2n = dot(v2,vn);
%                         a2m = dot(v2,vm); a3m = dot(v3,vm); a3n = dot(v3,vn);
%                         anm = dot(vn,vm);

                        a11 = dot(p(A,:),p(A,:)); a12 = dot(p(A,:),p(B,:)); a13 = dot(p(A,:),p(C,:));
                        a22 = dot(p(B,:),p(B,:)); a23 = dot(p(B,:),p(C,:)); a33 = dot(p(C,:),p(C,:));
                        a1n = dot(p(A,:),vn); a1m = dot(p(A,:),vm); a2n = dot(p(B,:),vn);
                        a2m = dot(p(B,:),vm); a3m = dot(p(C,:),vm); a3n = dot(p(C,:),vn);
                        anm = dot(vn,vm);
                        
                        if j==1
                        BasisAnalytic(i) =  (DistAnalytic(i))*(...
                                                L1L1*(a11-2*a12+a22)+...
                                                +L1L2*(a11-a13-a12+a23)+...
                                                +LI*(-a11+a1n+a12-a2n)+...
                                                LI*(-a11+a1m+a12-a2m)+...
                                                a11-a1n-a1m+anm)...
                                                + BasisAnalytic(i); 
                        elseif j==2
                        BasisAnalytic(i) =  (DistAnalytic(i))*(...
                                                L1L1*(a11-2*a13+a33)...
                                                +L1L2*(a11-a12-a13+a23)...
                                                +LI*(-a11+a1n+a13-a3n)+...
                                                LI*(-a11+a1m+a13-a3m))...
                                                + BasisAnalytic(i); 
                        end
                    end
                end
            end
            end
%             for j = 1:2
%                 for i=1:length(EdgeList)
%                     First = EdgeList(i,1);
%                     Second = EdgeList(i,2);
%                 
%                     v1 = p(EdgeList(i,1),:);
%                     v2 = p(EdgeList(i,2),:);
%                     if j ==1
%                         Third = EdgeList(i,3);
%                         v3 = p(EdgeList(i,3),:);
%                     elseif j==2
%                         Third = EdgeList(i,4);
%                         v3 = p(EdgeList(i,4),:);
%                     end
%                     TriPoints = [First, Second, Third];
%                     for l=1:3
%                         vm = p(TriPoints(l),:);
%                         for o=1:3
%                         vn = p(TriPoints(o),:);
%                     %%
%                     a = dot((v3-v1),(v3-v1));
%                     b = dot((v3-v1),(v3-v2));
%                     c = dot((v3-v2),(v3-v2));
%                     
%                     l1 = sqrt(c);
%                     l2 = sqrt(a);
%                     l3 = sqrt(a-2*b+c);
%                     
%                     ln1 = log(((l1+l2)^2-l3^2)/(l2^2-(l3-l1)^2));
%                     ln2 = log(((l2+l3)^2-l1^2)/(l3^2-(l1-l2)^2));
%                     ln3 = log(((l3+l1)^2-l2^2)/(l1^2-(l2-l3)^2));
%               
%                     I11 = 1/(20*l1)*ln1 + (l1^2+5*l2^2-l3^2)/(120*l2^3)*ln2...
%                     +(l1^2-l2^2+5*l3^2)/(120*l3^3)*ln3+(l3-l1)/(60*l2^2)+(l2-l1)/60*l3^2;
%                      
%                     I12 = (3*l1^2+l2^2-l3^2)/(80*l1^2)*ln1+(l1^2+3*l2^2-l3^2)/(80*l2^3)*ln2...
%                     +1/(40*l3)*ln3+(l3-l2)/(40*l1^2)+(l3-l1)/(40*l2^2);
%                 
%                     I = 1/(8*l2)*ln1 + (l1^2+5*l2^2-l3^2)/(48*l2^3)*ln2...
%                     +(l1^2-l2^2+5*l3^2)/(48*l3^3)*ln3 + (l3-l1)/(24*l2^2)+(l2-l1)/(24*l3^2);
%                     %%
%                     a11 = dot(v1,v1); a12 = dot(v1,v2); a13 = dot(v1,v3);
%                     a22 = dot(v2,v2); a23 = dot(v2,v3); a33 = dot(v3,v3);
%                     a1n = dot(v1,vn); a1m = dot(v1,vm); a2n = dot(v2,vn);
%                     a2m = dot(v2,vm); a3m = dot(v3,vm); a3n = dot(v3,vn);
%                     anm = dot(vn,vm);
%                     
%                     TriPoint = sort(TriPoints);
%                     Triangle1 = t(:,1)==TriPoint(1);
%                     Triangle2 = t(:,2)==TriPoint(2);
%                     Triangle3 = t(:,3)==TriPoint(3);
%                     Triangle = find(Triangle1.*Triangle2.*Triangle3);
%                     
%                     BasisAnalytic(i,j) =  (DistAnalytic(Triangle))*(...
%                                                 I11*(a11-2*a12+a22)+...
%                                                 I11*(a11-2*a13+a33)...
%                                                 +I12*(a11-a13-a12+a23)+...
%                                                 I12*(a11-a12-a13+a23)...
%                                                 +I*(-a11+a1n+a12-a2n)+...
%                                                 I*(-a11+a1n+a13-a3n)+...
%                                                 I*(-a11+a1m+a12-a2m)+...
%                                                 I*(-a11+a1m+a13-a3m)+...
%                                                 a11-a1n-a1m+anm)...
%                                                 + BasisAnalytic(i,j);
%                         end
%                     end
%                 end
%             end
%             BasisAnalytic = [];
%             DistAnalytic = [];
%             for i=1:length(t)
%                 Plus = find( PlusTri - i == 0);
%                 Minus = find( MinusTri - i == 0);
%                 for j=1:length(Plus)
%                     for self=1:3
% %                     Integration for selfterms
%                     v1 = p(t(i,1),self); v2 = p(t(i,2),self); v3 = p(t(i,3),self);
%                     rp = @(al, be) (1-al-be).*v1+al.*v2+be.*v3;
%                     v1 = p(t(i,1),self); v2 = p(t(i,2),self); v3 = p(t(i,3),self);
%                     rm = @(al, be) (1-al-be).*v1+al.*v2+be.*v3;
%                     
%                     RhoPlus = Basis{PlusTri(i), self};
%                     RhoPlus = @(be,al) RhoPlus(rp(al, be));
%                     RhoMinus = Basis{MinusTri, self+3};
%                     RhoMinus = @(be,al) RhoMinus(rm(al, be));
%                                     
%                     gp = @(be,al) exp(1i*k*rp(al, be))./(rp(al, be));
%                     gm = @(be,al) exp(1i*k*rm(al, be))./(rm(al, be));
%                     bem = @(al) 1-al;
%                     
%                     funcPP = @(be,al) RhoPlus(be,al).*gp(be,al);
%                     funcMP = @(be,al) RhoMinus(be,al).*gp(be,al);
%                     funcPM = @(be,al) RhoPlus(be,al).*gm(be,al);
%                     funcMM = @(be,al) RhoMinus(be,al).*gm(be,al);
%                     
%                     BasisAnalytic(i,self) =(integral2( funcPP , 0, 1, 0, bem)...
%                         +integral2( funcMP , 0, 1, 0, bem));
%                     BasisAnalytic(i,self) = (integral2( funcPM , 0, 1, 0, bem)...
%                         +integral2( funcMM , 0, 1, 0, bem));
%                 
% %                     Integration for selfterms
%                     DistAnalytic(i) = (integral2( gp, 0, 1, 0, bem)...
%                         +integral2( gp, 0, 1, 0, bem));
%                     DistAnalytic(i) = (integral2( gm , 0, 1, 0, bem)...
%                         +integral2( gm , 0, 1, 0, bem));   
%                     end
%                 end
%                 for j=1:length(Minus) 
%                     for self=1:3
% %                     Integration for selfterms
%                     v1 = p(t(i,1),self); v2 = p(t(i,2),self); v3 = p(t(i,3),self);
%                     rp = @(al, be) (1-al-be).*v1+al.*v2+be.*v3;
%                     v1 = p(t(i,1),self); v2 = p(t(i,2),self); v3 = p(t(i,3),self);
%                     rm = @(al, be) (1-al-be).*v1+al.*v2+be.*v3;
%                     
%                     RhoPlus = Basis{PlusTri, self};
%                     RhoPlus = @(be,al) RhoPlus(rp(al, be));
%                     RhoMinus = Basis{MinusTri, self+3};
%                     RhoMinus = @(be,al) RhoMinus(rm(al, be));
%                                     
%                     gp = @(be,al) exp(1i*k*rp(al, be))./(rp(al, be));
%                     gm = @(be,al) exp(1i*k*rm(al, be))./(rm(al, be));
%                     bem = @(al) 1-al;
%                     
%                     funcPP = @(be,al) RhoPlus(be,al).*gp(be,al);
%                     funcMP = @(be,al) RhoMinus(be,al).*gp(be,al);
%                     funcPM = @(be,al) RhoPlus(be,al).*gm(be,al);
%                     funcMM = @(be,al) RhoMinus(be,al).*gm(be,al);
%                     
%                     BasisAnalytic(i,self) =(integral2( funcPP , 0, 1, 0, bem)...
%                         +integral2( funcMP , 0, 1, 0, bem));
%                     BasisAnalytic(i,self) = (integral2( funcPM , 0, 1, 0, bem)...
%                         +integral2( funcMM , 0, 1, 0, bem));
%                 
% %                     Integration for selfterms
%                     DistAnalytic(i) = (integral2( gp, 0, 1, 0, bem)...
%                         +integral2( gp, 0, 1, 0, bem));
%                     DistAnalytic(i) = (integral2( gm , 0, 1, 0, bem)...
%                         +integral2( gm , 0, 1, 0, bem));   
%                     end
%                 end
%             end
        end
           
        function [Jface] = CurrentCalc(t, EdgeList, a, BasisLA, RhoP, RhoM)
            Jface = zeros(size(t));
            
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            for n=1:length(t)
                Plus     =  find(PlusTri-n==0);
                Minus    =  find(MinusTri-n==0);
                
                for y=1:length(Plus)
                    Jface(n,:) = a(Plus(y))*BasisLA(Plus(y),1)*RhoP(Plus(y),:) + Jface(n,:);    
                end
                for y=1:length(Minus)
                    Jface(n,:) = a(Minus(y))*BasisLA(Minus(y),3)*RhoM(Minus(y),:) + Jface(n,:);    
                end    
            end
        end
        
        function [Exy, Exz, Ezy, xrange, yrange, zrange, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = EField(Center, w, mu, k0, J,...
                xmin, xmax,  ymin, ymax,zmin, zmax, steps, Area, Reflect, xsurf, n, lambda)
       
            kR = 2*pi/(lambda*n);
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            zrange = linspace(zmin, zmax, steps);
            
            Exy = zeros(steps, steps); 
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
                
               for i=1:length(Center)
                        if j == 1
                            %xy
                            Rx = repmat(rx(i,:)',1,steps);
                            Ry = repmat(ry(i,:),steps,1);
                            Rz = repmat(rz(i),steps,steps);
                            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                            surfside = find(rx(i,:)>=xsurf);
                            k = zeros(steps,steps);
                            k(:,:) = k0; 
                            if Reflect
                                k(surfside,:) = kR;
                            end
                            
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k.*r)-1./(k.*r).^2);
                            G2 = (1+3i./(k.*r)-3./(k.*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rx;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rx;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Exyx = Exyx + 1i.*w.*mu.*(Gxx .* J(i,1) + Gxy .* J(i,2) + Gxz .* J(i,3))*Area(i);
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Exyy = Exyy + 1i.*w.*mu.*(Gyx .* J(i,1) + Gyy .* J(i,2) + Gyz .* J(i,3))*Area(i);
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Exyz = Exyz + 1i.*w.*mu.*(Gzx .* J(i,1) + Gzy .* J(i,2) + Gzz .* J(i,3))*Area(i);
                        elseif j==2
                            %xz
                            Rx = repmat(rx(i,:)',1,steps);
                            Ry = repmat(ry(i,:),steps,steps);
                            Rz = repmat(rz(i,:),steps,1);
                            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                            surfside = find(rx(i,:)>=xsurf);
                            k = zeros(steps,steps);
                            k(:,:) = k0; 
                            if Reflect
                                k(surfside,:) = kR;
                            end
                            
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k.*r)-1./(k.*r).^2);
                            G2 = (1+3i./(k.*r)-3./(k.*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Rx.*Ry;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rx.*Rz;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Exzx = Exzx + 1i.*w.*mu.*(Gxx .* J(i,1) + Gxy .* J(i,2) + Gxz .* J(i,3))*Area(i);
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Exzy = Exzy + 1i.*w.*mu.*(Gyx .* J(i,1) + Gyy .* J(i,2) + Gyz .* J(i,3))*Area(i);
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Exzz = Exzz + 1i.*w.*mu.*(Gzx .* J(i,1) + Gzy .* J(i,2) + Gzz .* J(i,3))*Area(i);
                        else
                            %yz
                            Rx = repmat(rx(i,:),steps,steps);
                            Ry = repmat(ry(i,:),steps,1);
                            Rz = repmat(rz(i,:)',1,steps);
                            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                            surfside = find(rx(i,:)>=xsurf);
                            k = zeros(steps,steps);
                            k(:,:) = k0; 
                            if Reflect
                                k(surfside,:) = kR;
                            end
                            
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k.*r)-1./(k.*r).^2);
                            G2 = (1+3i./(k.*r)-3./(k.*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Rx.*Ry;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rx.*Rz;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Eyzx = Eyzx + 1i.*w.*mu.*(Gxx .* J(i,1) + Gxy .* J(i,2) + Gxz .* J(i,3))*Area(i);
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Eyzy = Eyzy + 1i.*w.*mu.*(Gyx .* J(i,1) + Gyy .* J(i,2) + Gyz .* J(i,3))*Area(i);
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Eyzz = Eyzz + 1i.*w.*mu.*(Gzx .* J(i,1) + Gzy .* J(i,2) + Gzz .* J(i,3))*Area(i);
                        end  
                end 
            end
            Exy = sqrt(Exyx.^2+Exyy.^2+Exyz.^2);
            Exz = sqrt(Exzx.^2+Exzy.^2+Exzz.^2);
            Ezy = sqrt(Eyzx.^2+Eyzy.^2+Eyzz.^2);
        end
             
        function [Exy, Exz, Ezy, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = PointSourceEmmision(Center, k, p,...
                xmin, xmax, zmin, zmax, ymin, ymax, steps)
       
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            zrange = linspace(zmin, zmax, steps);
            
            Exy = zeros(steps, steps); 
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
                
                if j == 1
                            %xy
                            Rx = repmat(rx(:),1,steps);
                            Ry = repmat(ry(:)',steps,1);
                            Rz = repmat(rz,steps,steps);
                            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                                                        
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k*r)-1./(k*r).^2);
                            G2 = (1+3i./(k*r)-3./(k*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rx;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rx;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Exyx = Exyx + (Gxx .* p(1) + Gxy .* p(2) + Gxz .* p(3));
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Exyy = Exyy + (Gyx .* p(1) + Gyy .* p(2) + Gyz .* p(3));
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Exyz = Exyz + (Gzx .* p(1) + Gzy .* p(2) + Gzz .* p(3));
                        elseif j==2
                            %xz
                            Rx = repmat(rx(:)',steps,1);
                            Ry = repmat(ry(:),steps,steps);
                            Rz = repmat(rz(:),1,steps);
                            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                            
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k*r)-1./(k*r).^2);
                            G2 = (1+3i./(k*r)-3./(k*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Rx.*Ry;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rx.*Rz;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Exzx = Exzx + (Gxx .* p(1) + Gxy .* p(2) + Gxz .* p(3));
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Exzy = Exzy + (Gyx .* p(1) + Gyy .* p(2) + Gyz .* p(3));
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Exzz = Exzz + (Gzx .* p(1) + Gzy .* p(2) + Gzz .* p(3));
                        else
                            %yz
                            Rx = repmat(rx(:),steps,steps);
                            Ry = repmat(ry(:)',steps,1);
                            Rz = repmat(rz(:),1,steps);
                            r = sqrt(Rx.^2+Ry.^2+Rz.^2);
                                
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k*r)-1./(k*r).^2);
                            G2 = (1+3i./(k*r)-3./(k*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Rx.*Ry;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rx.*Rz;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Eyzx = Eyzx + (Gxx .* p(1) + Gxy .* p(2) + Gxz .* p(3));
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Eyzy = Eyzy + (Gyx .* p(1) + Gyy .* p(2) + Gyz .* p(3));
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Eyzz = Eyzz + (Gzx .* p(1) + Gzy .* p(2) + Gzz .* p(3));
                end    
            end
            Exy = sqrt(Exyx.^2+Exyy.^2+Exyz.^2);
            Exz = sqrt(Exzx.^2+Exzy.^2+Exzz.^2);
            Ezy = sqrt(Eyzx.^2+Eyzy.^2+Eyzz.^2);
        end
   
        function [Esc, EscPhi, EscTheta] = AngularFarField(w, mu, k, r, Center, J, steps, Area)
            phi = pi/2;% linspace(-pi, 2*pi, steps)';
            theta = linspace(-pi, 2*pi, steps)';
            xH = [1, 0, 0];
            yH = [0, 1, 0];
            zH = [0, 0, 1];
            phiH = -xH.*sin(phi)+yH.*cos(phi);
            thetaH = (xH.*cos(phi)+yH.*sin(phi)).*cos(theta)-zH.*sin(theta);
            rHat = (xH.*cos(phi)+yH.*sin(phi)).*sin(theta)+zH.*cos(theta);
            
            EscTheta = 1:steps; EscTheta(:) =0;
            EscPhi = 1:steps; EscPhi(:) = 0;
           for i=1:length(Center)
                           
                DirectGreens = J(i,:).*(exp(1i*k*r)/(4*pi*r))...
                .*exp(-1i*k*dot(rHat,repmat(Center(i,:),length(rHat),1),2));
            
                DirectGreensTheta = dot(DirectGreens,thetaH,2);
                
                DirectGreensPhi = dot(DirectGreens,repmat(phiH,steps,1),2);
                
                EscTheta = -1i*w*mu*Area(i).*DirectGreensTheta.' + EscTheta;
            
                EscPhi = -1i*w*mu*Area(i).*DirectGreensPhi.' + EscPhi;
            end
            
            Esc = EscPhi+EscTheta;
            
            figure(4)
            plot(theta, 1/2*abs(Esc).^2*r^2)
            xlabel('Theta')
            ylabel('E')
        end
        
        function [Esc] = AngularFarFieldSurf(w, mu, k, r, Center, J, steps, Area, dist, eps2, eps1, n)
            
                phi = pi/2;
                theta = linspace(-pi, 2*pi, steps)';
            
                under = theta>=pi/2;
                over = theta<=3/2*pi;
                UseTrans = logical(under.*over);
            %%
            eps2=eps2*eps1;
            k1 = k;
            k2 = k*n;
            kz1 = k1*cos(theta);
            krho = k1*sin(theta);
            krho2 = k2*sin(theta);
            kz2 = k2*cos(theta);
            
            refS = (kz1-kz2)./(kz1+kz2);
            refP = (eps2.*kz1-eps1.*kz2)./(eps2.*kz1+eps1.*kz2);
            traS = 1+refS;
            traP = 1+refP;
            %%
            xH = [1, 0, 0];
            yH = [0, 1, 0];
            zH = [0, 0, 1];
            phiH = -xH.*sin(phi)+yH.*cos(phi);
            
            thetaHD = (xH.*cos(phi)+yH.*sin(phi)).*cos(theta)-zH.*sin(theta);
            
            thetaHT = xH.*cos(phi).*(cos(theta).*kz1.*kz2/k1^2-sin(theta).*kz1.*krho2/k1^2)...
                +yH.*sin(phi).*(cos(theta).*kz1.*kz2/k1^2-sin(theta).*kz1.*krho2/k1^2)...
                +zH.*(-sin(theta).*krho2.^2/k1^2+cos(theta).*kz2.*krho2/k1^2);
            
            rHatT = xH.*cos(phi).*(sin(theta).*kz1.*kz2/k1^2 ...
            +cos(theta).*kz1.*krho2/k1^2)...
            +yH.*sin(phi).*(sin(theta).*kz1.*kz2/k1^2+cos(theta).*kz1.*krho2/k1^2) ...
            +zH.*(cos(theta).*krho2.^2/k1^2+sin(theta).*kz2.*krho2/k1^2);
            
            thetaHI = zH.*krho/k1+xH.*cos(phi).*kz1/k1+yH.*sin(phi).*kz1/k1;
        
            rHat = (xH.*cos(phi)+yH.*sin(phi)).*sin(theta)+zH.*cos(theta);
            
            rhoH = xH.*cos(theta)+yH.*sin(theta);
    %%        
            EscThetaD = 1:steps; EscThetaD(:) =0;
            EscPhiD = 1:steps; EscPhiD(:) = 0;

            EscThetaI = 1:steps; EscThetaI(:) =0;
            EscPhiI = 1:steps; EscPhiI(:) = 0;
           
            EscThetaT = 1:steps; EscThetaT(:) =0;
            EscPhiT = 1:steps; EscPhiT(:) = 0;
            EscrT = 1:steps; EscrT(:) = 0;
            
            EscD = EscThetaD; EscI = EscThetaD;
            EscT = EscThetaD; Esc = EscThetaD;
            %%
           for i=1:length(Center)
               z = Center(i,3);
               
                DirectGreens = (exp(1i*k*r)/(4*pi*r))...
                .*exp(-1i*k*dot(rHat,repmat(Center(i,:),length(rHat),1),2));
            DirectGreensTheta = DirectGreens.*thetaHD; 
            DirectGreensPhi = DirectGreens.*phiH;
            
                DirectGreensTheta = dot(DirectGreensTheta, repmat(J(i,:),steps,1),2);
                DirectGreensPhi = dot(DirectGreensPhi, repmat(J(i,:),steps,1),2);
              
%%
                    IDGreensBase = exp(1i*k1*r)/(4*pi*r)...
                    .*exp(-1i*krho.*dot(rhoH,repmat(Center(i,:),length(rhoH),1),2)).*exp(1i*kz1.*z);
                
                IDGreensBaseTheta = IDGreensBase.*thetaHI;
                IDGreensBasePhi = IDGreensBase.*phiH;
                
                    IndirectGreensPhi = dot(IDGreensBasePhi,repmat(J(i,:),steps,1),2).*refS;
                
                    IndirectGreensTheta =  -refP.*dot(IDGreensBaseTheta,repmat(J(i,:),steps,1),2);
                %%
                TransGreensBase = kz2./kz1.*exp(1i*k2*r)/(4*pi*r).*exp(1i*kz1.*z).*...
                    exp(-1i*krho2.*dot(rhoH,repmat(Center(i,:),length(rhoH),1),2));
                TransGreensPhi = TransGreensBase .*phiH;
                TransGreensTheta = TransGreensBase .*thetaHT;
                TransGreensR = TransGreensBase .*rHatT;
                
                TransmitGreensPhi =  dot(TransGreensPhi,repmat(J(i,:),steps,1),2).*traS;
                
                TransmitGreensTheta = dot(TransGreensTheta,repmat(J(i,:),steps,1),2).*traP.*eps1/eps2;
                
                TransmitGreensR = dot(TransGreensR,repmat(J(i,:),steps,1),2).*traP.*eps1/eps2;
               %% 
                EscThetaD = -1i*w*mu*Area(i).*DirectGreensTheta.' + EscThetaD;
                EscPhiD = -1i*w*mu*Area(i).*DirectGreensPhi.' + EscPhiD;
                
                EscThetaI = -1i*w*mu*Area(i).*IndirectGreensTheta.' + EscThetaI;
                EscPhiI = -1i*w*mu*Area(i).*IndirectGreensPhi.' + EscPhiI;
                
                EscThetaT = -1i*w*mu*Area(i).*TransmitGreensTheta.' + EscThetaT;
                EscPhiT = -1i*w*mu*Area(i).*TransmitGreensPhi.' + EscPhiT;
                EscrT = -1i*w*mu*Area(i).*TransmitGreensR.' + EscrT;
                
           end
                EscD = EscPhiD+EscThetaD;
                EscI = EscPhiI+EscThetaI;
%                 Esc = EscD+EscI;
                Esc = abs(EscPhiD+EscPhiI).^2+abs(EscThetaD+EscThetaI).^2;
                EscT = abs(EscPhiT).^2+abs(EscThetaT).^2+abs(EscrT).^2;
        
                Esc(UseTrans) = EscT(UseTrans);
            
            figure(6)
            plot(theta, 1/2*Esc*r^2)
            xlabel('Theta')
            ylabel('E')
        end
        
        
        function [GIxx, GIxy, GIxz, GIyx, GIyy, GIyz, GIzx, GIzy, GIzz] = IDGreens(k0, dist, ant_length, ant_width, dx, Nz, lambda, n, epsL2, eps1, Center, SubTri)
                
            SubAmount = size(SubTri);
            %% Soender
            k1i2=(2*pi/(lambda*n))^2;
            epsL2=epsL2*eps1;
            
            kz1f=@(krho) sqrt(k0^2*eps1-krho.*krho);
            kzL2f=@(krho) sqrt(k0^2*epsL2-krho.*krho);
    
            rpf=@(krho) (epsL2*kz1f(krho)-eps1*kzL2f(krho))./(epsL2*kz1f(krho)+eps1*kzL2f(krho));
            rsf=@(krho) (kz1f(krho)-kzL2f(krho))./(kz1f(krho)+kzL2f(krho));
            
            ellipse_length = k0*5;
            ellipse_height = k0*0.2;
            
            krhof=@(alpha) (1+cos(alpha))*ellipse_length/2+1i*ellipse_height*sin(alpha);
            dkrho_dalphaf=@(alpha) -sin(alpha)*ellipse_length/2+1i*ellipse_height*cos(alpha);
            
            rhomin = 0;
            rhomax = ant_length; % maybe half is enough due to symetri?
            
            rhotabv = linspace(rhomin,rhomax,Nz*3);
            Nrho=length(rhotabv);

            GItabzz=zeros(Nrho,Nz); GItabzr=GItabzz;
            GItabrr=GItabzz; GItabpp=GItabzz;
            tab_z=GItabzz; tab_r=GItabzz;
            
            for jrho=1:Nrho
                rho=rhotabv(jrho);
                for jz=1:Nz
                    z=jz*dx;
                    dJmf = @(krho) -besselj(1,krho*rho)./(krho*rho);
                    if jrho==1
                        dJmf = @(krho) -0.5;
                    end
                    Jmmf = @(krho) -0.5*(besselj(0,krho*rho)-besselj(2,krho*rho));
                    
                    Gizzf = @(krho) 1i/(4*pi)*rpf(krho).*besselj(0,krho*rho).*(krho.^2)/k1i2.*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
                    Gizrf = @(krho) 1/(4*pi)*rpf(krho).*krho.*kz1f(krho)/k1i2.*(-besselj(1,krho*rho)).*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
                    Gippf = @(krho) 1i/(4*pi)*(rpf(krho).*dJmf(krho).*(kz1f(krho).^2)/k1i2-rsf(krho).*Jmmf(krho)).*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
                    Girrf = @(krho) 1i/(4*pi)*(rpf(krho).*Jmmf(krho).*(kz1f(krho).^2)/k1i2-rsf(krho).*dJmf(krho)).*exp(1i*kz1f(krho)*z).*krho./kz1f(krho);
  
                    % Integral going into the complex plane to avoid poles
                    Gzz = integral(@(alpha) Gizzf(krhof(alpha)).*dkrho_dalphaf(alpha), -pi, 0, 'ArrayValued', 1);
                    Gzr = integral(@(alpha) Gizrf(krhof(alpha)).*dkrho_dalphaf(alpha), -pi, 0, 'ArrayValued', 1);
                    Gpp = integral(@(alpha) Gippf(krhof(alpha)).*dkrho_dalphaf(alpha), -pi, 0, 'ArrayValued', 1);
                    Grr = integral(@(alpha) Girrf(krhof(alpha)).*dkrho_dalphaf(alpha), -pi, 0, 'ArrayValued', 1);
    
                    % The rest of the integrals
                    Gzz = Gzz + integral(@(krho) Gizzf(krho), ellipse_length, inf, 'ArrayValued', 1);
                    Gzr = Gzr + integral(@(krho) Gizrf(krho), ellipse_length, inf, 'ArrayValued', 1);
                    Gpp = Gpp + integral(@(krho) Gippf(krho), ellipse_length, inf, 'ArrayValued', 1);
                    Grr = Grr + integral(@(krho) Girrf(krho), ellipse_length, inf, 'ArrayValued', 1);                

                    GItabzz(jrho,jz) = Gzz;
                    GItabzr(jrho,jz) = Gzr;
                    GItabpp(jrho,jz) = Gpp;
                    GItabrr(jrho,jz) = Grr;
                    tab_z(jrho,jz) = z;
                    tab_r(jrho,jz) = rho;
                end
            end

            %% Conversion to cartesian components
            TempCenters = permute(Center, [3 2 1]);
            TempCenters = repmat(TempCenters, [SubAmount(1) 1 1]);

            for j=1:length(Center)
                TriDist = TempCenters(:,:,j)-SubTri;
                x = TriDist(:,1,:);
                y = TriDist(:,2,:);
                z = 2*dist+TempCenters(:,3,j)+SubTri(:,3,:);
                 
                rho = sqrt(x.^2+y.^2);
                phi = atan2(y,x);
                    
                if max(max(tab_z)) < max(max(max(z)))
                   hej=4;
                end
                Gzz = interp2(tab_z, tab_r, GItabzz, z, rho, 'spline'); 
                Gzr = interp2(tab_z, tab_r, GItabzr, z, rho, 'spline'); 
                Grr = interp2(tab_z, tab_r, GItabrr, z, rho, 'spline'); 
                Gpp = interp2(tab_z, tab_r, GItabpp, z, rho, 'spline');                
                Grz = -Gzr;

                GIxx(:,j,:) = ((sin(phi)).^2).*Gpp+((cos(phi)).^2).*Grr;
                GIxy(:,j,:) = sin(phi).*cos(phi).*(Grr-Gpp);
                GIxz(:,j,:) = cos(phi).*Grz;
                GIyx(:,j,:) = sin(phi).*cos(phi).*(Grr-Gpp);
                GIyy(:,j,:) = ((cos(phi)).^2).*Gpp+((sin(phi)).^2).*Grr;
                GIyz(:,j,:) = sin(phi).*Grz;
                GIzx(:,j,:) = cos(phi).*Gzr;
                GIzy(:,j,:) = sin(phi).*Gzr;
                GIzz(:,j,:) = Gzz;
            end
        end
     
        end
end

