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
         
        function [EdgeList, Basis, BasisLA] = BasisFunc(p, t, ConnectCell)
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
%                     SubP = reshape(SubP, 3, []).';
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
%                     SubM = reshape(SubM, 3, []).';
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
            end
        end
            
        function [Ei, EdgeV] = VoltageFeed(t, p, Center, FeedPos, v, EdgeList, BasisLA) 
            
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            Ei = zeros(size(Center));
            EdgeV = 1:length(EdgeList);
            EdgeV(:) = 0;
            
            %% Based on EdgeCenter
            PlusCoords = Center(PlusTri, :);
            MinusCoords = Center(MinusTri, :);
            
%             EdgeCenter = (PlusCoords+MinusCoords)/2;
%             FeedEdges = abs(EdgeCenter(:,2))<=1e-18+FeedPos(2);
            
            %% straight edge selection mb?
            FeedEdges = [];
            counter=1;
            PointsOnFeed = find(abs(p(:,2))<=1e-18+FeedPos(2));
            for i = 1:length(PointsOnFeed)
                    StartPoint = EdgeList(:,1)==PointsOnFeed(i);
                for j=i+1:length(PointsOnFeed)
                    EndPoint = EdgeList(:,2)==PointsOnFeed(j);
                    if sum(StartPoint.*EndPoint)
                        FeedEdges(counter) = find(StartPoint.*EndPoint);
                        counter = counter+1;
                    end
                end
            end
%             v=v*Diam;
            EdgeV(FeedEdges)=BasisLA(FeedEdges,2)*v;
            if ~isempty(FeedEdges)
                Ei(PlusTri(FeedEdges),2) = v/length(FeedEdges);
                Ei(MinusTri(FeedEdges),2) = v/length(FeedEdges);
            end
        end
        
        function [Z, b, a] = MoM(w, mu, t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Center, k, SubTri, x, y, z, Point, Ei,...
                distx, Reflector, InEps, strip_length, strip_width, dx, Nz, lambda, n, eps1)
            % alocating space
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
            if ~Point
                Ei(:,1) = x.*exp(-1i*k.*(Center(:,1)));
                Ei(:,2) = y.*exp(-1i*k.*(Center(:,1)));
                Ei(:,3) = z.*exp(-1i*k.*(Center(:,1)));
            end
            if Reflector
                RefCoef = (1-n)/(1+n);
                Ei(:,1) = Ei(:,1)+x.*exp(1i*k.*(Center(:,1))).*RefCoef;
                Ei(:,2) = Ei(:,2)+y.*exp(1i*k.*(Center(:,1))).*RefCoef;
                Ei(:,3) = Ei(:,3)+z.*exp(1i*k.*(Center(:,1))).*RefCoef;
            end
            b = 1:length(EdgeList); b(:) =0; b=b';
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            
            if Reflector
            [GIx, GIy, GIz] = ArbitraryAntenna.IDGreens(k, distx, strip_length, strip_width, dx, Nz, lambda, n, InEps, eps1, Center, SubTri);
            end
                        
            SubAmount = size(SubTri);
            Quad = SubAmount(1);
            % Outer loop over triangles
            for y=1:length(t)
                PO = find(PlusTri - y ==0);
                MO = find(MinusTri - y ==0);
                % Inner triangle loop
                SO = SubTri(:,:,y);
                for h=1:length(t)
                    PI = find(PlusTri - h ==0);
                    MI = find(MinusTri - h ==0);
                    SI = SubTri(:,:,h);
                   
                    IoO = sqrt(sum((Center(y,:)-SI).^2,2));
                    OoI = sqrt(sum((Center(h,:)-SO).^2,2));
 
                     if Reflector
                        GIppo = GIx(:,:,y)+GIy(:,:,y)+GIz(:,:,y);
                        GImpo = GIx(:,:,y)+GIy(:,:,y)+GIz(:,:,y);
                        GIpmo = GIx(:,:,y)+GIy(:,:,y)+GIz(:,:,y);
                        GImmo = GIx(:,:,y)+GIy(:,:,y)+GIz(:,:,y);
                    
                        GIppi = GIx(:,:,h)+GIy(:,:,h)+GIz(:,:,h);
                        GImpi = GIx(:,:,h)+GIy(:,:,h)+GIz(:,:,h);
                        GIpmi = GIx(:,:,h)+GIy(:,:,h)+GIz(:,:,h);
                        GImmi = GIx(:,:,h)+GIy(:,:,h)+GIz(:,:,h);
                    end
                    %Loop over outer triangles basis functions
                    for i = 1:length(PO)
                        % Intermediate loading of plus and minus functions
                        % evaluated in center points, _ denotes sub
                        % triangle
                        zMP = (RhoP(PO(i),:));
                        zMP_ = (RhoP_(:,:,PO(i)));
                        zMPR = repmat(zMP ,Quad,1);
                        
                        % Subtriangles for the inner triangles basis
                        % functions
                        for j = 1:length(PI)
                            % Intermediate loading of plus and minus functions
                            % evaluated in center points, _ denotes sub
                            % triangle
                            zNP = (RhoP(PI(j),:));
                            zNP_ = (RhoP_(:,:,PI(j)));           
                            zNPR = repmat(zNP, Quad,1);
                        
                            gPPo = exp(1i.*k.*IoO)./IoO;
                            gPPi = exp(1i.*k.*OoI)./OoI;
                                                  
                            Z(PO(i), PI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *sum((dot(zMPR, zNP_,2)/4-1/k^2) .* gPPo ...
                                +(dot(zMP_, zNPR,2)/4-1/k^2) .* gPPi)/Quad...
                                + Z(PO(i), PI(j));

                            if Reflector
                                Z(PO(i), PI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(PI(j),2))...
                                *sum((dot(repmat(zMP,[Quad,1]).* GIppo, zNP_.* GIppo,2)) /Quad)...
                                + Z(PO(i), PI(j));
                            
                                Z(PO(i), PI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(PI(j),2))...
                                *sum((dot(zMP_.*GIppi, repmat(zNP,[Quad,1]).*GIppi,2)) /Quad)...
                                + Z(PO(i), PI(j));
                            end
                        end
                        for j=1:length(MI)
                            zNM = (RhoM(MI(j),:));
                            zNM_ = (RhoM_(:,:,MI(j)));        
                            zNMR = repmat(zNM ,Quad,1);
                        
                            gPMo = exp(1i.*k.*IoO)./IoO;
                            gMPi = exp(1i.*k.*OoI)./OoI;
                          
                            Z(PO(i), MI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *sum((dot(zMPR, zNM_ ,2)/4+1/k^2) .* gPMo...
                                + (dot(zMP_, zNMR,2)/4+1/k^2) .* gMPi)/Quad...
                                + Z(PO(i), MI(j));
                        
                          if Reflector
                            Z(PO(i), MI(j)) = ...
                            (BasisLA(PO(i),2)*BasisLA(MI(j),2))...
                            *sum((dot(repmat(zMP,[Quad,1]).* GIpmo, zNM_.* GIpmo,2)) /Quad)...
                            + Z(PO(i), MI(j));
                            
                            Z(PO(i), MI(j)) = ...
                            (BasisLA(PO(i),2)*BasisLA(MI(j),2))...
                            *sum((dot(zMP_.*GImpi, repmat(zNM,[Quad,1]).*GImpi,2)) /Quad)...
                            + Z(PO(i), MI(j));
                        
                          end
                        end
                    end
                    for i=1:length(MO)
                        zMM = (RhoM(MO(i),:));
                        zMM_ = (RhoM_(:,:,MO(i)));
                        zMMR = repmat(zMM ,Quad,1);
                          
                        for j=1:length(PI)
                            zNP = (RhoP(PI(j),:));
                            zNP_ = (RhoP_(:,:,PI(j)));   
                            zNPR = repmat(zNP ,Quad,1);
                        
                            gMPo = exp(1i.*k.*IoO)./IoO;
                            gPMi = exp(1i.*k.*OoI)./OoI;
                                
                            Z(MO(i), PI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *sum((dot(zMMR, zNP_,2)/4+1/k^2) .* gMPo ...
                                + (dot(zMM_, zNPR,2)/4+1/k^2) .* gPMi)/Quad...
                                + Z(MO(i), PI(j));
                            
                        if Reflector
                            Z(MO(i), PI(j)) = ...
                            (BasisLA(MO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                            *sum((dot(repmat(zMM,[Quad,1]).*GImpo, zNP_.*GImpo,2))  /Quad)...
                            + Z(MO(i), PI(j));
                            
                            Z(MO(i), PI(j)) = ...
                            (BasisLA(MO(i),2)*BasisLA(PI(j),2))...
                            *sum((dot(zMM_.*GIpmi, repmat(zNP,[Quad,1]).*GIpmi,2)) /Quad)...
                            + Z(MO(i), PI(j));   
                        end
                        end
                        for j=1:length(MI)
                            zNM = (RhoM(MI(j),:));
                            zNM_ = (RhoM_(:,:,MI(j)));   
                            zNMR = repmat(zNM ,Quad,1);
                            
                            gMMo = exp(1i.*k.*IoO)./IoO;
                            gMMi = exp(1i.*k.*OoI)./OoI;
                                             
                                Z(MO(i), MI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *sum((dot(zMMR, zNM_,2)/4-1/k^2) .* gMMo...
                                +(dot(zMM_, zNMR,2)/4-1/k^2) .* gMMi)/Quad...
                                + Z(MO(i), MI(j));
                            
                        if Reflector
                            Z(MO(i), MI(j)) = ...
                            (BasisLA(MO(i),2)*BasisLA(MI(j),2))...
                            *sum((dot(repmat(zMM,[Quad,1]).* GImmo, zNM_.* GImmo,2)) /Quad)...
                            + Z(MO(i), MI(j));
                            
                            Z(MO(i), MI(j)) = ...
                            (BasisLA(MO(i),2)*BasisLA(MI(j),2))...
                            *sum((dot(zMM_.* GImmi, repmat(zNM,[Quad,1]).* GImmi,2)) /Quad)...
                            + Z(MO(i), MI(j));    
                        end
           
                        end
                     end      
                end       
            end
            for y=1:length(t)
                PO = find(PlusTri - y ==0);
                MO = find(MinusTri - y ==0);
                for i=1:length(PO)
                    b(PO(i)) = -1i/(w*mu)*sum(dot(repmat(Ei(y,:),Quad,1),RhoP_(:,:,PO(i)),2).*BasisLA(PO(i),2)/Quad)/2 + b(PO(i));
                end
                for i=1:length(MO)
                    b(MO(i)) = -1i/(w*mu)*sum(dot(repmat(Ei(y,:),Quad,1),RhoM_(:,:,MO(i)),2).*BasisLA(MO(i),2)/Quad)/2 + b(MO(i));
                end
            end
            % Z\b is a newer faster version of inv(Z)*b
            a = Z\b;
        end
         
        function [Z, b, a] = MoMSergey(w, mu, t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Center, k, SubTri, x, y, z, Point, Ei,...
                distx, Reflector, InEps, strip_length, strip_width, dx, Nz, lambda, n, eps0)
            
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
            if ~Point
                Ei(:,1) = x.*exp(-1i*k.*(Center(:,2)));
                Ei(:,2) = y.*exp(-1i*k.*(Center(:,1)));
                Ei(:,3) = z.*exp(-1i*k.*(Center(:,1)));
            end
            if Reflector
                RefCoef = (1-n)/(1+n);
                Ei(:,1) = Ei(:,1)+x.*exp(-1i*k.*(Center(:,2))).*RefCoef;
                Ei(:,2) = Ei(:,2)+y.*exp(-1i*k.*(Center(:,1))).*RefCoef;
                Ei(:,3) = Ei(:,3)+z.*exp(-1i*k.*(Center(:,1))).*RefCoef;
             
                [GIx, GIy, GIz] = ArbitraryAntenna.IDGreens(k, distx, strip_length, strip_width, dx, Nz, lambda, n, InEps, eps0, Center, SubTri);
                
                GI = permute(GIx+GIy+GIz,[2,1,3]);
                GI = repmat(GI,3,1);
            end
            
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
                    if Reflector
                        gmPnP = exp(-1i*k*mPdist(:,:,PlusTri(n)))./mPdist(:,:,PlusTri(n))+GI(:,:,PlusTri(n));
                        gmMnP = exp(-1i*k*mMdist(:,:,PlusTri(n)))./mMdist(:,:,PlusTri(n))+GI(:,:,PlusTri(n));
                
                        gmPnM = exp(-1i*k*mPdist(:,:,MinusTri(n)))./mPdist(:,:,MinusTri(n))+GI(:,:,MinusTri(n));
                        gmMnM = exp(-1i*k*mMdist(:,:,MinusTri(n)))./mMdist(:,:,MinusTri(n))+GI(:,:,MinusTri(n));
                    else
                        gmPnP = exp(-1i*k*mPdist(:,:,PlusTri(n)))./mPdist(:,:,PlusTri(n));
                        gmMnP = exp(-1i*k*mMdist(:,:,PlusTri(n)))./mMdist(:,:,PlusTri(n));
                
                        gmPnM = exp(-1i*k*mPdist(:,:,MinusTri(n)))./mPdist(:,:,MinusTri(n));
                        gmMnM = exp(-1i*k*mMdist(:,:,MinusTri(n)))./mMdist(:,:,MinusTri(n));
                    end
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
        
        function [Z, a, b] = MoMVectorized(w,mu,t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, Center, k, SubTri, x, y, z, Point, Ei,...
                distx, Reflector, InEps, strip_length, strip_width, dx, Nz, lambda, n, eps0)
            % alocating space
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
            
            SubAmount = size(SubTri);
            Quad = SubAmount(1);
            if ~Point
                Ei(:,1) = x.*exp(1i*k.*(Center(:,2)));
                Ei(:,2) = y.*exp(1i*k.*(Center(:,1)));
                Ei(:,3) = z.*exp(1i*k.*(Center(:,1)));
            end
            if Reflector
                RefCoef = (1-n)/(1+n);
                if ~Point
                    Ei(:,1) = Ei(:,1)+x.*exp(-1i*k.*(Center(:,2))).*RefCoef;
                    Ei(:,2) = Ei(:,2)+y.*exp(-1i*k.*(Center(:,1)+distx)).*RefCoef;
                    Ei(:,3) = Ei(:,3)+z.*exp(-1i*k.*(Center(:,1))).*RefCoef;
                else
                    Ei(:,1) = Ei(:,1)+Ei(:,1).*RefCoef;
                    Ei(:,2) = Ei(:,2)+Ei(:,2).*RefCoef;
                    Ei(:,3) = Ei(:,3)+Ei(:,3).*RefCoef;
                end
                [GIx, GIy, GIz] = ArbitraryAntenna.IDGreens(k, distx, strip_length, strip_width, dx, Nz, lambda, n, InEps, eps0, Center, SubTri);
  
                GI = GIx+GIy+GIz;

            end
            
            EdgesTotal = length(EdgeList);
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);

            for m=1:EdgesTotal
                mPdist = sqrt(sum((Center(PlusTri(m),:)-SubTri(:,:,:)).^2,2));
                mMdist = sqrt(sum((Center(MinusTri(m),:)-SubTri(:,:,:)).^2,2));
                rhomP = repmat(RhoP(m,:),length(EdgeList),1);
                rhomM = repmat(RhoM(m,:),length(EdgeList),1);
                    
                rhonP_ = RhoP_(:,:,:);
                rhonM_ = RhoM_(:,:,:);
                if Reflector 
                    gmPnPR = exp(-1i*k*(mPdist(:,:,PlusTri)+distx))./(mPdist(:,:,PlusTri)+distx);
                    gmMnPR = exp(-1i*k*(mMdist(:,:,PlusTri)+distx))./(mMdist(:,:,PlusTri)+distx);
                
                    gmPnMR = exp(-1i*k*(mPdist(:,:,MinusTri)+distx))./(mPdist(:,:,MinusTri)+distx);
                    gmMnMR = exp(-1i*k*(mMdist(:,:,MinusTri)+distx))./(mMdist(:,:,MinusTri)+distx);
                    
                    gmPnP = exp(-1i*k*mPdist(:,:,PlusTri))./mPdist(:,:,PlusTri);%+gmPnPR;%+GI(:,:,PlusTri)%
                    gmMnP = exp(-1i*k*mMdist(:,:,PlusTri))./mMdist(:,:,PlusTri);%+gmMnPR;%+GI(:,:,PlusTri);%
                
                    gmPnM = exp(-1i*k*mPdist(:,:,MinusTri))./mPdist(:,:,MinusTri);%+gmPnMR;%GI(:,:,MinusTri);
                    gmMnM = exp(-1i*k*mMdist(:,:,MinusTri))./mMdist(:,:,MinusTri);%+gmMnMR;%+GI(:,:,MinusTri)
                else
                    gmPnP = exp(-1i*k*mPdist(:,:,PlusTri))./mPdist(:,:,PlusTri);
                    gmMnP = exp(-1i*k*mMdist(:,:,PlusTri))./mMdist(:,:,PlusTri);
                
                    gmPnM = exp(-1i*k*mPdist(:,:,MinusTri))./mPdist(:,:,MinusTri);
                    gmMnM = exp(-1i*k*mMdist(:,:,MinusTri))./mMdist(:,:,MinusTri);
                end    
                
                AmnP = mu/(4*pi)*(BasisLA(:,2).*permute(sum(rhonP_.*gmPnP/(2*Quad)),[3 2 1])+BasisLA(:,2).*permute(sum(rhonM_.*gmPnM/(2*Quad)),[3 2 1]));
                AmnM = mu/(4*pi)*(BasisLA(:,2).*permute(sum(rhonP_.*gmMnP/(2*Quad)),[3 2 1])+BasisLA(:,2).*permute(sum(rhonM_.*gmMnM/(2*Quad)),[3 2 1]));
            
                PhiP = -1/(4*pi*1i*w*eps0)*(BasisLA(:,2).*permute(sum(gmPnP),[3 2 1])/Quad-BasisLA(:,2).*permute(sum(gmPnM),[3 2 1])/Quad);
                PhiM = -1/(4*pi*1i*w*eps0)*(BasisLA(:,2).*permute(sum(gmMnP),[3 2 1])/Quad-BasisLA(:,2).*permute(sum(gmMnM),[3 2 1])/Quad);
                
                if Reflector
                    AmnP = mu/(4*pi)*(BasisLA(:,2).*permute(sum(rhonP_.*GI(:,:,PlusTri)/(2*Quad)),[3 2 1])+BasisLA(:,2).*permute(sum(rhonM_.*GI(:,:,PlusTri)/(2*Quad)),[3 2 1]))...
                        +AmnP;
                    AmnM = mu/(4*pi)*(BasisLA(:,2).*permute(sum(rhonP_.*GI(:,:,MinusTri)/(2*Quad)),[3 2 1])+BasisLA(:,2).*permute(sum(rhonM_.*GI(:,:,MinusTri)/(2*Quad)),[3 2 1]))...
                        +AmnM;
            
                    PhiP = -1/(8*pi*1i*w*eps0)*(BasisLA(:,2).*permute(sum(dot(rhonP_,GI(:,:,PlusTri),2)),[3 2 1])/Quad-BasisLA(:,2).*permute(sum(dot(rhonM_,GI(:,:,PlusTri),2)),[3 2 1])/Quad)...
                        +PhiP;
                    PhiM = -1/(8*pi*1i*w*eps0)*(BasisLA(:,2).*permute(sum(dot(rhonP_,GI(:,:,MinusTri),2)),[3 2 1])/Quad-BasisLA(:,2).*permute(sum(dot(rhonM_,GI(:,:,MinusTri),2)),[3 2 1])/Quad)...
                        +PhiM;   
                end
            
                Z(m,:) = BasisLA(m,2).*(1i*w*(dot(AmnP,rhomP,2)/2+dot(AmnM,rhomM,2)/2)+PhiM-PhiP);  
            end
            
            b = BasisLA(:,2).*(dot(Ei(PlusTri,:),RhoP,2)/2+dot(Ei(MinusTri,:),RhoM,2)/2);
            
            %System solution
            a=Z\b;
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
   
        function [Esc, EscPhi, EscTheta] = AngularFarField(w, mu, k, r, Center, J, steps)
            phi = linspace(-pi, 2*pi, steps)';
            theta = linspace(-pi, 2*pi, steps)';
            xH = [1, 0, 0];
            yH = [0, 1, 0];
            zH = [0, 0, 1];
            phiH = -xH.*sin(phi)+yH.*cos(phi);
            thetaH = (xH.*cos(phi)+yH.*sin(phi)).*cos(theta)-zH.*sin(theta);
            rHat = (xH.*cos(phi)+yH.*sin(phi)).*sin(theta)+zH.*cos(theta);
            EscTheta = zeros(steps);
            EscPhi = zeros(steps);
           for i=1:length(Center)
                IntegralTerm = 1i*w*mu*(exp(1i*k*r)/(4*pi*r))...
                *exp(-1i*k*sum(rHat.*Center(i,:),2)).*J(i,:); 
                           
                EscTheta = dot(thetaH,IntegralTerm,2) + EscTheta;
            
                EscPhi = dot(phiH,IntegralTerm,2) + EscPhi;
            end
            
            Esc = EscPhi.^2+EscTheta.^2;
            
            figure(4)
            plot(phi, abs(Esc).^2*r^2)
            xlabel('Phi')
            ylabel('E')
            title('EscPhi^2+EscTheta^2')
        end
        
        function [Esc, EscPhi, EscTheta] = AngularFarFieldSurf(w, mu, k, r, Center, J, steps, xdist)
            phi = linspace(-pi, 2*pi, steps)';
            theta = linspace(-pi, 2*pi, steps)';
            
            SurfAngle = acos(xdist);
            under = phi<=SurfAngle;
            over = phi>=2*pi-SurfAngle;
            UseTrans = logical(under+over);
            
            xH = [1, 0, 0];
            yH = [0, 1, 0];
            zH = [0, 0, 1];
            phiH = -xH.*sin(phi)+yH.*cos(phi);
            thetaH = (xH.*cos(phi)+yH.*sin(phi)).*cos(theta)-zH.*sin(theta);
            rHat = (xH.*cos(phi)+yH.*sin(phi)).*sin(theta)+zH.*cos(theta);
            EscTheta = zeros(steps);
            EscPhi = zeros(steps);
           for i=1:length(Center)
                NormalGreens = 1i*w*mu*(exp(1i*k*r)/(4*pi*r))...
                *exp(-1i*k*sum(rHat.*Center(i,:),2)).*J(i,:); 
                
                IndirectGreens = exp(1i*k*r)/(4*pi*r)*exp(-1i*krho*rhohat*rmark)*exp(1i*kz1*zmark)*...
                    (refS(krho)*phiH*phiH-refP(krho)*thetaH*(zH*krho/k1+rhoH*kz1/k1));
                
                TransmitGreens = 0;
                
                Greens = NormalGreens+IndirectGreens;
                Greens(UseTrans,:) = TransmitGreens;
                
                EscTheta = dot(thetaH,Greens,2) + EscTheta;
            
                EscPhi = dot(phiH,Greens,2) + EscPhi;
            end
            
            Esc = EscPhi.^2+EscTheta.^2;
            
            figure(4)
            plot(phi, abs(Esc).^2*r^2)
            xlabel('Phi')
            ylabel('E')
            title('EscPhi^2+EscTheta^2')
        end
        
        function [GIx, GIy, GIz] = IDGreens(k0, distx, ant_length, ant_width, dx, Nz, lambda, n, epsL2, eps1, Center, SubTri)
                
            SubAmount = size(SubTri);
%             SubTri = permute(SubTri,[2 1 3]);
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
%             rhomin = -ant_length/2-ant_width/2;
            rhomax = ant_length/2;%sqrt((ant_length)^2+(distx-ant_width)^2);
%             rhomax = ant_length/2+ant_width/2;
%             rhov1 = [rhomin rhomin+dx rhomin+dx*sqrt(2) rhomin+dx*2 rhomin+dx*sqrt(5) rhomin+dx*sqrt(8)];
            rhov1 = [0 dx dx*sqrt(2) dx*2 dx*sqrt(5) dx*sqrt(8)];
            
            if rhomax>dx*3
                rhov2=linspace(dx*3,rhomax,ceil((rhomax-dx*3)/dx)+1);
            else
                rhov2=[];
            end
            rhotabv=[rhov1 rhov2];
            Nrho=length(rhotabv);
    
            GItabzz=zeros(Nrho,2*Nz-1); GItabzr=GItabzz;
            GItabrr=GItabzz; GItabpp=GItabzz;
            tab_z=GItabzz; tab_r=GItabzz;
            
            for jrho=1:Nrho
                rho=rhotabv(jrho);
                for jz=1:Nz*2-1
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
            AllSubDists = TempCenters+SubTri+distx;

            x = AllSubDists(:,1,:);
            y = AllSubDists(:,2,:);
            z = sum(AllSubDists(:,3,:).^2,2);
                 
            for i=1:SubAmount(1)
                rho = sqrt(x.*x+y.*y);
                phi = atan2(y,x);

                Gzz = permute(interp2(tab_z, tab_r, GItabzz, permute(z(i,:,:), [3 2 1]), permute(rho(i,:,:), [3 2 1]), 'spline'), [3 2 1]); 
                Gzr = permute(interp2(tab_z, tab_r, GItabzr, permute(z(i,:,:), [3 2 1]), permute(rho(i,:,:), [3 2 1]), 'spline'), [3 2 1]); 
                Grr = permute(interp2(tab_z, tab_r, GItabrr, permute(z(i,:,:), [3 2 1]), permute(rho(i,:,:), [3 2 1]), 'spline'), [3 2 1]); 
                Gpp = permute(interp2(tab_z, tab_r, GItabpp, permute(z(i,:,:), [3 2 1]), permute(rho(i,:,:), [3 2 1]), 'spline'), [3 2 1]);                
                Grz = -Gzr;
                
                GIxx(i,:,:) = ((sin(phi(i,:,:))).^2).*Gpp+((cos(phi(i,:,:))).^2).*Grr;
                GIxy(i,:,:) = sin(phi(i,:,:)).*cos(phi(i,:,:)).*(Grr-Gpp);
                GIxz(i,:,:) = cos(phi(i,:,:)).*Grz;
                GIyx(i,:,:) = sin(phi(i,:,:)).*cos(phi(i,:,:)).*(Grr-Gpp);
                GIyy(i,:,:) = ((cos(phi(i,:,:))).^2).*Gpp+((sin(phi(i,:,:))).^2).*Grr;
                GIyz(i,:,:) = sin(phi(i,:,:)).*Grz;
                GIzx(i,:,:) = cos(phi(i,:,:)).*Gzr;
                GIzy(i,:,:) = sin(phi(i,:,:)).*Gzr;
                GIzz(i,:,:) = Gzz;
            end
            GIxx = reshape(GIxx, SubAmount(1), [],SubAmount(3));
            GIxy = reshape(GIxy, SubAmount(1), [],SubAmount(3));
            GIxz = reshape(GIxz, SubAmount(1), [],SubAmount(3));
            GIx = [GIxx GIxy GIxz];
                    
            GIyx = reshape(GIyx, SubAmount(1),[],SubAmount(3));
            GIyy = reshape(GIyy, SubAmount(1),[],SubAmount(3));
            GIyz = reshape(GIyz, SubAmount(1),[],SubAmount(3));
            GIy = [GIyx GIyy GIyz];
            
            GIzx = reshape(GIzx, SubAmount(1),[],SubAmount(3));
            GIzy = reshape(GIzy, SubAmount(1),[],SubAmount(3));
            GIzz = reshape(GIzz, SubAmount(1),[],SubAmount(3));
            GIz = [GIzx GIzy GIzz];
        end
    
        function [] = DistJPlot(Center, J)
%             unq = unique(sqrt(Center(:,1).^2+Center(:,3).^2));
unq = unique(Center(:,1));
unq(1) = [];
            figure(96)
            hold on
            for i=1:length(unqi)
%                 index = find(sqrt(Center(:,1).^2+Center(:,3).^2)-unq(i) == 0);
index = find(Center(:,1)-unq(i)==0);     
if length(index)>10
plot(Center(index,2),abs(J(index,2)),'*')
end
            end
        end
        end
end

