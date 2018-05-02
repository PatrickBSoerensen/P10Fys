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
                    [aa1] = ArbitraryAntenna.SubTrianglesIttMet(C1,C5,r1, a1);
                    [aa2] = ArbitraryAntenna.SubTrianglesIttMet(C1,C2,M, a2);           
                    [aa3] = ArbitraryAntenna.SubTrianglesIttMet(C2,C3,r2, a3); 
                    [aa4] = ArbitraryAntenna.SubTrianglesIttMet(C2,C3,M, a4);                    
                    [aa5] = ArbitraryAntenna.SubTrianglesIttMet(C3,C4,M, a5); 
                    [aa6] = ArbitraryAntenna.SubTrianglesIttMet(C1,C5,M, a6);
                    [aa7] = ArbitraryAntenna.SubTrianglesIttMet(C5,C6,M, a7); 
                    [aa8] = ArbitraryAntenna.SubTrianglesIttMet(C4,C6,M, a8);
                    [aa9] = ArbitraryAntenna.SubTrianglesIttMet(C4,C6,r3, a9); 
                    
                    SubTri(:,:,i)=...
                    [aa1 aa2 aa3 aa4 aa5 aa6 aa7 aa8 aa9];
                else
                    SubTri(:,:,i)=...
                    [a1 a2 a3 a4 a5 a6 a7 a8 a9];
                end
                
            end
        end
        
        function [SubTri, SubTriArea] = SubTrianglesIttMet(r1,r2,r3, Center)       
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
        
        function I2 = SelfTerm(p, t)
            % Method for calculating self coupling terms
            TotTri = length(t);
            I2 = 1:TotTri;
            % Looping through triangles
            for i=1:TotTri
                % vertices coordinates
                v1 = p(t(i,1),:);
                v2 = p(t(i,2),:);
                v3 = p(t(i,3),:);
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
            end
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
                    SubP = reshape(SubP, 3, []).';
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
                    SubM = reshape(SubM, 3, []).';
                    RhoM_(:,:,TM) = BasisM(SubM);
                end        
            end
        end
             
        function [PlusTri, MinusTri] = PMTri(t, EdgeList)
            for y=1:length(EdgeList)
                %Finding the points of the plus and minus triangle
                FaceEdge(1:2) = EdgeList(y,1:2);
                FaceEdgeP = [FaceEdge EdgeList(y,3)];
                FaceEdgeP = sort(FaceEdgeP,2);
                FaceEdgeM = [FaceEdge EdgeList(y,4)];
                FaceEdgeM = sort(FaceEdgeM,2);
                % Finding the index of plus and minus triangles for
                % the current basis for the outer
                PlusTri(y) = find(sum(abs(t-FaceEdgeP),2)==0);
                MinusTri(y) = find(sum(abs(t-FaceEdgeM),2)==0); 
            end
        end
        
        function [Z, b, a] = MoM(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, x, y, z)
            % alocating space
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
            Ei(:,1) = x.*exp(-1i*k.*(Center(:,1)));
            Ei(:,2) = y.*exp(-1i*k.*(Center(:,1)));
            Ei(:,3) = z.*exp(-1i*k.*(Center(:,1)));
            b1 = 1:length(EdgeList);
            b2 = 1:length(EdgeList);
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            
            SubAmount = size(SubTri);
            Quad = SubAmount(2)/3;
            % Outer loop over triangles
            for y=1:length(t)
                PO = find(PlusTri - y ==0);
                MO = find(MinusTri - y ==0);
                % Inner triangle loop
                SPO = reshape(SubTri(:,:,y),3,[]).';
                SMO = reshape(SubTri(:,:,y),3,[]).';
                for h=1:length(t)
                    PI = find(PlusTri - h ==0);
                    MI = find(MinusTri - h ==0);
                    SPI = reshape(SubTri(:,:,h),3,[]).';
                    SMI = reshape(SubTri(:,:,h),3,[]).';
                               
                    ppo = sqrt(sum((Center(y,:)-SPI).^2,2));
                    mpo = sqrt(sum((Center(y,:)-SPI).^2,2));
                    pmo = sqrt(sum((Center(y,:)-SMI).^2,2));
                    mmo = sqrt(sum((Center(y,:)-SMI).^2,2));
                           
                    ppi = sqrt(sum((Center(h,:)-SPO).^2,2));
                    mpi = sqrt(sum((Center(h,:)-SPO).^2,2));
                    pmi = sqrt(sum((Center(h,:)-SMO).^2,2));
                    mmi = sqrt(sum((Center(h,:)-SMO).^2,2));
                    %Loop over outer triangles basis functions
                    for i = 1:length(PO)
                        % Intermediate loading of plus and minus functions
                        % evaluated in center points, _ denotes sub
                        % triangle
                        zMP = (RhoP(PO(i),:));
                        zMP_ = (RhoP_(:,:,PO(i)));
                        % Subtriangles for the inner triangles basis
                        % functions
                        for j = 1:length(PI)
                            % Intermediate loading of plus and minus functions
                            % evaluated in center points, _ denotes sub
                            % triangle
                            zNP = (RhoP(PI(j),:));
                            zNP_ = (RhoP_(:,:,PI(j)));                   
                            if PO(i)==PI(j)
                                g = I2(y);
                                Z(PO(i), PI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *((dot(zMP, zNP)/4-1/k^2) * g)...
                                + Z(PO(i), PI(j));
                            else                                
                                gPPo = exp(1i.*k.*ppo)./ppo;
                                gPPi = exp(1i.*k.*ppi)./ppi;
                                
                                Z(PO(i), PI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *sum((dot(repmat(zMP,[Quad,1]), zNP_,2)/4-1/k^2) .* gPPo/Quad)...
                                + Z(PO(i), PI(j));
                            
                                Z(PO(i), PI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *sum((dot(zMP_, repmat(zNP,[Quad,1]),2)/4-1/k^2) .* gPPi/Quad)...
                                + Z(PO(i), PI(j));
                            end
                        end
                        for j=1:length(MI)
                            zNM = (RhoM(MI(j),:));
                            zNM_ = (RhoM_(:,:,MI(j)));                             
                            if PO(i)==MI(j)
                                g = I2(y);
                                Z(PO(i), MI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *((dot(zMP, zNM)/4-1/k^2) * g)...
                                + Z(PO(i), MI(j));
                            else   
                                gPMo = exp(1i.*k.*pmo)./pmo;
                                gMPi = exp(1i.*k.*mpi)./mpi;
                                
                                Z(PO(i), MI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *sum((dot(repmat(zMP,[Quad,1]), zNM_ ,2)/4+1/k^2) .* gPMo/Quad)...
                                + Z(PO(i), MI(j));
                                
                                Z(PO(i), MI(j)) = ...
                                (BasisLA(PO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *sum((dot(zMP_, repmat(zNM,[Quad,1]),2)/4+1/k^2) .* gMPi/Quad)...
                                + Z(PO(i), MI(j));
                            end    
                        end  
                        b1(PO(i),:) = sum(sum(Ei(y,:).*RhoP_(:,:,PO(i)).*BasisLA(PO(i),2)/2,2)/Quad);
                    end
                    for i=1:length(MO)
                        zMM = (RhoM(MO(i),:));
                        zMM_ = (RhoM_(:,:,MO(i)));
                        for j=1:length(PI)
                            if MO(i)==PI(j)
                                g = I2(y); 
                                Z(MO(i), PI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *((dot(zMM, zNP)/4-1/k^2) * g)...
                                + Z(MO(i), PI(j));   
                            else
                                gPMi = exp(1i.*k.*pmi)./pmi;
                                gMPo = exp(1i.*k.*mpo)./mpo;
                                   
                                Z(MO(i), PI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *sum((dot(repmat(zMM,[Quad,1]), zNP_,2)/4+1/k^2) .* gMPo/Quad)...
                                + Z(MO(i), PI(j));
                            
                                Z(MO(i), PI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(PI(j),2))/(4*pi)...
                                *sum((dot(zMM_, repmat(zNP,[Quad,1]),2)/4+1/k^2) .* gPMi/Quad)...
                                + Z(MO(i), PI(j));  
                            end 
                        end
                        for j=1:length(MI)
                            if MO(i)==MI(j)
                                g = I2(y);
                                Z(MO(i), MI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *((dot(zMM, zNM)/4-1/k^2) * g)...
                                + Z(MO(i), MI(j));        
                            else
                                gMMo = exp(1i.*k.*mmo)./mmo;
                                gMMi = exp(1i.*k.*mmi)./mmi;
                                    
                                Z(MO(i), MI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *sum((dot(repmat(zMM,[Quad,1]), zNM_,2)/4-1/k^2) .* gMMo/Quad)...
                                + Z(MO(i), MI(j));
                            
                                Z(MO(i), MI(j)) = ...
                                (BasisLA(MO(i),2)*BasisLA(MI(j),2))/(4*pi)...
                                *sum((dot(zMM_, repmat(zNM,[Quad,1]),2)/4-1/k^2) .* gMMi/Quad)...
                                + Z(MO(i), MI(j));
                            end 
                        end
                        b2(MO(i),:) = sum(sum(Ei(y,:).*RhoM_(:,:,MO(i)).*BasisLA(MO(i),2)/2,2)/Quad);
                    end      
                end       
            end
            b = b1+b2;
            % Z\b is a newer faster version of inv(Z)*b
            a = Z\b;
        end
      
        function [Z, a, b] = MoMVectorized(t, EdgeList, BasisLA, RhoP, RhoM, RhoP_, RhoM_, I2, Center, k, SubTri, x, y, z)
            % alocating space
            Z = zeros(length(EdgeList),length(EdgeList))+1i*zeros(length(EdgeList),length(EdgeList));
            TotTri = length(t);
            SubTri = permute(reshape(SubTri, 3, [], TotTri),[2 1 3]);
            SubAmount = size(SubTri);
            Ei = zeros(size(t));
            Ei(:,1) = x.*exp(-1i*k.*(Center(:,1)));
            Ei(:,2) = y.*exp(-1i*k.*(Center(:,1)));
            Ei(:,3) = z.*exp(-1i*k.*(Center(:,1)));
            EdgesTotal = length(EdgeList);
                            
            for m=1:EdgesTotal
                RhoPRep(:,:,m)=repmat(RhoP(m,:),[SubAmount(1) 1]);   %[3 9 EdgesTotal]
                RhoMRep(:,:,m)=repmat(RhoM(m,:),[SubAmount(1) 1]);  %[3 9 EdgesTotal]
            end
           
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);

            for y=1:TotTri
                Plus     =find(PlusTri-y==0);
                Minus    =find(MinusTri-y==0);
                
                D=SubTri-permute(reshape(repmat(Center(y,:),[1 SubAmount(1) TotTri]), 3, [], TotTri),[2 1 3]); %[9 3 TrianglesTotal]
                R=sqrt(sum(D.*D,2));                               %[9 1 TrianglesTotal]
      
                D1=SubTri(:,:,y)-reshape(Center,[1, 3, TotTri]); %[9 3 TrianglesTotal]
                R1=sqrt(sum(D1.*D1,2));                               %[9 1 TrianglesTotal]
                
                %Block for self-coupling terms
                Index=1:TotTri;
                Index(y)=[];
                g(:,:,Index) = exp(1i*k*R(:,:,Index))./R(:,:,Index);
                g(:,:,y)     = I2(y);
                
                g1(:,:,Index) = exp(1i*k*R1(:,:,Index))./R1(:,:,Index);
                g1(:,:,y)     = I2(y);
           
                gP=g(:,:,PlusTri);           %[9 1 EdgesTotal]
                gM=g(:,:,MinusTri);          %[9 1 EdgesTotal]
                
                gP1=g1(:,:,PlusTri);         %[9 1 EdgesTotal]
                gM1=g1(:,:,MinusTri);        %[9 1 EdgesTotal]
        
                for i=1:length(Plus)
                    n=Plus(i);
                    L = BasisLA(n,2).*BasisLA(:,2)/(4*pi);
                    pp = sum(sum(RhoPRep(:,:,n).*RhoP_/4-1/(k^2).*gP/SubAmount(1)));
                    pm = sum(sum(RhoPRep(:,:,n).*RhoM_/4+1/(k^2).*gM/SubAmount(1)));
                    Z(:,n)=Z(:,n)+L.*reshape(pp+pm,EdgesTotal,1);
                    
                    pp = sum(sum((RhoPRep.*RhoP_(:,:,n)/4-1/(k^2)).*gP1/SubAmount(1)));
                    pm = sum(sum((RhoMRep.*RhoP_(:,:,n)/4+1/(k^2)).*gM1/SubAmount(1)));
                    Z(:,n)=Z(:,n)+(L.*reshape(pp+pm,EdgesTotal,1));
                end
                for i=1:length(Minus)
                    n=Minus(i);
                    L = BasisLA(n,2).*BasisLA(:,2)/(4*pi);
                    mp = sum(sum((RhoMRep(:,:,n).*RhoP_/4+1/(k^2)).*gP/SubAmount(1)));
                    mm = sum(sum((RhoMRep(:,:,n).*RhoM_/4-1/(k^2)).*gM/SubAmount(1)));
                    Z(:,n)=Z(:,n)+L.*reshape(mp+mm,EdgesTotal,1);
                    
                    mp = sum(sum((RhoPRep.*RhoM_(:,:,n)/4+1/(k^2)).*gP1/SubAmount(1)));
                    mm = sum(sum((RhoMRep.*RhoM_(:,:,n)/4-1/(k^2)).*gM1/SubAmount(1)));
                    Z(:,n)=Z(:,n)+(L.*reshape(mp+mm,EdgesTotal,1));
                end 
            end

            for m=1:EdgesTotal
                b1 =sum(sum(Ei(PlusTri(m),:).*RhoP_(:,:,m),2)/SubAmount(1));
                b2 =sum(sum(Ei(MinusTri(m),:).*RhoM_(:,:,m),2)/SubAmount(1));
                b(m)=BasisLA(m,2)*(b1+b2)/2;
            end
            b = b.';
            %System solution
            a=Z\b;
          end
          
        function [Jface] = CurrentCalc(t, EdgeList, w, mu0, a, BasisLA, RhoP, RhoM, RhoP_, RhoM_, sub)
            Jface = zeros(size(t));
            
            [PlusTri, MinusTri] = ArbitraryAntenna.PMTri(t, EdgeList);
            for n=1:length(t)
                Plus     =find(PlusTri-n==0);
                Minus    =find(MinusTri-n==0);
                Quad = size(RhoP_);
                Quad = Quad(1);
                for y=1:length(Plus)
                    if sub
                        Jface(n,:) = sum(1i.*w.*mu0*a(Plus(y))*BasisLA(Plus(y),2)*RhoP_(:,:,Plus(y))/(2*Quad)) + Jface(n,:);
                    else 
                        Jface(n,:) = 1i.*w.*mu0*a(Plus(y))*BasisLA(Plus(y),1)*RhoP(Plus(y),:)/2 + Jface(n,:);    
                    end
                end
                for y=1:length(Minus)
                    if sub
                        Jface(n,:) = sum(1i.*w.*mu0*a(Minus(y))*BasisLA(Minus(y),2)*RhoM_(:,:,Minus(y))/(2*Quad)) + Jface(n,:);
                    else
                        Jface(n,:) = 1i.*w.*mu0*a(Minus(y))*BasisLA(Minus(y),3)*RhoM(Minus(y),:)/2 + Jface(n,:);    
                    end
                end    
            end
        end
        
        function [Exy, Exz, Ezy, xrange, yrange, zrange, Exyx, Exzx, Eyzx, Exyy, Exzy, Eyzy, Exyz, Exzz, Eyzz] = EField(Center, k, J,...
                xmin, xmax, zmin, zmax, ymin, ymax, steps, Area)
       
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
                                                        
                            g = exp(1i.*k.*r)./(4*pi*r);
                            G1 = (1+1i./(k*r)-1./(k*r).^2);
                            G2 = (1+3i./(k*r)-3./(k*r).^2);
                            
                            RR = Rx.*Rx;
                            Gxx = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rx;
                            Gxy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rx;
                            Gxz = (-(RR./r.^2).*G2).*g;
                            Exyx = Exyx + (Gxx .* J(i,1) + Gxy .* J(i,2) + Gxz .* J(i,3)).*Area(i);
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Exyy = Exyy + (Gyx .* J(i,1) + Gyy .* J(i,2) + Gyz .* J(i,3)).*Area(i);
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Exyz = Exyz + (Gzx .* J(i,1) + Gzy .* J(i,2) + Gzz .* J(i,3)).*Area(i);
                        elseif j==2
                            %xz
                            Rx = repmat(rx(i,:)',1,steps);
                            Ry = repmat(ry(i,:),steps,steps);
                            Rz = repmat(rz(i,:),steps,1);
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
                            Exzx = Exzx + (Gxx .* J(i,1) + Gxy .* J(i,2) + Gxz .* J(i,3)).*Area(i);
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Exzy = Exzy + (Gyx .* J(i,1) + Gyy .* J(i,2) + Gyz .* J(i,3)).*Area(i);
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Exzz = Exzz + (Gzx .* J(i,1) + Gzy .* J(i,2) + Gzz .* J(i,3)).*Area(i);
                        else
                            %yz
                            Rx = repmat(rx(i,:),steps,steps);
                            Ry = repmat(ry(i,:),steps,1);
                            Rz = repmat(rz(i,:)',1,steps);
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
                            Eyzx = Eyzx + (Gxx .* J(i,1) + Gxy .* J(i,2) + Gxz .* J(i,3)).*Area(i);
                            
                            RR = Ry.*Rx;
                            Gyx = (-(RR./r.^2).*G2).*g;
                            RR = Ry.*Ry;
                            Gyy = (G1-(RR./r.^2).*G2).*g;
                            RR = Ry.*Rz;
                            Gyz = (-(RR./r.^2).*G2).*g;
                            Eyzy = Eyzy + (Gyx .* J(i,1) + Gyy .* J(i,2) + Gyz .* J(i,3)).*Area(i);
                            
                            RR = Rz.*Rx;
                            Gzx = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Ry;
                            Gzy = (-(RR./r.^2).*G2).*g;
                            RR = Rz.*Rz;
                            Gzz = (G1-(RR./r.^2).*G2).*g;
                            Eyzz = Eyzz + (Gzx .* J(i,1) + Gzy .* J(i,2) + Gzz .* J(i,3)).*Area(i);
                        end
                end 
            end
            Exy = sqrt(Exyx.^2+Exyy.^2+Exyz.^2);
            Exz = sqrt(Exzx.^2+Exzy.^2+Exzz.^2);
            Ezy = sqrt(Eyzx.^2+Eyzy.^2+Eyzz.^2);
        end
        
        function [] = PolarPlot(E, xmin, xmax,  ymin, ymax, steps, dist)
            xrange = linspace(xmin, xmax, steps);
            yrange = linspace(ymin, ymax, steps);
            Rx = repmat(xrange',1,steps);
            Ry = repmat(yrange,steps,1);
            [theta,rho] = cart2pol(Rx,Ry);
            angles = linspace(0, 2*pi, steps);
            Look = angles*dist;
            PlotScheme = find(round(sqrt(Rx.^2+Ry.^2)-Look)==0);
            EPlot = E(PlotScheme);
            angles = linspace(0, 2*pi, length(EPlot));
            figure(2)
            polarplot(angles,abs(EPlot));
        end
    end    
end

