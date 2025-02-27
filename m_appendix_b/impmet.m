function [Z]=       impmet( EdgesTotal,TrianglesTotal,...
                            EdgeLength,K,...
                            Center,Center_,...
                            TrianglePlus,TriangleMinus,...
                            RHO_P,RHO_M,...
                            RHO__Plus,RHO__Minus,...
                            FactorA,FactorFi,Integral)   
%IMPMET Standard impedance matrix (metal surface)
%
%	Returns the complex impedance matrix [EdgesTotal x EdgesTotal]
%	Uses 9 integration points for every triangle 
%   (barycentric subdivision)
%
%   The impedance matrix is calculated as a sum of the contributions
%   due to separate triangles (similar to the "face-pair" method). 
%   See Appendix B for a detailed algorithm.
% 
%   A 9-point quadrature is used for all integrals, including 
%   the self-coupling terms (Chapters 1-10). The alternative source 
%   code with the analytical approximation of the self-coupling terms 
%   is given here. The difference between two methods is not significant. 
%
%   Copyright 2002 AEMM. Revision 2002/03/26 
%   Chapter 2/Appendix B

%Memory allocation
Z   =zeros  (EdgesTotal,EdgesTotal)+1i*zeros(EdgesTotal,EdgesTotal);

%Loop over integration triangles
for p=1:TrianglesTotal
    
    Plus     =find(TrianglePlus-p==0);
    Minus    =find(TriangleMinus-p==0);
    
    D=Center_-repmat(Center(:,p),[1 9 TrianglesTotal]); %[3 9 EdgesTotal]
    R=sqrt(sum(D.*D));                               %[1 9 TrianglesTotal]
      
    %Block for self-coupling terms
    Index=1:TrianglesTotal;
    Index(p)=[];
    g(:,:,Index) = exp(-K*R(:,:,Index))./R(:,:,Index);
    g(:,:,p)     = -K+Integral(p);
           
    gP=g(:,:,TrianglePlus);                         %[1 9 EdgesTotal]
    gM=g(:,:,TriangleMinus);                        %[1 9 EdgesTotal]
        
    Fi=sum(gP)-sum(gM);                             %[1 1 EdgesTotal]
    ZF= FactorFi.*reshape(Fi,EdgesTotal,1);
        
    for k=1:length(Plus)
        n=Plus(k);
        RP=repmat(RHO__Plus(:,:,n),[1 1 EdgesTotal]);         %[3 9 EdgesTotal]
        A=sum(gP.*sum(RP.*RHO_P))+sum(gM.*sum(RP.*RHO_M));
        Z1= FactorA.*reshape(A,EdgesTotal,1);
        Z(:,n)=Z(:,n)+EdgeLength(n)*(Z1+ZF);
    end
    for k=1:length(Minus)
        n=Minus(k);
        RP=repmat(RHO__Minus(:,:,n),[1 1 EdgesTotal]);        %[3 9 EdgesTotal]
        A=sum(gP.*sum(RP.*RHO_P))+sum(gM.*sum(RP.*RHO_M));
        Z1= FactorA.*reshape(A,EdgesTotal,1);    
        Z(:,n)=Z(:,n)+EdgeLength(n)*(Z1-ZF); 
    end
end