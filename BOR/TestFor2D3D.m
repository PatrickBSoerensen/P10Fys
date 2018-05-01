clear t2
r=0.0031;
t1=linspace(0,r*pi/2,20).';
for i=1:20
    t2(i) = r*sin(t1(i)/r)*pi/2;
end
t2=t2.';

r_vecf =[cos(t2)/(r*sin(t1/r))*r*sin(t1/r), ...
        sin(t2/(r*sin(t1/r)))*r*sin(t1/r), ...
        r*cos(t1/r)];
   
    x = @(theta,phi) r*sin(theta)*cos(phi);
    y = @(theta,phi) r*sin(theta)*sin(phi);
    z = @(theta) r*cos(theta);
    
%     theta = linspace(0,pi/2,100);
%     phi = linspace(0,pi/2,100);
%     x = r*sin(theta).*cos(phi);
%     y = r*sin(theta).*sin(phi);
%     z = r*cos(theta);


%    rect = [3; 4; t2(1); t2(1); ...
%                 t2(end); t2(end); ...
%                 t1(1); t1(end); ...
%                 t1(end); t1(1);];
%    Circ = [1; ant.Length/2-ant.Radii; ant.Centre(2); ant.Radii];
%    Circ = [UpperCirc;zeros(length(rect) - length(UpperCirc),1)];
%             gdc = [rect, Circ];
%             ns = char('R', 'C');
%             ns = ns';
%             sf = 'R+C';
%             [dl, bt] = decsg(gdc,sf, ns);
% %             [dl, bt] = csgdel(dl,bt);
%         
%             solver = createpde;
% %             geometryFromEdges(solver,dl);
%             %0.004 er højeste for at få trekanter i enderne
%             generateMesh(solver, 'Hmin', 0.000001, 'GeometricOrder', 'quadratic');
%             figure(1)
% %             pdeplot(solver)