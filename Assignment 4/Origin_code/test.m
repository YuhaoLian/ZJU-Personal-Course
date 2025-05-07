clc;
clear;

errval = zeros(size(2:10));
for k = 2:10
    %   generate the mesh file
    rectangle3d(k);
    %   solve the equation by FEM
    [coordinates,u]=fem3dlinear;
    %   the exact solution
    uex = coordinates(:,1).^2+coordinates(:,2).^2+coordinates(:,3).^2;
    %   evaluate the error
    errval(k-1)=norm(uex - u)/norm(uex);
end
% plot the error
figure(); plot(2:10,errval,'*-')

function rectangle3d(n)
%Generate the mesh data file for a cuboid
%n=10;
xmin=0;xmax=1;
ymin=0;ymax=1;zmin=0;zmax=1;
nx=n; ny=n; nz=n;
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
z=linspace(zmin,zmax,nz);
np=0; % no. of vertices
X=zeros(2*nx*ny*nz,1);
Y=zeros(2*nx*ny*nz,1);
Z=zeros(2*nx*ny*nz,1);
for l=1:nz
    for k=1:ny
        for j=1:nx
            np = np+ 1;
            X(np) = x(j);
            Y(np) = y(k);
            Z(np) = z(l);
        end
    end
end
for l=1:nz-1
    for k=1:ny-1
        for j=1:nx-1
            np = np+ 1;
            X(np) = (x(j)+x(j+1))/2;
            Y(np) = (y(k)+y(k+1))/2;
            Z(np) = (z(l)+z(l+1))/2;
        end
    end
end
tri=delaunayTriangulation([X(1:np),Y(1:np),Z(1:np)]);
coordinates = tri.Points; % all coordinate points
%element4 = tri.ConnectivityList; % all tetrahedrons
%e=tri.edges;
% all edges
eb=tri.freeBoundary;   % boundary faces
%tetramesh(tri,'FaceAlpha',0);
% figure();
% trisurf(eb,coordinates(:,1),coordinates(:,2),coordinates(:,3),1);
nt = size(tri,1); % no of triangles
fprintf(1,'Number of points = %d\n', np);
fprintf(1,'Number of triangles = %d\n', nt);

n_eb = size(eb,1); % no of boundary facets
ebd = zeros(size(eb));
ebn = zeros(size(eb));
n_ebd = 0;
n_ebn = 0;
tp = [1;1;1];
for j =1:n_eb
    t1 = sum(abs(X(eb(j,:))-tp));
    t2 = sum(abs(Y(eb(j,:))-tp));
    t3 = sum(abs(Z(eb(j,:))-tp));
    if(t1==0 || t2==0 ||t3 ==0)
        n_ebn = n_ebn+1;
        ebn(n_ebn,:) = eb(j,:);
    else
        n_ebd = n_ebd+1;
        ebd(n_ebd,:) = eb(j,:);
    end
end
ebn = ebn(1:n_ebn,:).';
ebd = ebd(1:n_ebd,:).';

fid=fopen('coordinates.dat','w');
fprintf(fid,'%24.14e %24.14e %24.14e\n',tri.Points');
fclose(fid);
fid=fopen('elements4.dat','w');
fprintf(fid,'%8d %8d %8d %8d\n',tri.ConnectivityList');
fclose(fid);
fid=fopen('dirichlet.dat','w');
fprintf(fid,'%8d %8d %8d\n', ebd);
fclose(fid);
fid=fopen('neumann.dat','w');
fprintf(fid,'%8d %8d %8d\n', ebn);
fclose(fid);
end
function value = u_d ( u )
value = u(:,1).^2+u(:,2).^2+u(:,3).^2;
end
function M = stima4(vertices)
% Given the vertices of a tetrahedron, evaluate the stiffness matrix element
G = [ones(1,4);vertices'] \ [zeros(1,3);eye(3)];
M = det([ones(1,4);vertices']) * G * G'/6;
end
function value = g ( u )
value = 2*ones(size(u,1),1);
end
function [coordinates,u]=fem3dlinear
% FEM3D three-dimensional finite element method for Laplacian.
% Initialization
load coordinates.dat;
load elements4.dat;
eval('load neumann.dat;','neumann=[];');
load dirichlet.dat;
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse(size(coordinates,1),size(coordinates,1));
b = sparse(size(coordinates,1),1);
% Assembly
for j = 1:size(elements4,1)
    A(elements4(j,:),elements4(j,:)) = A(elements4(j,:), ...
        elements4(j,:)) + stima4(coordinates(elements4(j,:),:));
end
% Volume Forces
for j = 1:size(elements4,1)
    b(elements4(j,:)) = b(elements4(j,:)) + ...
        det([1,1,1,1;coordinates(elements4(j,:),:)']) ...
        * f(sum(coordinates(elements4(j,:),:))/4) / 24;
end
% Neumann conditions
if ( ~isempty(neumann) )
    for j = 1 : size(neumann,1)
        b(neumann(j,:)) = b(neumann(j,:)) + ...
            norm(cross(coordinates(neumann(j,3),:)- ...
            coordinates(neumann(j,1),:), ...
            coordinates(neumann(j,2),:) - coordinates(neumann(j,1),:))) ...
            * g(sum(coordinates(neumann(j,:),:))/3)/6;
    end
end
% Dirichlet conditions
u = sparse(size(coordinates,1),1);
u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));
b = b - A * u;
% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);
% Graphic representation
figure();
trisurf([dirichlet;neumann],coordinates(:,1),coordinates(:,2),coordinates(:,3),full(u),'facecolor','interp');
colorbar;
end
function value = f ( u )
value = -6*ones(size(u,1),1);
end

