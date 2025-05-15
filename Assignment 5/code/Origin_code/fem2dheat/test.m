clc;
clear;

errval = zeros(size(2:10));
for k = 2:10
    %   generate the mesh file
    rectangle2d(k);
    %   solve the equation by FEM
    [coordinates,u]=fem2dheatlinear;
    %    [coordinates,u]=fem2dheatlinearfd; % non-stable version
    %   the exact solution
    uex = (coordinates(:,1).^2-1+coordinates(:,2).^2).*exp(-1);
    %   evaluate the error
    errval(k-1)=norm(uex-u)/norm(uex);
end
% plot the error
figure(); plot(2:10,errval,'*-')

function value = u_d( u,t )
value = (u(:,1).^2-1+u(:,2).^2).*exp(-t);
end


function value = u_0( u )
value = u(:,1).^2-1+u(:,2).^2;
end

function M = stima3 ( vertices )
G = [ ones(1,3); vertices' ] \ [ zeros(1,2); eye(2) ];
M = det ( [ ones(1,3); vertices' ] ) * G * G' / 2;
end

function rectangle2d(n)
%n=10;
xmin=-1;xmax=1; ymin=0;ymax=1;
nx=2*n-1; ny=n;
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
np=0; % no of vertices
X=zeros(2*nx*ny,1);
Y=zeros(2*nx*ny,1);
for k=1:ny
    for j=1:nx
        np  = np + 1;
        X(np) = x(j);
        Y(np) = y(k);
    end
end
for k=1:ny-1
    for j=1:nx-1
        np  = np + 1;
        X(np) = (x(j)+x(j+1))/2;
        Y(np) = (y(k)+y(k+1))/2;
    end
end
tri=delaunayTriangulation([X(1:np),Y(1:np)]);
%triplot(tri)
e=tri.edges;           % all edges
eb=tri.freeBoundary;   % boundary edges
nt   = size(tri,1); % no of triangles
fprintf(1,'Number of points         = %d\n', np);
fprintf(1,'Number of triangles      = %d\n', nt);

n_e  = size(e ,1); % no of all edges
n_eb = size(eb,1); % no of boundary edges
ebd = zeros(size(eb));
ebn = zeros(size(eb));
n_ebd = 0;
n_ebn = 0;
for j =1:n_eb
    yp1 = Y(eb(j,1));
    yp2 = Y(eb(j,2));
    if(yp1==1 && yp2==1)
        n_ebn = n_ebn+1;
        ebn(n_ebn,1) = eb(j,1);
        ebn(n_ebn,2) = eb(j,2);
    else
        n_ebd = n_ebd+1;
        ebd(n_ebd,1) = eb(j,1);
        ebd(n_ebd,2) = eb(j,2);
    end
end

ebn = ebn(1:n_ebn,:);
ebd = ebd(1:n_ebd,:);

fid=fopen('coordinates.dat','w');
fprintf(fid,'%24.14e %24.14e\n', tri.Points');
fclose(fid);

fid=fopen('elements3.dat','w');
fprintf(fid,'%8d %8d %8d\n', tri.ConnectivityList');
fclose(fid);

fid=fopen('dirichlet.dat','w');
fprintf(fid,'%8d %8d\n', ebd');
fclose(fid);

fid=fopen('neumann.dat','w');
fprintf(fid,'%8d %8d\n', ebn');
fclose(fid);
end

function value = g ( u,t )
value = 2*ones(size(u,1),1).*exp(-t);
end

function [coordinates,u]=fem2dheatlinear()
fprintf ( 1, 'A program to demonstrate the 2D finite element method.\n' );
fprintf ( 1, 'Heat equation with backward Euler for time discretization.\n' );
%  Read the nodal coordinate data file.
load coordinates.dat;
%  Read the triangular element data file.
eval ( 'load elements3.dat;', 'elements3=[];' );
%  Read the Neumann boundary condition data file.
eval ( 'load neumann.dat;', 'neumann=[];' );
%  Read the Dirichlet boundary condition data file.
eval ( 'load dirichlet.dat;', 'dirichlet=[];' );
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse ( size(coordinates,1), size(coordinates,1) );
B = sparse(size(coordinates,1),size(coordinates,1));
T = 1; dt = 0.01; N = T/dt;
U = zeros(size(coordinates,1),N+1);
b = sparse ( size(coordinates,1), 1 );
%  Assembly.
for j = 1 : size(elements3,1)
    A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
        + stima3(coordinates(elements3(j,:),:));
end
for j = 1:size(elements3,1)
    B(elements3(j,:),elements3(j,:)) = B(elements3(j,:), ...
        elements3(j,:)) + det([1,1,1;coordinates(elements3(j,:),:)'])...
        *[2,1,1;1,2,1;1,1,2]/24;
end
% Initial Condition
U(:,1) = u_0(coordinates);
% time steps
for n = 2:N+1
    b = sparse(size(coordinates,1),1);
    % Volume Forces
    for j = 1:size(elements3,1)
        b(elements3(j,:)) = b(elements3(j,:)) + ...
            det([1,1,1; coordinates(elements3(j,:),:)']) * ...
            dt*f(sum(coordinates(elements3(j,:),:))/3,(n-1)*dt)/6;
    end
    % Neumann conditions
    for j = 1 : size(neumann,1)
        b(neumann(j,:)) = b(neumann(j,:)) + ...
            norm(coordinates(neumann(j,1),:)-coordinates(neumann(j,2),:))*...
            dt*g(sum(coordinates(neumann(j,:),:))/2,(n-1)*dt)/2;
    end
    % previous timestep
    b = b + B * U(:,n-1);
    % Dirichlet conditions
    u = sparse(size(coordinates,1),1);
    u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:),n*dt);
    b = b - (dt * A + B) * u;
    % Computation of the solution
    u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
        B(FreeNodes,FreeNodes))\b(FreeNodes);
    U(:,n) = u;
end
%  Graphic representation.
figure();trisurf ( elements3, coordinates(:,1), coordinates(:,2), full(U(:,N+1)));
end

function value = f (u,t)
value = -(u(:,1).^2-1+u(:,2).^2).*exp(-t)-4*ones(size(u,1),1).*exp(-t);
end
