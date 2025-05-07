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

