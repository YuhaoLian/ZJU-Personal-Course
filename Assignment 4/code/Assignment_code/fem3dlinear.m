function [coordinates, u] = fem3dlinear
load coordinates.dat;
load elements4.dat;
eval('load neumann.dat;', 'neumann=[];');
load dirichlet.dat;

FreeNodes = setdiff(1:size(coordinates,1), unique(dirichlet));
A = sparse(size(coordinates,1), size(coordinates,1));
b = sparse(size(coordinates,1), 1);

% Assemble stiffness and mass matrices
for j = 1:size(elements4,1)
    K = stima4(coordinates(elements4(j,:), :));
    M = stima4_mass(coordinates(elements4(j,:), :));
    A(elements4(j,:), elements4(j,:)) = A(elements4(j,:), elements4(j,:)) + K + 3*M;
end

% Apply Dirichlet boundary conditions
u = zeros(size(coordinates,1), 1);
u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet), :));
b = b - A * u;

% Solve linear system
u(FreeNodes) = A(FreeNodes, FreeNodes) \ b(FreeNodes);

% Visualization
figure();
trisurf(dirichlet, coordinates(:,1), coordinates(:,2), coordinates(:,3), u, 'FaceColor', 'interp');
colorbar;
end