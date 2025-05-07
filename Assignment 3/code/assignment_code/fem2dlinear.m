function [coordinates,u] = fem2dlinear
fprintf(1, 'Solving -Δu + 2u = 0 with Dirichlet boundary u = e^(x+y)\n');
load coordinates.dat;
eval('load elements3.dat;', 'elements3=[];');
eval('load neumann.dat;', 'neumann=[];');
eval('load dirichlet.dat;', 'dirichlet=[];');

A = sparse(size(coordinates,1), size(coordinates,1));
b = sparse(size(coordinates,1), 1);

% Assemble stiffness matrix (Laplacian + 2*Mass)
for j = 1:size(elements3,1)
    verts = coordinates(elements3(j,:), :);
    A(elements3(j,:), elements3(j,:)) = A(elements3(j,:), elements3(j,:)) ...
        + stima3(verts) + 2 * stima_mass(verts);
end

% Apply Dirichlet conditions
u = sparse(size(coordinates,1), 1);
BoundNodes = unique(dirichlet);
u(BoundNodes) = u_d(coordinates(BoundNodes,:));
b = b - A * u;

% Solve
FreeNodes = setdiff(1:size(coordinates,1), BoundNodes);
u(FreeNodes) = A(FreeNodes, FreeNodes) \ b(FreeNodes);

% Plot solution
figure();
trisurf(elements3, coordinates(:,1), coordinates(:,2), full(u));
title('Computed Solution');
end