function [coordinates, u] = fem3dlinear
load coordinates.dat;
load elements4.dat;
load dirichlet.dat;
neumann = [];

FreeNodes = setdiff(1:size(coordinates,1), unique(dirichlet));
A = sparse(size(coordinates,1), size(coordinates,1));
b = sparse(size(coordinates,1), 1);

% 组装刚度矩阵 + 3*质量矩阵
for j = 1:size(elements4, 1)
    verts = coordinates(elements4(j,:), :);
    stiff = stima4(verts);
    mass = stima4_mass(verts);
    A(elements4(j,:), elements4(j,:)) = A(elements4(j,:), elements4(j,:)) + stiff + 3 * mass;
end

% Dirichlet边界处理
u = sparse(size(coordinates,1), 1);
dirichlet_nodes = unique(dirichlet);
u(dirichlet_nodes) = u_d(coordinates(dirichlet_nodes, :));
b = b - A * u;

% 求解
u(FreeNodes) = A(FreeNodes, FreeNodes) \ b(FreeNodes);

% 绘图
figure();
trisurf(dirichlet, coordinates(:,1), coordinates(:,2), coordinates(:,3), full(u), 'FaceColor', 'interp');
colorbar;
end