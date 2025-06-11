function [coordinates,u]=fem2dlinear
% 加载网格数据
coordinates = load('coordinates.dat');
elements6 = load('elements6.dat');
dirichletNodes = load('dirichletNodes.dat');  % 边界节点

% 总节点数
nNodes = size(coordinates, 1);

% 初始化全局刚度矩阵和载荷向量
A = sparse(nNodes, nNodes);
b = sparse(nNodes, 1);

% 遍历所有单元
nElements = size(elements6, 1);
for el = 1:nElements
    % 获取单元节点
    nodes = elements6(el, :);
    verts = coordinates(nodes, :);
    
    % 计算单元刚度矩阵和载荷向量
    [Ke, Fe] = elemStiffness(verts);
    
    % 组装到全局矩阵
    A(nodes, nodes) = A(nodes, nodes) + Ke;
    b(nodes) = b(nodes) + Fe;
end

u = zeros(nNodes, 1);
u_exact = exact_solution(coordinates);
u(dirichletNodes) = u_exact(dirichletNodes);  % 直接设置边界值

% 分离自由节点和边界节点
freeNodes = setdiff(1:nNodes, dirichletNodes);

% 修正右端项
b(freeNodes) = b(freeNodes) - A(freeNodes, dirichletNodes) * u(dirichletNodes);

% 求解自由节点上的解
u(freeNodes) = A(freeNodes, freeNodes) \ b(freeNodes);

% %  Graphic representation.
figure();
trisurf ( elements6, coordinates(:,1), coordinates(:,2), full(u) );
end
