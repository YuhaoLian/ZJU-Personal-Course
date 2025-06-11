function [err] = error_calculate(u, u_exact)
coordinates = load('coordinates.dat');
dirichletNodes = load('dirichletNodes.dat');  % 边界节点
nNodes = size(coordinates, 1);

err_nodes = setdiff(1:nNodes, dirichletNodes);
num_free = length(err_nodes);

err_numerator = 0;
err_denominator = 0;

for i = 1:num_free
    node_idx = err_nodes(i);
    err_numerator = err_numerator + (u(node_idx) - u_exact(node_idx))^2;
    err_denominator = err_denominator + u_exact(node_idx)^2;
end

err = sqrt(err_numerator) / sqrt(err_denominator);
end

