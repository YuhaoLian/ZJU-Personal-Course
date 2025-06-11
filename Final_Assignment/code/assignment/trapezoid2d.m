function trapezoid2d(k)
% 四边形顶点定义
P1 = [-2, 0];   % 左下
P2 = [ 2, 0];   % 右下
P3 = [ 1, 1];   % 右上
P4 = [-1, 1];   % 左上

m = 2^k;        % 每边分割数
n_xi = m + 1;   % xi方向节点数
n_eta = m + 1;  % eta方向节点数

% 生成参数域网格
xi = linspace(0, 1, n_xi);
eta = linspace(0, 1, n_eta);

% 坐标生成(双线性插值)
np = n_xi * n_eta;
X = zeros(np,1);
Y = zeros(np,1);
xi_param = zeros(np,1);
eta_param = zeros(np,1);

idx = 1;
for j = 1:n_eta
    for i = 1:n_xi
        xi_i = xi(i);
        eta_j = eta(j);
        
        % 双线性插值
        X(idx) = (1-xi_i)*(1-eta_j)*P1(1) + xi_i*(1-eta_j)*P2(1) + ...
            xi_i*eta_j*P3(1) + (1-xi_i)*eta_j*P4(1);
        Y(idx) = (1-xi_i)*(1-eta_j)*P1(2) + xi_i*(1-eta_j)*P2(2) + ...
            xi_i*eta_j*P3(2) + (1-xi_i)*eta_j*P4(2);
        
        % 存储参数坐标
        xi_param(idx) = xi_i;
        eta_param(idx) = eta_j;
        
        idx = idx + 1;
    end
end

% 生成三角形连接性
triangles = zeros(2*m*m, 3);
tri_idx = 1;
for j = 1:m
    for i = 1:m
        ll = (j-1)*n_xi + i;
        lr = ll + 1;
        ul = ll + n_xi;
        ur = ul + 1;
        
        triangles(tri_idx,:) = [ll, lr, ul];
        triangles(tri_idx+1,:) = [lr, ur, ul];
        tri_idx = tri_idx + 2;
    end
end

% 创建三角剖分对象
t = triangulation(triangles, X, Y);

% 获取所有边界边
edges = t.edges();
tol = 1e-10;

% 生成唯一边并计算中点
unique_edges = containers.Map('KeyType','char','ValueType','any');
for edge_idx = 1:size(edges, 1)
    p1 = edges(edge_idx, 1);
    p2 = edges(edge_idx, 2);
    if p1 > p2
        [p1, p2] = deal(p2, p1);
    end
    key = sprintf('%d-%d', p1, p2);
    if ~isKey(unique_edges, key)
        unique_edges(key) = struct('p1', p1, 'p2', p2);
    end
end

% 计算中点并扩展顶点列表
num_edges = size(edges, 1);
X_new = [X; zeros(num_edges, 1)];
Y_new = [Y; zeros(num_edges, 1)];
xi_new = [xi_param; zeros(num_edges, 1)];
eta_new = [eta_param; zeros(num_edges, 1)];
midpoint_map = containers.Map('KeyType','char','ValueType','double');

current_midpoint_index = np + 1;
for key = unique_edges.keys()
    key_str = key{1};
    data = unique_edges(key_str);
    p1 = data.p1;
    p2 = data.p2;
    
    % 计算中点坐标
    mid_x = (X(p1) + X(p2)) / 2;
    mid_y = (Y(p1) + Y(p2)) / 2;
    
    % 计算中点参数坐标
    mid_xi = (xi_param(p1) + xi_param(p2)) / 2;
    mid_eta = (eta_param(p1) + eta_param(p2)) / 2;
    
    % 存储中点
    midpoint_map(key_str) = current_midpoint_index;
    X_new(current_midpoint_index) = mid_x;
    Y_new(current_midpoint_index) = mid_y;
    xi_new(current_midpoint_index) = mid_xi;
    eta_new(current_midpoint_index) = mid_eta;
    
    current_midpoint_index = current_midpoint_index + 1;
end

% 构建六节点三角形单元
triangles_new = zeros(size(triangles, 1), 6);
for tri_idx = 1:size(triangles, 1)
    v1 = triangles(tri_idx, 1);
    v2 = triangles(tri_idx, 2);
    v3 = triangles(tri_idx, 3);
    
    % 获取各边中点
    m12 = midpoint_map(sprintf('%d-%d', min(v1,v2), max(v1,v2)));
    m23 = midpoint_map(sprintf('%d-%d', min(v2,v3), max(v2,v3)));
    m31 = midpoint_map(sprintf('%d-%d', min(v3,v1), max(v3,v1)));
    
    % 按逆时针顺序排列
    triangles_new(tri_idx, :) = [v1, m12, v2, m23, v3, m31];
end

% 确定边界节点
boundaryNodes = find((abs(xi_new) < tol) | (abs(xi_new - 1) < tol) | ...
    (abs(eta_new) < tol) | (abs(eta_new - 1) < tol));


figure;
triplot(t);
title(sprintf('Quadrilateral Mesh with %d Triangles', size(triangles,1)));


% 保存文件
fid = fopen('coordinates.dat','w');
fprintf(fid, '%24.16e %24.16e\n', [X_new, Y_new]');
fclose(fid);

fid = fopen('elements6.dat','w');
fprintf(fid, '%6d %6d %6d %6d %6d %6d\n', triangles_new');
fclose(fid);

fid = fopen('dirichletNodes.dat','w');
fprintf(fid, '%6d\n', boundaryNodes);
fclose(fid);

fprintf('细化级别 %d: 节点数 = %d, 单元数 = %d\n', k, length(X_new), size(triangles_new,1));
end