% 设置公共参数
a1 = 0; a2 = 1; % 区间 [0,1]
u0 = 0; u1 = 1; % 边界条件
m = 10; % 细化次数

% 初始化误差数组
err_linear = zeros(m, 1);
err_quadratic = zeros(m, 1);

% 线性基函数的误差计算
for k = 1:m
    n_linear = k*10 + 1; % 线性节点数
    xn_linear = linspace(a1, a2, 2*(n_linear - 1) + 1).';
    un_linear = fem1dlinear(u0, u1, 2*(n_linear - 1) + 1, xn_linear); % 假设返回内部节点解
    % 精确解在内部节点
    uex_linear = exp(1)/(exp(2)-1) * (exp(xn_linear(2:end-1)) - exp(-xn_linear(2:end-1)));
    err_linear(k) = norm(un_linear - uex_linear) / norm(uex_linear);
end

% 二次基函数的误差计算
for k = 1:m
    n_linear = k*10 + 1; % 对应的线性节点数
    n_quadratic = 2*(n_linear - 1) + 1; % 二次节点数
    xn_quadratic = linspace(a1, a2, n_quadratic).';
    un_quadratic = fem1dquadratic(u0, u1, n_quadratic, xn_quadratic); % 返回所有节点解
    % 精确解在内部节点（排除边界）
    uex_quadratic = exp(1)/(exp(2)-1) * (exp(xn_quadratic(2:end-1)) - exp(-xn_quadratic(2:end-1)));
    err_quadratic(k) = norm(un_quadratic(2:end-1) - uex_quadratic) / norm(uex_quadratic);
end

% 绘制对比图
figure;
plot((1:m)*20 + 1, err_linear, '-*r', 'DisplayName', 'Linear Basis', 'LineWidth', 2, 'MarkerSize', 11);
hold on;
plot((1:m)*20 + 1, err_quadratic, '-ob', 'DisplayName', 'Quadratic Basis', 'LineWidth', 2, 'MarkerSize', 11);
hold off;
title('FEM Error Comparison: Linear vs Quadratic Basis');
xlabel('Number of Nodes');
ylabel('Relative Error');
legend('show');
grid on;

% 使用对数坐标以更清晰显示收敛阶
figure;
semilogy((1:m)*20 + 1, err_linear, '-*r', 'DisplayName', 'Linear Basis', 'LineWidth', 2, 'MarkerSize', 11);
hold on;
semilogy((1:m)*20 + 1, err_quadratic, '-ob', 'DisplayName', 'Quadratic Basis', 'LineWidth', 2, 'MarkerSize', 11);
hold off;
title('FEM Error Comparison (Log Scale)');
xlabel('Number of Nodes');
ylabel('Relative Error (log scale)');
legend('show');
grid on;