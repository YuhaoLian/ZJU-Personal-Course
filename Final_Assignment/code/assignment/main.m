clc;
clear;

% 初始化误差数组
N = 8; % 细化次数
errors = zeros(N, 1);

for k = 1:N
    fprintf('正在处理第 %d/%d 层细化...\n', k, N);
    % 生成网格
    trapezoid2d(k);
    % 计算数值解
    [coordinates,u]=fem2dlinear();
    % 计算精确解
    u_exact = exact_solution(coordinates);
    errors(k) = error_calculate(u, u_exact);
end


% 绘制误差曲线
figure('Position', [100, 100, 800, 500]); % 设置图形窗口位置和大小

% 自定义颜色和标记样式
h = plot(1:N, errors, 'o-', ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.4660, 0.6740, 0.1880], ... % 绿色填充
    'Color', [0.2000, 0.3500, 0.6100], ...           % 深蓝色线条
    'MarkerEdgeColor', 'k');                          % 黑色标记边框

% 设置坐标轴
grid on;
grid minor; % 添加次要网格线
set(gca, ...
    'LineWidth', 1.2, ...
    'FontSize', 12, ...
    'GridAlpha', 0.3, ...
    'MinorGridAlpha', 0.15);

ylabel('相对L2误差');

% 设置x轴
xlabel('细化次数');
xticks(1:N); % 设置x轴刻度为整数

% 添加标题和图例
title('误差随网格细化的收敛性分析', 'FontSize', 14, 'FontWeight', 'bold');
legend('数值误差', 'Location', 'northeast');

