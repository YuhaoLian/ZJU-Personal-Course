clc;
clear;

errval = zeros(9,1);
for k = 2:10
    rectangle2d(k);
    [coordinates,u] = fem2dlinear;
    uex = exp(coordinates(:,1) + coordinates(:,2));
    errval(k-1) = norm(uex - u)/norm(uex);
end

figure('Position', [100, 100, 800, 600]); % 设置图形窗口大小
semilogy(2:10, errval, 'o-', ...          % 使用对数坐标和带填充的圆点
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.3, 0.6, 0.9], ...
    'MarkerEdgeColor', 'k');

xlabel('Mesh Refinement Level', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Error', 'FontSize', 12, 'FontWeight', 'bold');
title('Finite Element Method Convergence', 'FontSize', 14, 'FontWeight', 'bold');

grid on;
grid minor; % 添加主次网格线
set(gca, 'FontSize', 11, ...             % 坐标轴字体设置
    'XColor', [0.3, 0.3, 0.3], ...       % 深灰色坐标轴
    'YColor', [0.3, 0.3, 0.3], ...
    'GridLineStyle', '--', ...           % 虚线网格
    'GridAlpha', 0.4);

xticks(2:10); % 明确设置x轴刻度
xlim([2, 10]); % 调整x轴范围

% 优化图形边距
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));