clc;
clear;

errval = zeros(size(2:10));
for k = 2:10
    % 生成网格文件
    rectangle3d(k);
    % 有限元求解
    [coordinates,u] = fem3dlinear;
    % 精确解计算
    uex = coordinates(:,1).^2 + coordinates(:,2).^2 + coordinates(:,3).^2;
    % 误差评估（相对误差）
    errval(k-1) = norm(uex - u) / norm(uex);
end

% 绘图美化
figure('Position', [100 100 800 500]);  % 设置图窗位置和尺寸（宽800，高500）
plot(2:10, errval, '*-', ...
    'LineWidth', 1.5, ...          % 线条宽度
    'MarkerSize', 10, ...          % 标记大小
    'MarkerEdgeColor', [0 0.447 0.741], ...  % 标记边缘颜色（蓝色系）
    'MarkerFaceColor', [0.85 0.325 0.098], ... % 标记填充颜色（橙色系）
    'Color', [0.466 0.674 0.188]);       % 线条颜色（绿色系）

% 图表标题与坐标轴标签
title('有限元误差随网格细化参数k的变化', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('网格细化参数k', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('相对误差', 'FontSize', 12, 'FontWeight', 'normal');

% 坐标轴设置
xlim([1.8 10.2]);         % x轴范围略微扩展
xticks(2:10);             % 显式设置x轴刻度
ytickformat('%.2e');      % y轴采用科学计数法显示（保留2位小数）
grid on;                  % 主网格线
grid minor;               % 次网格线（更细腻的网格）
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.5);  % 网格样式（虚线，透明度50%）

% 图例与字体增强
legend('数值解误差', 'Location', 'northeast');  % 图例位置（右上角）
% set(findall(gcf,'Type','text'), 'FontName', 'Arial');  % 统一字体为Arial（更专业）
    