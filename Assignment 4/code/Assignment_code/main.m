clc;
clear;

errval = zeros(size(2:10));
for k = 2:10
    rectangle3d(k);
    [coordinates, u] = fem3dlinear;
    uex = exp(coordinates(:,1) + coordinates(:,2) + coordinates(:,3));
    errval(k-1) = norm(full(u) - uex) / norm(uex);
end

% 创建1行2列的子图布局
figure('Position', [100 100 1200 500]);  % 宽度调整为1200以容纳两个子图

% 第一个子图：线性坐标
subplot(1, 2, 1);
plot(2:10, errval, '*-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', [0 0.447 0.741], ...
    'MarkerFaceColor', [0.85 0.325 0.098], ...
    'Color', [0.466 0.674 0.188]);
title('线性坐标：误差随k变化', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('网格细化参数k', 'FontSize', 12);
ylabel('相对误差', 'FontSize', 12);
xlim([1.8 10.2]);
xticks(2:10);
ytickformat('%.2e');  % 科学计数法显示
grid on;
grid minor;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.5);
legend('数值解误差', 'Location', 'northeast');

% 第二个子图：对数坐标
subplot(1, 2, 2);
semilogy(2:10, errval, '*-', ...  % 关键修改：对数坐标绘图函数
    'LineWidth', 1.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', [0 0.447 0.741], ...
    'MarkerFaceColor', [0.85 0.325 0.098], ...
    'Color', [0.466 0.674 0.188]);
title('对数坐标：误差随k变化', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('网格细化参数k', 'FontSize', 12);
ylabel('相对误差（对数坐标）', 'FontSize', 12);
xlim([1.8 10.2]);
xticks(2:10);
grid on;
grid minor;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.5);
legend('数值解误差', 'Location', 'northeast');
    