clc;
clear;


rectangle3d(10); % 生成10分网格
[coordinates,u] = fem3dheatlineargif;


errval = zeros(size(2:10));
for k = 2:10
    % 生成网格文件
    rectangle3d(k);
    % 有限元求解
    [coordinates,u] = fem3dheatlinear;
    % 精确解计算
    uex = (coordinates(:,1).^2 + coordinates(:,2).^2 + coordinates(:,3).^2).*exp(-1);
    % 误差评估（相对误差）
    errval(k-1) = norm(uex - u) / norm(uex);
end

% 创建带两个子图的图形窗口
figure('Position', [100 100 1200 500]);  % 加宽窗口显示两个子图

% ================= 子图1：线性坐标 =================
subplot(1,2,1);
plot(2:10, errval, '*-', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', [0 0.447 0.741], ...
    'MarkerFaceColor', [0.85 0.325 0.098], ...
    'Color', [0.466 0.674 0.188]);

title('线性坐标下的误差趋势', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('网格细化参数k', 'FontSize', 12);
ylabel('相对误差', 'FontSize', 12);
xlim([1.8 10.2]);
xticks(2:10);
grid on;
grid minor;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.5);

% ================= 子图2：对数坐标 =================
subplot(1,2,2);
semilogy(2:10, errval, '*-', ...  % 使用semilogy函数
    'LineWidth', 1.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', [0 0.447 0.741], ...
    'MarkerFaceColor', [0.85 0.325 0.098], ...
    'Color', [0.466 0.674 0.188]);

title('对数坐标下的误差趋势', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('网格细化参数k', 'FontSize', 12);
ylabel('相对误差（对数坐标）', 'FontSize', 12);
xlim([1.8 10.2]);
xticks(2:10);
yticks(10.^(-8:1:0));  % 设置对数刻度范围
ytickformat('10^{%.0f}');  % 指数格式显示
grid on;
grid minor;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.5, ...
    'YMinorTick', 'on', 'YMinorGrid', 'on');  % 启用次刻度

% ================= 公共图例设置 =================
legend('数值解误差', 'Location', 'northeast');