clc;
clear;

%% 测试单个点误差 比较一次基函数和二次基函数

x = 1.113;
n = 10;
error1 = [];
error2 = [];
n_values = [];

threshold = 1E5;

while n <= threshold
    xn = linspace(0, 5, n);
    un = log(xn + 1);

    % 这里假设linearrep和Integral_liner_rep函数已经定义
    [u1,f1,f2] = linearrep(n,xn,un,x);
    [u2,f3,f4] = Integral_liner_rep(n,xn,un,x);

    u_true = log(x + 1);

    error1 = [error1; abs(u_true - u1)];
    error2 = [error2; abs(u_true - u2)];
    n_values = [n_values; n];
    
    n = n * 10;
end

% 绘制误差曲线，使用loglog函数让横轴和纵轴都为对数坐标，并加粗曲线
figure;
loglog(n_values, error1, 'b-o', 'DisplayName', '一次基函数误差', 'LineWidth', 2, 'MarkerSize', 11);
hold on;
loglog(n_values, error2, 'r-s', 'DisplayName', '二次基函数误差', 'LineWidth', 2, 'MarkerSize', 11);
hold off;

% 添加标题和标签
title('信噪比函数基函数有限元计算误差曲线');
xlabel('n', 'FontName', 'Times New Roman');
ylabel('误差');

% 添加图例
legend;

% 网格线
grid on;

%% 绘制对应区间的曲线 比较误差
Sample_num = 512;

x = linspace(0, 5, Sample_num);
u_true = log(x + 1);
u1 = zeros(1, Sample_num);
u2 = zeros(1, Sample_num);
phi11 = zeros(1, Sample_num);
phi12 = zeros(1, Sample_num);
phi21 = zeros(1, Sample_num);
phi22 = zeros(1, Sample_num);
phi23 = zeros(1, Sample_num);

n = 5; % 可调参数
xn = linspace(0, 5, n);
un = log(xn + 1);

for i = 1:length(x) 
    [u1(i),phi11(i),phi12(i)] = linearrep(n,xn,un,x(i));
    [u2(i),phi21(i),phi22(i),phi23(i)] = Integral_liner_rep(n,xn,un,x(i));
end

% 绘制三条曲线
figure;
plot(x, u_true, 'b-', 'DisplayName', '真实值', 'LineWidth', 2, 'MarkerSize', 11);
hold on;
plot(x, u1, 'r--', 'DisplayName', '一次基函数计算值', 'LineWidth', 2, 'MarkerSize', 11);
plot(x, u2, 'g-.', 'DisplayName', '二次基函数计算值', 'LineWidth', 2, 'MarkerSize', 11);
hold off;

% 添加标题和标签
title('不同方法计算结果对比曲线');
xlabel('x', 'FontName', 'Times New Roman');
ylabel('函数值');

% 添加图例
legend;

% 显示网格线
grid on;
    
% 绘制 phi11 和 phi12 曲线
figure;
plot(x, phi11, 'm-', 'DisplayName', 'phi11', 'LineWidth', 2, 'MarkerSize', 11);
hold on;
plot(x, phi12, 'c--', 'DisplayName', 'phi12', 'LineWidth', 2, 'MarkerSize', 11);
hold off;
title('phi11 和 phi12 曲线');
xlabel('phi', 'FontName', 'Times New Roman');
ylabel('函数值');
legend;
grid on;

% 绘制 phi21、phi22 和 phi23 曲线
figure;
plot(x(:,257:end), phi21(:,257:end), 'y-', 'DisplayName', 'phi21', 'LineWidth', 2, 'MarkerSize', 11);
hold on;
plot(x(:,257:end), phi22(:,257:end), 'k--', 'DisplayName', 'phi22', 'LineWidth', 2, 'MarkerSize', 11);
plot(x(:,257:end), phi23(:,257:end), '--', 'DisplayName', 'phi23', 'LineWidth', 2, 'MarkerSize', 11);
hold off;
title('phi21、phi22 和 phi23 曲线');
xlabel('phi', 'FontName', 'Times New Roman');
ylabel('函数值');
legend;
grid on;

