function [coordinates,u]=fem3dheatlineargif
% FEM3D three-dimensional finite element method for Laplacian.
% 初始化部分保持不变
fprintf ( 1, 'A program to demonstrate the 3D finite element method.\n' );
fprintf ( 1, 'Heat equation with backward Euler for time discretization.\n' );
load coordinates.dat;
load elements4.dat;
eval('load neumann.dat;','neumann=[];');
load dirichlet.dat;
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse(size(coordinates,1),size(coordinates,1));
B = sparse(size(coordinates,1),size(coordinates,1));
T = 1; dt = 0.01; N = T/dt; % 调整时间步长以提高计算效率
U = zeros(size(coordinates,1),N+1);
b = sparse(size(coordinates,1),1);

% 组装矩阵（保持不变）
for j = 1:size(elements4,1)
    A(elements4(j,:),elements4(j,:)) = A(elements4(j,:), ...
        elements4(j,:)) + stima4(coordinates(elements4(j,:),:));
end

for j = 1:size(elements4,1)
    B(elements4(j,:),elements4(j,:)) = B(elements4(j,:), ...
        elements4(j,:)) + det([1,1,1,1;coordinates(elements4(j,:),:)'])...
        *[2,1,1,1; 1,2,1,1; 1,1,2,1; 1,1,1,2]/120;
end

% 初始条件（保持不变）
U(:,1) = u_0(coordinates);

% 新增：动图参数设置
gif_filename = 'heat3d.gif';
frame_interval = 5;    % 每5帧捕获一次
delay_time = 0.2;      % 每帧延迟时间
caxis_range = [0 1.2]; % 固定颜色轴范围
view_angle = [30,30];  % 初始视角

% 初始化图形窗口
fig = figure('Position', [200 200 800 600]);
colormap(jet);

% 时间推进循环
for n = 2:N+1
    % 原有计算步骤保持不变
    b = sparse(size(coordinates,1),1);
    
    % Volume forces
    for j = 1:size(elements4,1)
        b(elements4(j,:)) = b(elements4(j,:)) + ...
            det([1,1,1,1; coordinates(elements4(j,:),:)']) * ...
            dt*f(sum(coordinates(elements4(j,:),:))/4,(n-1)*dt)/24;
    end
    
    % Neumann conditions
    for j = 1 : size(neumann,1)
        b(neumann(j,:)) = b(neumann(j,:)) + ...
            norm(cross(coordinates(neumann(j,3),:)- ...
            coordinates(neumann(j,1),:), ...
            coordinates(neumann(j,2),:) - coordinates(neumann(j,1),:))) ...
            * dt * g(sum(coordinates(neumann(j,:),:))/3, (n-1)*dt)/6;
    end
    
    % Previous timestep
    b = b + B * U(:,n-1);
    
    % Dirichlet conditions
    u = sparse(size(coordinates,1),1);
    u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:),n*dt);
    b = b - (dt * A + B) * u;
    
    % Computation of the solution
    u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
        B(FreeNodes,FreeNodes))\b(FreeNodes);
    U(:,n) = u;
    
    % 动态绘图部分
    if mod(n,frame_interval) == 0
        clf;
        trisurf([dirichlet;neumann],...
            coordinates(:,1),coordinates(:,2),coordinates(:,3),...
            full(u),'facecolor','interp');
        title(sprintf('3D Heat Equation Simulation\nTime = %.2f s',(n-1)*dt));
        colorbar;
        caxis(caxis_range);
        view(view_angle);
        view_angle = view_angle + [0.5, 0]; % 每帧稍微旋转视角
        drawnow;
        
        % 捕获帧并写入GIF
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == frame_interval
            imwrite(imind,cm,gif_filename,'gif',...
                'Loopcount',inf,'DelayTime',delay_time);
        else
            imwrite(imind,cm,gif_filename,'gif',...
                'WriteMode','append','DelayTime',delay_time);
        end
    end
end

% 最终静态图
figure('Position', [200 200 800 600]);
trisurf([dirichlet;neumann],...
    coordinates(:,1),coordinates(:,2),coordinates(:,3),...
    full(u),'facecolor','interp');
title('Final Temperature Distribution');
colorbar;
caxis(caxis_range);
view(3);
end