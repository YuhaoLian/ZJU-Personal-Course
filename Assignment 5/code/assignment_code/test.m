clc;
clear;

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

function value = u_d ( u,t )
value = (u(:,1).^2+u(:,2).^2+u(:,3).^2).*exp(-t);
end
function value = u_0( u )
value = u(:,1).^2+u(:,2).^2+u(:,3).^2;
end
function M = stima4(vertices)
% Given the vertices of a tetrahedron, evaluate the stiffness matrix element
G = [ones(1,4);vertices'] \ [zeros(1,3);eye(3)];
M = det([ones(1,4);vertices']) * G * G'/6;
end
function rectangle3d(n)
%Generate the mesh data file for a cuboid
%n=10;
xmin=0;xmax=1;
ymin=0;ymax=1;zmin=0;zmax=1;
nx=n; ny=n; nz=n;
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
z=linspace(zmin,zmax,nz);
np=0; % no. of vertices
X=zeros(2*nx*ny*nz,1);
Y=zeros(2*nx*ny*nz,1);
Z=zeros(2*nx*ny*nz,1);
for l=1:nz
    for k=1:ny
        for j=1:nx
            np = np+ 1;
            X(np) = x(j);
            Y(np) = y(k);
            Z(np) = z(l);
        end
    end
end
for l=1:nz-1
    for k=1:ny-1
        for j=1:nx-1
            np = np+ 1;
            X(np) = (x(j)+x(j+1))/2;
            Y(np) = (y(k)+y(k+1))/2;
            Z(np) = (z(l)+z(l+1))/2;
        end
    end
end
tri=delaunayTriangulation([X(1:np),Y(1:np),Z(1:np)]);
coordinates = tri.Points; % all coordinate points
%element4 = tri.ConnectivityList; % all tetrahedrons
%e=tri.edges;
% all edges
eb=tri.freeBoundary;   % boundary faces
%tetramesh(tri,'FaceAlpha',0);

%%  绘图
% figure();
% trisurf(eb,coordinates(:,1),coordinates(:,2),coordinates(:,3),1);

%%
nt = size(tri,1); % no of triangles
fprintf(1,'Number of points = %d\n', np);
fprintf(1,'Number of triangles = %d\n', nt);

n_eb = size(eb,1); % no of boundary facets
ebd = zeros(size(eb));
ebn = zeros(size(eb));
n_ebd = 0;
n_ebn = 0;
tp = [1;1;1];
for j =1:n_eb
    t1 = sum(abs(X(eb(j,:))-tp));
    t2 = sum(abs(Y(eb(j,:))-tp));
    t3 = sum(abs(Z(eb(j,:))-tp));
    if(t1==0 || t2==0 ||t3 ==0)
        n_ebn = n_ebn+1;
        ebn(n_ebn,:) = eb(j,:);
    else
        n_ebd = n_ebd+1;
        ebd(n_ebd,:) = eb(j,:);
    end
end
ebn = ebn(1:n_ebn,:).';
ebd = ebd(1:n_ebd,:).';

fid=fopen('coordinates.dat','w');
fprintf(fid,'%24.14e %24.14e %24.14e\n',tri.Points');
fclose(fid);
fid=fopen('elements4.dat','w');
fprintf(fid,'%8d %8d %8d %8d\n',tri.ConnectivityList');
fclose(fid);
fid=fopen('dirichlet.dat','w');
fprintf(fid,'%8d %8d %8d\n', ebd);
fclose(fid);
fid=fopen('neumann.dat','w');
fprintf(fid,'%8d %8d %8d\n', ebn);
fclose(fid);
end
function value = g ( u,t )
value = 2*ones(size(u,1),1).*exp(-t);
end
function [coordinates,u]=fem3dheatlinear
% FEM3D three-dimensional finite element method for Laplacian.
% Initialization
fprintf ( 1, 'A program to demonstrate the 3D finite element method.\n' );
fprintf ( 1, 'Heat equation with backward Euler for time discretization.\n' );
load coordinates.dat;
load elements4.dat;
eval('load neumann.dat;','neumann=[];');
load dirichlet.dat;
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse(size(coordinates,1),size(coordinates,1));
B = sparse(size(coordinates,1),size(coordinates,1)); %% 不一样 时间导数项的离散化
T = 1; dt = 0.001; N = T/dt; %% 不一样
U = zeros(size(coordinates,1),N+1); %% 不一样 解向量集合
b = sparse(size(coordinates,1),1);

% Assembly
for j = 1:size(elements4,1)
    A(elements4(j,:),elements4(j,:)) = A(elements4(j,:), ...
        elements4(j,:)) + stima4(coordinates(elements4(j,:),:));
end

for j = 1:size(elements4,1)
    B(elements4(j,:),elements4(j,:)) = B(elements4(j,:), ...
        elements4(j,:)) + det([1,1,1,1;coordinates(elements4(j,:),:)'])...
        *[2,1,1,1; 1,2,1,1; 1,1,2,1; 1,1,1,2]/120;
end

% Initial Condition
U(:,1) = u_0(coordinates);
% time steps
for n = 2:N+1
    b = sparse(size(coordinates,1),1);
    % Volume Forces
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
    % previous timestep
    b = b + B * U(:,n-1);
    % Dirichlet conditions
    u = sparse(size(coordinates,1),1);
    u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:),n*dt);
    b = b - (dt * A + B) * u;
    % Computation of the solution
    u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
        B(FreeNodes,FreeNodes))\b(FreeNodes);
    U(:,n) = u;
end
%  Graphic representation.
figure();
trisurf([dirichlet;neumann],coordinates(:,1),coordinates(:,2),coordinates(:,3),full(u),'facecolor','interp');
colorbar;
end
function value = f ( u,t )
value = -(u(:,1).^2+u(:,2).^2+u(:,3).^2).*exp(-t)-6*ones(size(u,1),1).*exp(-t);
end
