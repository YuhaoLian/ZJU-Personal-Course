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
