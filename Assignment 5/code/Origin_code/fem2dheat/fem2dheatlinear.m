function [coordinates,u]=fem2dheatlinear
% FEM two-dimensional finite element method for 2D Heat equation.
% Contact: laijun6@zju.edu.cn

fprintf ( 1, 'A program to demonstrate the 2D finite element method.\n' );
fprintf ( 1, 'Heat equation with backward Euler for time discretization.\n' );
%  Read the nodal coordinate data file.
load coordinates.dat;
%  Read the triangular element data file.
eval ( 'load elements3.dat;', 'elements3=[];' );
%  Read the Neumann boundary condition data file.
eval ( 'load neumann.dat;', 'neumann=[];' );
%  Read the Dirichlet boundary condition data file.
eval ( 'load dirichlet.dat;', 'dirichlet=[];' );
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse ( size(coordinates,1), size(coordinates,1) );
B = sparse(size(coordinates,1),size(coordinates,1)); %% 不一样 时间导数项的离散化
T = 1; dt = 0.01; N = T/dt; %% 不一样
U = zeros(size(coordinates,1),N+1); %% 不一样 解向量集合
% b = sparse ( size(coordinates,1), 1 );
%  Assembly.
for j = 1 : size(elements3,1)
    A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
        + stima3(coordinates(elements3(j,:),:));
end
%% 不一样
for j = 1:size(elements3,1)
    B(elements3(j,:),elements3(j,:)) = B(elements3(j,:), ...
        elements3(j,:)) + det([1,1,1;coordinates(elements3(j,:),:)'])...
        *[2,1,1;1,2,1;1,1,2]/24;
end
%% 不一样
% Initial Condition
U(:,1) = u_0(coordinates);
% time steps
for n = 2:N+1
    b = sparse(size(coordinates,1),1);
    % Volume Forces
    for j = 1:size(elements3,1)
        b(elements3(j,:)) = b(elements3(j,:)) + ...
            det([1,1,1; coordinates(elements3(j,:),:)']) * ...
            dt*f(sum(coordinates(elements3(j,:),:))/3,(n-1)*dt)/6;
    end
    % Neumann conditions
    for j = 1 : size(neumann,1)
        b(neumann(j,:)) = b(neumann(j,:)) + ...
            norm(coordinates(neumann(j,1),:)-coordinates(neumann(j,2),:))*...
            dt*g(sum(coordinates(neumann(j,:),:))/2,(n-1)*dt)/2;
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
figure();trisurf ( elements3, coordinates(:,1), coordinates(:,2), full(U(:,N+1)));