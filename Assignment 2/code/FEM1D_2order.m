% FEM code for 1D PDE: -d^2u/dx^2+u = 0 using quadratic basis
a1 = 0; a2 = 1; % interval [0,1]
u0 = 0; u1 = 1; % boundary conditions
m = 10; % No. of refinement times
err = zeros(m,1); % Error

for k = 1:m
    n_linear = k*10 +1; % Linear nodes
    n = 2*(n_linear -1)+1; % Quadratic nodes
    xn = linspace(a1,a2,n).'; 
    un = fem1dquadratic(u0,u1,n,xn);
    % Exact solution at internal nodes
    uex = exp(1)./(exp(2)-1).*(exp(xn(2:n-1)) - exp(-xn(2:n-1)));
    % Calculate error
    err(k) = norm(un(2:n-1)-uex)/norm(uex);
end

% Plot the error
figure; plot((1:m)*20+1, err, '-*');
title('Error of FEM in 1D with Quadratic Elements');
xlabel('No. of nodes');
ylabel('Error');