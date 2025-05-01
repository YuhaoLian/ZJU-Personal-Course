% FEM code for 1D PDE: -d^2u/dx^2+u = 0
% using linear basis
a1 = 0; a2 = 1; % the interval [0,1]
u0 = 0; u1 = 1; % the boundary condition
m = 10; % No. of refinement times
err = zeros(m,1); % Error
for k = 1:m
    n = k*10+1; % No. of nodes
    xn = linspace(a1,a2,n).'; % partition the interval
    un = fem1dlinear(u0,u1,n,xn);
    % exact solution
    uex = exp(1)./(exp(2)-1).*(exp(xn(2:n-1))-exp(-xn(2:n-1)));
    % calculate the error
    err(k) = norm(un-uex)/norm(uex);
end
% plot the error
figure; plot((1:m)*10+1, err, '-*');
title('Error of FEM in 1D');
xlabel('No. of nodes');
ylabel('Error');