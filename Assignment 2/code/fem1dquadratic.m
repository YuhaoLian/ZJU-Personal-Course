function un = fem1dquadratic(u0, u1, n, xn)
% Assemble the stiffness matrix A for quadratic elements
A = zeros(n, n);
% Reference matrices K_ref and M_ref for quadratic elements
K_ref = [7/6, -4/3, 1/6;
         -4/3, 8/3, -4/3;
         1/6, -4/3, 7/6];
M_ref = [4/15,  2/15, -1/15;
          2/15, 16/15, 2/15;
         -1/15, 2/15, 4/15];

num_elements = (n - 1)/2; % Number of quadratic elements

for k = 1:num_elements
    left = 2*k -1;
    mid = 2*k;
    right = 2*k +1;
    h = xn(right) - xn(left);
    % Compute local matrix contribution
    local_A = (2/h) * K_ref + (h/2) * M_ref;
    % Assemble into global matrix
    A(left, left) = A(left, left) + local_A(1,1);
    A(left, mid) = A(left, mid) + local_A(1,2);
    A(left, right) = A(left, right) + local_A(1,3);
    A(mid, left) = A(mid, left) + local_A(2,1);
    A(mid, mid) = A(mid, mid) + local_A(2,2);
    A(mid, right) = A(mid, right) + local_A(2,3);
    A(right, left) = A(right, left) + local_A(3,1);
    A(right, mid) = A(right, mid) + local_A(3,2);
    A(right, right) = A(right, right) + local_A(3,3);
end

% Apply boundary conditions
internal_nodes = 2:n-1;
% Construct system matrix and RHS
A_system = A(internal_nodes, internal_nodes);
b_system = -A(internal_nodes, 1)*u0 - A(internal_nodes, n)*u1;
% Solve the linear system
un_system = A_system \ b_system;
% Combine the solution
un = zeros(n,1);
un(1) = u0;
un(end) = u1;
un(internal_nodes) = un_system;
end