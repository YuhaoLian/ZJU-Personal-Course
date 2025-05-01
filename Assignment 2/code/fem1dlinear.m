function un=fem1dlinear(u0,u1,n,xn)
% Assemble the stiffness matrix
A = zeros(n);
% special treatment for the first and last nodes
A(1,1) = intlinear(xn(1),xn(2),1)+intlinear(xn(1),xn(2),3);
A(1,2) = intlinear(xn(1),xn(2),2)+intlinear(xn(1),xn(2),4);
A(2,1) = A(1,2);
A(n,n) = intlinear(xn(n-1),xn(n),1)+intlinear(xn(n-1),xn(n),3);
for k = 2:n-1
    A(k,k) = intlinear(xn(k-1),xn(k),1)+intlinear(xn(k),xn(k+1),1)+...
        intlinear(xn(k-1),xn(k),3)+intlinear(xn(k),xn(k+1),3);
    A(k,k+1) = intlinear(xn(k),xn(k+1),2)+intlinear(xn(k),xn(k+1),4);
    A(k+1,k) = A(k,k+1);
end
b = zeros(n-2,1); % the right hand side
b(1) = -A(2,1)*u0;
b(end) = -A(n-1,n)*u1;
un = A(2:n-1,2:n-1)\b; % solve the system
end