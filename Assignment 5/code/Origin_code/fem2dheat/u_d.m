function value = u_d( u,t )

%*****************************************************************************80
%
% U_D evaluates the Dirichlet boundary conditions.
%  Parameters:
%
%    Input, real U(N,M), contains the M-dimensional coordinates of N points.
%
%    Output, VALUE(N), contains the value of the Dirichlet boundary
%    condition at each point.
%
    value = (u(:,1).^2-1+u(:,2).^2).*exp(-t);
end
