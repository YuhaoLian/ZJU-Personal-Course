function value = g ( u,t )

%*****************************************************************************80
%
% G evaluates the outward normal values assigned at Neumann boundary conditions.
%  Parameters:
%
%    Input, real U(N,M), contains the M-dimensional coordinates of N points.
%
%    Output, VALUE(N), contains the value of outward normal at each point
%    where a Neumann boundary condition is applied.
%
  value = 2*ones(size(u,1),1).*exp(-t);
end
