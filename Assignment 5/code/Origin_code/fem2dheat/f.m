function value = f (u,t)

%*****************************************************************************80
%
%% F evaluates the right hand side of heat equation.
%
%  Discussion:
%
%    This routine must be changed by the user to reflect a particular problem.
%
%  Parameters:
%
%    Input, real U(N,M), contains the M-dimensional coordinates of N points.
%
%    Output, VALUE(N), contains the value of the right hand side of Laplace's
%    equation at each of the points.
%
%  n = size ( u, 1 );

%  value(1:n) = 2.0 * pi * pi * sin ( pi * u(1:n,1) ) .* sin ( pi * u(1:n,2) );
  value = -(u(:,1).^2-1+u(:,2).^2).*exp(-t)-4*ones(size(u,1),1).*exp(-t);
end
