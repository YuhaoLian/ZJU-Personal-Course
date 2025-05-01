function [u,phi1,phi2] = linearrep(n,xn,un,x)
% Assume xn(1)<=x<=xn(n). Find the linear interpolation of u at x using linear basis function
k=max(find(x>=xn)); % Find the corresponding subinterval
if(k==n) k=k-1; end % In case x=xn(n)
xi = (x-xn(k))/(xn(k+1)-xn(k)); % Normalize to [0,1]
phi1 = 1-xi; % linear basis function phi_1
phi2 = xi; % linear basis function phi_2
u = un(k)*phi1+un(k+1)*phi2;
end
