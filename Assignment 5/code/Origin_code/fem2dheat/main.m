% FEM two-dimensional finite element method for 2D Heat equation.
% Contact: laijun6@zju.edu.cn

errval = zeros(size(2:10));
for k = 2:10
%   generate the mesh file    
    rectangle2d(k);
%   solve the equation by FEM    
    [coordinates,u]=fem2dheatlinear;
%    [coordinates,u]=fem2dheatlinearfd; % non-stable version
%   the exact solution    
    uex = (coordinates(:,1).^2-1+coordinates(:,2).^2).*exp(-1);
%   evaluate the error    
    errval(k-1)=norm(uex-u)/norm(uex);
end
% plot the error
figure(); plot(2:10,errval,'*-')
