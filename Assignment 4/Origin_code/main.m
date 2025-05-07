clc;
clear;

errval = zeros(size(2:10));
for k = 2:10
    %   generate the mesh file
    rectangle3d(k);
    %   solve the equation by FEM
    [coordinates,u]=fem3dlinear;
    %   the exact solution
    uex = coordinates(:,1).^2+coordinates(:,2).^2+coordinates(:,3).^2;
    %   evaluate the error
    errval(k-1)=norm(uex - u)/norm(uex);
end
% plot the error
figure(); plot(2:10,errval,'*-')