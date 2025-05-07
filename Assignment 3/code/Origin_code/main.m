errval = zeros(size(2:10));
for k = 2:10
    rectangle2d(k);
    [coordinates,u]=fem2dlinear;
    uex = coordinates(:,1).^2 - 1+coordinates(:,2).^2;
    errval(k-1)=norm(uex - u)/norm(uex);
end
plot(2:10,errval,'*-')
