function [u,phi1,phi2,phi3] = Integral_liner_rep(n,xn,un,x)
% 假设 xn(1) <= x <= xn(n)
k = max(find(x >= xn)); % 寻找到相关子区间
if(k == n) k = k - 2; end 
if(k == n - 1) k = k - 1; end % 处理边界条件 
Xevion = (x - xn(k)) / (xn(k + 2) - xn(k)); % 归一化到 [0,1]
phi1 = 2 * (Xevion - 1) * (Xevion - 0.5); % 二次基函数 phi_1
phi2 = 4 * Xevion * (1 - Xevion); % 二次基函数 phi_2
phi3 = 2 * Xevion * (Xevion - 0.5); % 二次基函数 phi_3    
u = un(k) * phi1 + un(k + 1) * phi2 + un(k + 2) * phi3;
end

