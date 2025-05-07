# 三维有限元方法实验报告

### —连宇昊

## 一、问题描述
求解三维椭圆型偏微分方程：
$$
\begin{cases}
-\Delta u + 3u = 0, & x \in \Omega \\
u = e^{x+y+z}, & x \in \partial\Omega
\end{cases}
$$
其中计算区域为立方体：
$$
\Omega = [-1,1]^3
$$
精确解为：
$$
u_{ex} = e^{x+y+z}
$$



## 二、实验原理
### 1. 弱形式推导
对原方程 $-Δu + 3u = 0$ 两边同乘测试函数 $v \in H_0^1(Ω)$，利用Green公式：

$$
\begin{aligned}
\int_Ω (-Δu)v \, dΩ + 3\int_Ω uv \, dΩ &= 0 \\
\int_Ω \nabla u \cdot \nabla v \, dΩ + 3\int_Ω uv \, dΩ &= \int_{∂Ω} \frac{∂u}{∂n}v \, ds
\end{aligned}
$$

因Dirichlet边界条件 $u|_{\partialΩ}=g$，取$v|_{\partialΩ}=0$，得到弱形式：

$$
a(u,v) = \int_Ω (\nabla u \cdot \nabla v + 3uv) \, dΩ = 0
$$

### 2. 有限元离散
采用线性单元基函数$\phi_i$，近似解$u_h = \sum_{j=1}^n u_j\phi_j$，得到代数方程组：

$$
(K + 3M)u = 0
$$

其中刚度矩阵$K$和质量矩阵$M$的元素为：
$$
K_{ij} = \int_Ω \nabla\phi_i \cdot \nabla\phi_j \, dΩ, \quad 
M_{ij} = \int_Ω \phi_i\phi_j \, dΩ
$$

## 三、方法实现

### 1. 有限元离散
```matlab
% 主程序框架
for k = 2:10
    rectangle3d(k);              % 生成三维网格
    [coordinates,u] = fem3dlinear;% 有限元求解
    uex = exp(sum(coordinates,2));% 计算精确解
    errval(k-1) = norm(uex - u)/norm(uex); % 计算相对误差
end
```

### 2. 刚度矩阵与质量矩阵

```matlab
function M = stima4(vertices)    % 刚度矩阵
d = det([ones(1,4); vertices']);
G = [ones(1,4); vertices'] \ [zeros(1,3); eye(3)];
M = abs(d) * G * G' / 6;

function M = stima4_mass(vertices) % 质量矩阵
d = det([ones(1,4); vertices']);
V = abs(d) / 6;
M = V * [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2] / 20;
```

### 3. 边界条件处理

```matlab
% 强制施加Dirichlet边界条件
dirichlet = tri.freeBoundary; % 获取自由边界
u(unique(dirichlet)) = exp(sum(coordinates(unique(dirichlet),:),2));
```

## 四、结果分析



![可视化](D:\Project\MATLAB\工业有限元\ZJU-Personal-Course\Assignment 4\image\可视化.png)

<center>图1 三维有限元可视化</center>



<img src="D:\Project\MATLAB\工业有限元\ZJU-Personal-Course\Assignment 4\image\误差分析.png" alt="误差分析" width="600" height="250">

<center>图2 三维有限元误差分析</center>

​	如上图所示，当网格密度从n=2增加到n=10时，相对误差不断降低，符合结果预期。