# 三维热传导方程有限元方法
**——连宇昊**

## 一、实验目的
1. 完成作业中三维热传导方程的有限元离散方法
2. 实现基于后向欧拉法的时间离散方案
3. 分析网格细化对数值解精度的影响

## 二、实验原理

### 1. 控制方程与定解条件
考虑三维长方体区域$\Omega=(0,1)^3$内的热传导方程：
```math
\begin{cases}
\frac{\partial u}{\partial t} - \Delta u = f(x,t) & \text{在}\ \Omega\times(0,T]\\
u(x,0) = u_0(x) & \text{初值条件}\\
u = g_D & \text{在}\ \Gamma_D\times(0,T] \ (\text{Dirichlet边界})\\
\frac{\partial u}{\partial n} = g_N & \text{在}\ \Gamma_N\times(0,T] \ (\text{Neumann边界})
\end{cases}
```

### 2. 弱形式推导
乘以测试函数$v\in H^1_0(\Omega)$并积分：
$$
\int_\Omega \frac{\partial u}{\partial t}v \,d\mathbf{x} + \int_\Omega \nabla u \cdot \nabla v \,d\mathbf{x} = \int_\Omega fv \,d\mathbf{x} + \int_{\Gamma_N} g_N v \,ds
$$

### 3. 时空离散格式
- ​**空间离散**：采用P1四面体单元
- ​**时间离散**：后向欧拉法

离散格式：
$$
\left(\frac{1}{\Delta t}\mathbf{M} + \mathbf{A}\right)\mathbf{U}^{n} = \frac{1}{\Delta t}\mathbf{M}\mathbf{U}^{n-1} + \mathbf{F}^n
$$

其中质量矩阵和刚度矩阵：
$$
\mathbf{M} = \begin{bmatrix}
\langle \phi_i,\phi_j \rangle 
\end{bmatrix},\quad 
\mathbf{A} = \begin{bmatrix}
\langle \nabla\phi_i,\nabla\phi_j \rangle 
\end{bmatrix}
$$

### 4. 单元矩阵计算
四面体单元$K$的局部矩阵：
$$
\mathbf{M}_K = \frac{|K|}{20}
\begin{bmatrix}
2 & 1 & 1 & 1 \\
1 & 2 & 1 & 1 \\
1 & 1 & 2 & 1 \\
1 & 1 & 1 & 2
\end{bmatrix},\quad
\mathbf{A}_K = \frac{|K|}{6}\mathbf{G}\mathbf{G}^T
$$

梯度矩阵$\mathbf{G}$通过坐标变换求得：
$$
\mathbf{G} = \begin{bmatrix}
1 & 1 & 1 & 1 \\
\mathbf{x}_1 & \mathbf{x}_2 & \mathbf{x}_3 & \mathbf{x}_4
\end{bmatrix}^{-1}
\begin{bmatrix}
0 & 0 & 0 \\
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}
$$

## 三、实验结果
<img src="D:\Project\MATLAB\工业有限元\ZJU-Personal-Course\Assignment 5\image\图片1.png" alt="有限元计算图像" width="750" style="margin: 0 auto;"/>

<div style="text-align:center">图1 有限元计算图像</div>


<img src="D:\Project\MATLAB\工业有限元\ZJU-Personal-Course\Assignment 5\image\误差曲线.png" alt="误差分析曲线" width="700" style="margin: 0 auto;"/>
<div style="text-align:center">图2 有误差分析曲线</div>

如上图所⽰，当⽹格密度从n=2增加到n=10时，相对误差不断降低，符合结果预期。

<img src="D:\Project\MATLAB\工业有限元\ZJU-Personal-Course\Assignment 5\image\heat3d.gif" alt="0-1s内三位热传导方程可视化" width="650" style="margin: 0 auto;"/>
<div style="text-align:center">图3 0-1s内三位热传导方程可视化</div>

如上图所⽰，设置网格为10，0-1s内热传导方程的可视化如上图gif所示。