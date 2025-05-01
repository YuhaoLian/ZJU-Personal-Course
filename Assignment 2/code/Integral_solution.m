clc;
clear;

%% <ph0',ph0'>
disp('计算 <ph0”,ph0”> 的不定积分和定积分：');
syms x;
f = (4*x - 3 / 2)^2;
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph0,ph0>
disp('计算 <ph0,ph0> 的不定积分和定积分：');
syms x;
f = (2*(x-1)*(x-1/2))^2;
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph0',ph1/2'>
disp('计算 <ph0”,ph1/2”> 的不定积分和定积分：');
syms x;
f = (4*x - 3 / 2) * (4 - 8 * x);
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph0,ph1/2>
disp('计算 <ph0,ph1/2> 的不定积分和定积分：');
syms x;
f = (2*(x-1)*(x-1/2))*(4 * x * (1 - x));
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph0',ph1'>
disp('计算 <ph0”,ph1”> 的不定积分和定积分：');
syms x;
f = (4*x - 3 / 2) * (4 * x - 1);
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph0,ph1>
disp('计算 <ph0,ph1> 的不定积分和定积分：');
syms x;
f = (2*(x-1)*(x-1/2)) * (2 * x * (x - 1/2));
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph1/2',ph1/2'>
disp('计算 <ph1/2”,ph1/2”> 的不定积分和定积分：');
syms x;
f = (4 - 8 * x)^2;
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph1/2,ph1/2>
disp('计算 <ph1/2,ph1/2> 的不定积分和定积分：');
syms x;
f = (4 * x * (1 - x))^2;
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph1/2',ph1'>
disp('计算 <ph1/2”,ph1”> 的不定积分和定积分：');
syms x;
f = (4 - 8 * x) * (4 * x - 1);
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph1/2,ph1>
disp('计算 <ph1/2,ph1> 的不定积分和定积分：');
syms x;
f = (4 * x * (1 - x)) * (2 * x * (x - 1/2));
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph1',ph1'>
disp('计算 <ph1”,ph1”> 的不定积分和定积分：');
syms x;
f = (4 * x - 1)^2;
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);

%% <ph1,ph1>
disp('计算 <ph1,ph1> 的不定积分和定积分：');
syms x;
f = (2 * x * (x - 1/2))^2;
y = int(f,x);
disp(y);
result = int(f, x, 0, 1);
disp(result);