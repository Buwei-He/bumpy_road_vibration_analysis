syms k c;
global A

M = [
    1,0,0,0;
    0,4,0,0;
    0,0,2,0
    0,0,0,7
    ]

C = [
    3,4,5,6;
    4,5,6,7;
    5,6,7,8;
    6,7,8,c
    ]

K = [
    1,1,0,0;
    1,1,1,0;
    0,1,1,1;
    0,0,1,k]

F = [1;2;3;4]

A = OdeMatrix(4,M,C,K,F)

%常微分数值解区间 interval
interval = [0,7];
%广义坐标初始位置
v_ini = [0,0,0,0];
%广义坐标初始速度
v_ini_dot = [0,0,0,0];

[t,x]=ode45(double(@Odefunc),interval,[v_ini;v_ini_dot])
%这里没做出来 似乎不让带参数计算

function A = OdeMatrix(num, M, C, K, F)
    % 确定原始方程变量数（3/4）；
    varN = num;
    order = 2;
    % 确定v和dv的变量数；此处阶数order固定为2；
    % 原始方程最高有二阶导项，所以v包含一阶导和原变量，dv有二阶导和一阶导
    dv = zeros(order*num,1)
    % 对M、C、K同时左乘inv(M)，得到归一化的Me、Ce、Ke
    % M需要是可逆矩阵…
    Me = eye(num);
    Ce = M\C
    Ke = M\K
    % 建立变换方程（我也不知道怎么叫）A
    % 表达式： (dv/dt) = Av + F(t)
    % 预期解：v = exp(A*t)-inv(A)*F
    % A矩阵应该是可逆的…
    A = zeros(order*varN,order*varN);
    for ie = 1:varN
        for je = 1:varN
            A(order*ie,order*je) = -Ce(ie,je);
            A(order*ie-1,order*je) = Me(ie,je);
            A(order*ie,order*je-1) = Me(ie,je);
        end
    end
end

function output = Odefunc(t,x)
    global A
    x_col = [x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8)]
    % 微分方程
    dx = A*x_col + F
end