clc; clear all;
close all
syms s t x...
     
cf=500000*2;
cb=500000*2;
c=500000;
kf=5000000*2;
kb=5000000*2;
k=5000000;

%% Vehicle Suspension
v =1.6*pi; % 小车近似的平动速度 拿Adams中轮子每秒2圈，配上模型中轮子周长算出来的 
g = 9.80665; % 重力加速度 check
m = 1964.110288522; % 小车总质量（除去Beam）check
J = 3263.8860842253; % 小车绕质心的转动惯量 单位是kg/平方米   3263.8860842253 kg-meter**2 check
% adams中小车质心坐标：-9.3360264826E-03, 0.4799051531, 2.4788406839E-02 y↑z→
% adams中左后轮轮毂中心坐标：-1.1921446343, 0.4090892672, 0.8906648264
% adams中左前轮轮毂中心坐标：-1.1948454679, 0.4073302891, -1.0928608775
Lb = 0.8906648264; % 小车质心到后轮的水平距离 check
Lf = 1.0928608775; % 小车质心到前轮的水平距离 check
theta0 =atan((0.879-0.4799051531)/0.5313848963)-0.04 ; %初始状态下小车质心到Beam转轴的连线与水平方向的夹角 弧度制 check
% adams中小车质心坐标：-9.3360264826E-03, 0.4799051531, 2.4788406839E-02 y↑z→
% adams中Beam转轴中心坐标：-0.18, 0.879, 0.5313848963
D = sqrt((0.879-0.4799051531)^2+(0.5313848963)^2); %小车质心到Beam转轴的直线长度 check
Dx = 0.5313848963; %小车质心到Beam转轴的水平长度 check
hb=0.02*(cos(v/0.1*(t-(Lb+Lf)/v))-1); % 小车后轮所在地面的高度
hf=0.02*(cos(v/0.1*t)-1); % 小车前轮所在地面的高度  考虑：是否用 线积分/轮子转速 计算前后轮相位差 更精确
hb_dot=diff(hb,t);
hf_dot=diff(hf,t);
% delta_yb = yc-Lb*theta-hb; % 小车后轮弹簧形变量
% delta_yf = yc+Lf*theta-hf; % 小车前轮弹簧形变量
% delta_yb_dot = yc_dot-Lb*theta_dot-diff(hb,t); % 小车后轮弹簧阻尼器形变的速度
% delta_yf_dot = yc_dot+Lf*theta_dot-diff(hf,t); % 小车前轮弹簧阻尼器形变的速度
% 以小车车身作为研究对象，yc的零点
% Tc=1/2*m*yc_dot^2+1/2*J*theta_dot^2; % 小车动能（平动能+转动能）
% Vc=m*g*yc; % 小车势能（重力势能）
% Lc=Tc-Vc; % 小车拉格朗日函数
% 在广义坐标yc下的拉格朗日方程
% m*yc_ddot - mg = -kb*delta_yb-kf*delta_yf-cb*delta_yb_dot-cf*delta_yf_dot
% 在广义坐标theta下的拉格朗日方程
% J*theta_ddot=kb*delta_yb*Lb+cb*delta_yb_dot*Lb-kf*delta_yf*Lf-cf*delta_yf_dot*Lf

%% Rigid  Beam

L1 = 0.36; % Beam旋转轴距离第一个阻尼器放置点的距离 check
L2 = 0.48; % Beam旋转轴距离第二个阻尼器放置点的距离 check
L3 = 0.60; % Beam旋转轴距离第三个阻尼器放置点的距离 check
L=L1; % Beam阻尼器弹力与阻力的力臂 等于L1/L2/L3
lb = 0.74; % Beam旋转轴距离Beam末端(back)的长度 check
lf = 2.12; % Beam旋转轴距离Beam前段(front)的长度 check 
lambda = 15.602; % 前段Beam的线密度 kg/m check
mb = lambda*lf+2.5*lambda*lb; % Beam的总质量  31.204 kg  36.9806179046 kg check
Jb = int(2.5*lambda*x^2,x,-lb,0)+int(lambda*x^2,x,0,lf); % Beam绕旋转轴旋转的转动惯量 单位是kg/平方米 check
lbc=(lambda*lf*lf/2+2.5*lambda*lb*(-lb/2))/mb; % Beam旋转轴到质心的长度(旋转轴x坐标为0) check
% yb=yc+D*(theta0-theta); % Beam旋转轴在世界坐标系下的y坐标
% yb_dot=yc_dot-D*theta_dot; % Beam旋转轴在世界坐标系y方向的速度
% ybc=yb+lbc*phi; % Beam质心在世界坐标系下的y坐标
% ybc_dot=yc_dot-D*theta_dot+lbc*phi_dot; % Beam质心在世界坐标系y方向的速度

% Trb=1/2*mb*ybc_dot^2 + 1/2*Jb*phi_dot^2; % Beam的动能（平动能与转动能）
% Vrb=mb*g*ybc; % Beam的势能（重力势能）
% Lrb=Trb-Vrb; % Beam的拉格朗日函数

%% Dynamic Beam
% Beam旋转轴后方视为rigid body，旋转轴前方视为dynamic body，旋转轴在积分时视作x原点

% PHI = (3*x^2*lf-x^3)/(2*lf^3); % 形状函数 PHI(x)
% PHI_prime = 3/lf^2*x - 3/(2*lf^3)*x^2; % 形状函数对x求一阶导
% PHI_pprime = diff(PHI_prime,x);  % 形状函数对x求二阶导
% syms q % q(t)
% syms q_dot % q(t)对t求一阶导数
% syms q_ddot % q(t)对t求二阶导数
% y = PHI * q; % y(x,t)=PHI(x)*q(t) 表示悬臂梁某一点在某一时刻下相对刚性梁的变形幅度
E = 2.07*10^11; % Beam的弹性模量 adams.steel check
Iz = 0.1*0.02^3/12; % Beam的惯性二次矩 单位m的四次方 check

% Tdb=1/2*int(lambda*(yb_dot+x*phi_dot)^2,x,-lb,0)+1/2*int(lambda*(yb_dot+x*phi_dot-PHI*q_dot)^2,x,-lb,0); % Dynamic Beam的动能，包含平动能和转动能
% Vdb=(lambda*lb)*g*(yb-lb*phi)+int(lambda*g*(yb+x*phi-y),x,0,lf)+1/2*int(E*q^2*diff(PHI,x,2)*Iz,x,0,lf); % Dynamic Beam的势能，包含重力势能，以及Strain Energy
% Ldb=Tdb-Vdb; % Dynamic Beam的拉格朗日函数

% Stain Energy 推导过程

% Ldb对q_dot求偏导，后对t求导
% int(lambda*(-PHI)*(phi_ddot*x+yc_ddot-D*theta_ddot-PHI*q_ddot),x,0,lf);
% Ldb对q求偏导
% -int(lambda*g*(-PHI),x,0,lf)-int(E*q*(PHI_pprime)^2*Iz,x,0,lf);
% 在广义坐标q下的广义力
% f_db = 0; % 旋转轴之前的Beam考虑为以旋转轴为端点的悬臂梁，形变由重力造成，故广义力为0
% 在广义坐标q下的拉格朗日方程
% int(lambda*(-PHI)*(phi_ddot*x+yc_ddot-d*theta_ddot-PHI*q_ddot),x,0,lf)+int(lambda*g*(-PHI),x,0,lf)+int(E*q*(PHI_pprime)^2*Iz,x,0,lf) = 0
% 手动求导后在广义坐标q下的拉格朗日方程 
% -11/40*lambda*phi_ddot*lf^2-3/8*lambda*lf*(yc_ddot-d*theta_ddot)+23/140*q_ddot*lf-3/8*lambda*g*lf+E*q*Iz*3/(lf)^3 = 0


% Ldb对phi_dot求偏导，后对t求导
% Ldb对phi求偏导
% 在广义坐标theta下的广义力矩
% tau_b;
% 在广义坐标theta下的拉格朗日方程
% 手动求导后在广义坐标theta下的拉格朗日方程 
% lambda*(1/2*phi_ddot*(lf^2-lb^2)+(y_ddot-d*theta_dot)*(lf+lb)-3/8*q_ddot*lf)+lb*lambda*g*(y+d*(theta0-theta))+1/2*lambda*g*lf^2 = tau_b
%% Rigid Beam 方程求解

M_R = [
m, 0, 0;
0, J, 0;
mb*lbc, -mb*lbc*D, (Jb+mb*lbc^2)
];

C_R = [
(cf+cb), (Lf*cf-Lb*cb), 0;
(Lf*cf-Lb*cb), (Lf^2*cf+Lb^2*cb), 0;
0, -(c*L^2), (c*L^2)
];

K_R = [
(kf+kb), (Lf*kf-Lb*kb), 0;
(Lf*kf-Lb*kb), (Lf^2*kf+Lb^2*kb), 0;
0, -(k*L^2), (k*L^2)
];

F_R = [(cb*hb_dot + cf*hf_dot + kb*hb +kf*hf + m*g);
    (-Lb*cb*hb_dot + Lf*cf*hf_dot - Lb*kb*hb + Lf*kf*hf);
    (-mb*g*lbc)];

%% Dynamic Beam 计算部分

%M矩阵(Dynamic)
%修订1: 更改了M32和M42的符号，逆时针为正…
M31= lambda*(1/2*lf^2-5/4*lb^2);
M32 = -D*M31;
M33 = lambda*(1/3*lf^3+5/6*lb^3); 
M34=-11/40*lf^2*lambda;
M41=-3/8*lf*lambda;
M42 = -D*M41;
M43=-11/40*lf^2*lambda;
M44= 33/140*lf*lambda;

M_D=[
m  0  0  0;
0   J   0  0;
M31 M32 M33 M34;
M41 M42 M43 M44
];

%C矩阵

c11=cf+cb;         
c12=cf*Lf-cb*Lb;                     
c13=0 ;
c21=cf*Lf-cb*Lb; 
c22=cf*Lf^2 + cb*Lb^2;   
c23=0;
c31=0;               
c32=-c*L^2;                          
c33=c*L^2;

C_D=[
c11 c12 c13 0;
c21 c22 c23 0;
0     c32 c33 0;
0  0  0  0
];
 
%K矩阵

% k11=kf+kb;         
% k12=kf*Lf-kb*Lb-k*L;                
% k13=k*L;
% k21=kf*Lf-kb*Lb; 
% k22=kf*Lf^2+kb*Lb^2+k*L*(L+Dx);   
% k23=-k*L*(L+Dx);
% k31=0;               
% k32=-k*L^2;                      
% k33=k*L^2;                               
% k44=E*Iz*3/lf^3;

k11=kf+kb;         
k12=kf*Lf-kb*Lb;                
k13=0;
k21=kf*Lf-kb*Lb; 
k22=kf*Lf^2+kb*Lb^2;   
k23=0;
k31=0;               
k32=-k*L^2;                      
k33=k*L^2;                               
k44=E*Iz*3/lf^3;

K_D=[
   k11 k12 k13 0  ;
   k21 k22 k23 0  ;
   k31 k32 k33 0  ;
   0   0   0   k44
];

% 预先变换
hfs = laplace(hf,t,s);
hbs = laplace(hb,t,s);
hf_dot_s = laplace(hf_dot,t,s);
hb_dot_s = laplace(hb_dot,t,s);

F_Ds=[
    kf*hfs + cf*hf_dot_s + kb*hbs + cb*hb_dot_s + m*g/s;
    cf*Lf*hf_dot_s-cb*Lb*hb_dot_s+kf*Lf*hfs-kb*Lb*hbs;
    lambda*g*(5/4*lb^2-1/2*lf^2)/s;
    3/8*lambda*g*lf/s
    ];


%X的拉普拉斯变换后化简后的式子。

Lx=(M_D*s^2+C_D*s+K_D)\F_Ds;

ILx=ilaplace(Lx,s,t);
ILx=vpa(ILx,200);
global yc theta psi q 
yc = ILx(1,1);
theta = ILx(2,1);
psi = ILx(3,1);
q = ILx(4,1);
Delta = yc + D*(theta0  - theta)+ lf*psi- q - 0.3981506469;

ANS = simplify(int(Delta^2*v,t,0,33.91593/v))

integrate = 0;
% subs(Delta,t,1)
for i = 0 : 0.01 : 7
    integrate = integrate + ((subs(Delta,t,i)+subs(Delta,t,i+0.01))/2)^2 * v * 0.01;
end
integrate

figure 
fsurf(q, [0, 10])
view(0, 0)
figure
fsurf(Delta, [0, 10])
view(0, 0)
% 
% function f = delta(t)
% global yc theta psi q 
%     D = sqrt((0.879-0.4799051531)^2+(0.5313848963)^2);
%     lf = 2.12;
%     theta0 = atan((0.879-0.4799051531)/0.5313848963)-0.04;
%     f = yc + D*(theta0  - theta)+ lf*psi- q - 0.3981506469;
% end


