%% Vehicle Suspension
% 在此项目中所有涉及的角度，例如theta phi等，角度值均较小，故sin(theta)近似为theta，cos(theta)近似为1
t % 小车运动的当前时刻
v % 小车近似的平动速度
g % 重力加速度
m % 小车总质量（除去Beam）
I % 小车绕质心的转动惯量
Lb % 小车质心到后轮的距离
Lf % 小车质心到前轮的距离
yc % 小车质心在世界坐标系下的y坐标，yc原点的定义为小车初始受力平衡位置
yc_dot % 小车质心在世界坐标系y方向的速度
yc_ddot % 小车质心在世界坐标系y方向的加速度
theta % 小车车身相对水平方向逆时针旋转的角度
theta_dot % theta的角速度
theta_ddot % theta的角速度
theta0 %初始状态下小车质心到Beam转轴的连线与水平方向的夹角
d %小车质心到Beam转轴的长度
kb % 小车后轮阻尼器k值
cb % 小车后轮阻尼器c值
kf % 小车前轮阻尼器k值
cf % 小车前轮阻尼器c值
h_b=20*cos(v/100*(t-L/v)) % 小车后轮所在地面的高度
h_f=20*cos(v/100*t) % 小车前轮所在地面的高度  考虑：是否用 线积分/轮子转速 计算前后轮相位差 更精确
delta_yb = yc-Lb*theta-h_b % 小车后轮弹簧形变量
delta_yf = yc+Lf*theta-h_f % 小车前轮弹簧形变量
delta_yb_dot = yc_dot-Lb*theta_dot-diff(h_b,t) % 小车后轮弹簧阻尼器形变的速度
delta_yf_dot = yc_dot+Lf*theta_dot-diff(h_f,t) % 小车前轮弹簧阻尼器形变的速度
% 以小车车身作为研究对象，yc的零点
T=1/2*m*yc_dot^2+1/2*I*theta_dot^2 % 小车动能（平动能+转动能）
V=m*g*yc % 小车势能（重力势能）
L=T-V % 小车拉格朗日函数
% 在广义坐标yc下的拉格朗日方程
m*yc_ddot - mg = -kb*delta_yb-kf*delta_yf-cb*delta_yb_dot-cf*delta_yf_dot
% 在广义坐标theta下的拉格朗日方程
I*theta_ddot=kb*delta_yb*Lb+cb*delta_yb_dot*Lb-kf*delta_yf*Lf-cf*delta_yf_dot*Lf

%% Rigid  Beam

k % Beam阻尼器的k值
c % Beam阻尼器的c值
L1 % Beam旋转轴距离第一个阻尼器放置点的距离
L2 % Beam旋转轴距离第二个阻尼器放置点的距离
L3 % Beam旋转轴距离第三个阻尼器放置点的距离
L % Beam阻尼器弹力与阻力的力臂 等于L1/L2/L3
phi % Beam相对水平方向逆时针旋转的角度
phi_dot % phi的角速度
phi_ddot % phi的角加速度
phi-theta % Beam相对小车车身逆时针旋转的角度
lb % Beam旋转轴距离Beam末端(back)的长度
lf % Beam旋转轴距离Beam前段(front)的长度
mb % Beam的总质量
Ib % Beam绕旋转轴旋转的转动惯量
lbc=1/2*(lf-lb) % Beam旋转轴到质心的长度
yb=yc+d*(theta0-theta) % Beam旋转轴在世界坐标系下的y坐标
yb_dot=yc_dot-d*theta_dot % Beam旋转轴在世界坐标系y方向的速度
ybc=yb+lbc*phi % Beam质心在世界坐标系下的y坐标
ybc_dot=yc_dot-d*theta_dot+lbc*phi_dot % Beam质心在世界坐标系y方向的速度

Tb=1/2*mb*ybc_dot^2 + 1/2*Ib*phi_dot^2 % Beam的动能（平动能与转动能）
Vb=mb*g*ybc % Beam的势能（重力势能）
Lb=Tb-Vb % Beam的拉格朗日函数

% L对phi_dot求偏导，后对t求导
mb*lbc*phi_dot
mb*lbc*phi_ddot+Ib*phi_ddot
% L对phi求偏导
-mb*g*lbc
% 在phi广义坐标下的广义力矩
tau_b=k*L^2*(phi-theta)+c*L*(phi_dot-theta_dot)
% 在广义坐标phi下的拉格朗日方程
mb*lbc*phi_ddot+Ib*phi_ddot+mb*g*lbc = tau_b


%% Dynamic Beam
% Beam旋转轴后方视为rigid body，旋转轴前方视为dynamic body，旋转轴在积分时视作x原点
PHI = (3*x^2*lf-x^3)/(2*lf^3) % 形状函数 PHI(x)
PHI_prime % 形状函数对x求一阶导
PHI_pprime % 形状函数对x求二阶导
q % q(t)
q_dot % q(t)对t求一阶导数
q_ddot % q(t)对t求二阶导数
y = PHI * q % y(x,t)=PHI(x)*q(t) 表示悬臂梁某一点在某一时刻下相对刚性梁的变形幅度
lambda = mb/(lb+lf) % Beam的线密度
E % Beam的弹性模量
Iz % Beam的惯性二次矩

Tdb=1/2*int(lambda*(yb_dot+x*phi_dot)^2,x,-lb,0)+1/2*int(lambda*(yb_dot+x*phi_dot-PHI*q_dot)^2,x,-lb,0) % Dynamic Beam的动能，包含平动能和转动能
Vdb=(lambda*lb)*g*(yb-lb*phi)+int(lambda*g*(yb+x*phi-y),x,0,lf)+1/2*int(E*q^2*diff(PHI,x,2)*Iz,x,0,lf) % Dynamic Beam的势能，包含重力势能，以及Strain Energy
Ldb=Tdb-Vdb % Dynamic Beam的拉格朗日函数

% Stain Energy 推导过程

% L对q_dot求偏导，后对t求导
int(lambda*(-PHI)*(phi_ddot*x+yc_ddot-d*theta_ddot-PHI*q_ddot),x,0,lf)
% L对q求偏导
-int(lambda*g*(-PHI),x,0,lf)-int(E*q*(PHI_pprime)^2*Iz,x,0,lf)
% 在广义坐标q下的广义力
f_db = 0 % 旋转轴之前的Beam考虑为以旋转轴为端点的悬臂梁，形变由重力造成，故广义力为0
% 在广义坐标q下的拉格朗日方程
int(lambda*(-PHI)*(phi_ddot*x+yc_ddot-d*theta_ddot-PHI*q_ddot),x,0,lf)+int(lambda*g*(-PHI),x,0,lf)+int(E*q*(PHI_pprime)^2*Iz,x,0,lf) = 0
% 手动求导后在广义坐标q下的拉格朗日方程 
-11/40*lambda*phi_ddot*lf^2-3/8*lambda*lf*(yc_ddot-d*theta_ddot)+23/140*q_ddot*lf-3/8*lambda*g*lf+E*q*Iz*3/(lf)^3 = 0


% L对phi_dot求偏导，后对t求导
% L对phi求偏导
% 在广义坐标theta下的广义力矩
tau_b
% 在广义坐标theta下的拉格朗日方程
% 手动求导后在广义坐标theta下的拉格朗日方程 
lambda*(1/2*phi_ddot*(lf^2-lb^2)+(y_ddot-d*theta_dot)*(lf+lb)-3/8*q_ddot*lf)+lb*lambda*g*(y+d*(theta0-theta))+1/2*lambda*g*lf^2 = tau_b

%% 
