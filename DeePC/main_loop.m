%% DeePC 仿真
clear,clc,close all
seed = 123;

%% 参数设置

% ===== 仿真参数： =====
dt = 1e-4;  % 仿真时长
T = 20;     % 仿真时间

% ===== 系统参数： =====
A = [0,1;-2,-3];
B = [0;1];
C = [1,0];
D = [0]; %#ok<*NBRAK2>

x_ini = [0;0];

% ===== 控制参数： =====
Ts = 0.02;  % 控制周期
method = 0; % 0:MPC    1:ID-MPC    2:DeePC
% ----- MPC 参数： -----
Q = diag([1]);
R = diag([5,5]);
N = 30;     % 预测步数
% ----- DeePC参数： -----
lambda_y = 1e5;     % 松弛正则因子
lambda_g = 200;     % 优化正则因子
initial_len = 20;    % 初始条数

% ===== 采集参数： =====
sampling_len = 200;  % 样本数目
sampling_mode = 1;  % 0:random    1:Guassian    2:sine    3:square
% ----- mode:0 随机参数 ----- 
range = [-1,1];
% ----- mode:1 Gaussian -----
sigma = 1;
mu = 0;
% ----- mode:2 正弦参数 & mode:3 方波参数 -----
angle_frequency = 2*pi;
amplitude = 1;

%% 初始化
% ===== 仿真变量 =====
tspan = 0:dt:T; % 仿真时间戳
% ===== 系统变量 =====
[n,m] = size(B);    % [系统维数,控制维数]
[p,~] = size(C);    % 输出维数
% ===== 激励信号 =====
rng(seed)
sampling_tspan = 0:dt:(sampling_len-1)*Ts;
switch sampling_mode
    case 0
        signal = @(t)range(1) + (range(2)-range(1))*rand(1);
    case 1
        signal = @(t)mu + sigma*randn(1);
    case 2
        signal = @(t)amplitude*sin(angle_frequency*t);
    case 3
        signal = @(t) amplitude*(2*(mod(round(t/Ts),round(2*pi/(angle_frequency*Ts))) ...
                 < round(2*pi/(angle_frequency*Ts))/2)-1);
end
% ===== 系统动态 =====
dynamic = @(x,u,t) A*x + B*u;
%% 数据生成
% ===== 初始化 =====
u = zeros(m,length(sampling_tspan));
x = zeros(n,length(sampling_tspan));
    x(:,1) = 0;
y = zeros(p,length(sampling_tspan));

u_data = zeros(m,sampling_len);
y_data = zeros(p,sampling_len);

% ===== 采集循环 =====
i = 1;  % 用于记录数据集
for k=1:length(sampling_tspan)
    t = sampling_tspan(k);
    if mod(k-1,round(Ts/dt))==0     % 更稳健的周期判断方法：用整数索引判断，而非直接浮点数取余
        u(:,k) = signal(t)+randn*0.0;
        u_data(:,i) = u(:,k);
        y_data(:,i) = C*x(:,k) + D*u(:,k);
        i = i+1;
    else
        u(:,k) = u(:,k-1);
    end
    x(:,k+1) = update_rk4(dynamic,x(:,k),u(:,k),t,dt);
    y(:,k) = C*x(:,k) + D*u(:,k);
end

figure(1)
% plot(sampling_tspan,u);
plot(u_data,'*-');
title("激励信号")

figure(2)
% plot(sampling_tspan,y);
plot(y_data,'x-');
title("采集输出")

[Hu,r,~] = HankelMat(initial_len+N+n,u_data);

s = svd(Hu);
figure(3)
semilogy(s,'o-');
grid on;
xlabel('奇异值序号');
ylabel('奇异值');

fprintf("The rank of Hankel Matrix is [%d].\n",r);

%% 控制循环

% ===== 初始化 =====
import casadi.*
x = zeros(n,length(tspan));
    x(:,1) = x_ini;
y = zeros(p,length(tspan));
u = zeros(m,length(tspan));

% ===== 仿真开始 =====
for k = 1:length(tspan)
    if mod(k-1,round(Ts/dt)) == 0
        % 进入控制周期
        switch method
            case 0 % MPC
                
        end
    else
    end
end
