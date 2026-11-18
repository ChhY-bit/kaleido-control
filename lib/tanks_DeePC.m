clear,clc,close all

%% 仿真设置
dt = 1e-3;      % 仿真步长
rng(123);

%% 参数设置

% 被控对象模型
yaml_filepath = "quadruple_tank.yaml";  % 参数路径
state_ini = [0.2;0.3;0.1;0.4];          % 初始状态（水位[m]）

% 控制器参数
Ts = 0.02;          % 控制周期
N = 20;             % 预测步数
T_ini = 10;         % 初始长度
Q = diag([35,35]);    % 输出权重
R = diag([1e-4,1e-4]);    % 输入权重
lambda_u = 1e6;
lambda_y = 1e4;
lambda_a = 1e2;

% 离线收集数据
S = 400;           % 数据样本数
exciting_mode = 0;  % 0:rand    1:freq_sweep    2:test
    amplitude = 1;  % 扫频最大值（信号范围[0,amplitude]）
    omega_min = 0;  % 扫频最小值
    omega_max = 10; % 扫频最大值
    phase = pi/2;   % 2个输入之间的相位差

% 在线控制
T = 2;

ur = [0;0];
% yr = [0.5;0.1;0.2;0.3];
yr = [0.65;0.77];

%% 离线收集数据

% 被控对象初始化
params = rpm.utils.load_quadruple_tank_params(yaml_filepath);
tanks = rpm.quadruple_tank(params,state_ini);
C = [eye(2),zeros(2)];  % 输出矩阵
% C = eye(4);

% 信号初始化
T_ext = S*Ts;   % 总激励时长
tspan_ext = 0:dt:T_ext;
u_ext = zeros(2,length(tspan_ext));
x_ext = zeros(4,length(tspan_ext));
    x_ext(:,1) = state_ini;
y_ext = zeros(size(C,1),length(tspan_ext));
    y_ext(:,1) = C*state_ini;

% 确定激励信号模式
omega_spd = (omega_max-omega_min)/T_ext;
switch exciting_mode
    case 0
        signal = @(t) amplitude*rand();
    case 1
        signal = @(t) amplitude.*(0.5*...
                      cos((omega_min+0.5*omega_spd.*t).*t)+0.5);
end


for i=1:length(tspan_ext)-1
    t = tspan_ext(i);       % 当前时刻
    if mod((i-1),round(Ts/dt)) == 0
        % 进入控制周期：更新
        u_ext(:,i) = [signal(t);...
                      signal(t+phase)]; % 输入信号样本
    else
        % 非控制周期：保持
        u_ext(:,i) = u_ext(:,i-1);
    end

    % 仿真推进
    tanks.send(u_ext(1,i),u_ext(2,i));
    tanks.step(dt,t);
    x_ext(:,i+1) = [tanks.x1;tanks.x2;tanks.x3;tanks.x4];
    y_ext(:,i+1) = C*x_ext(:,i+1);
end

% 预览数据
figure
plot(tspan_ext,u_ext)

figure
plot(tspan_ext,y_ext)

% 抽取数据
u_data = u_ext(:,1:round(Ts/dt):length(tspan_ext)-1);
y_data = y_ext(:,1:round(Ts/dt):length(tspan_ext)-1);

% 加噪
u_data = u_data + 0.02*rand(size(u_data,1),size(u_data,2));
y_data = y_data + 0.02*rand(size(y_data,1),size(y_data,2));

figure
plot(1:size(u_data,2),u_data)

figure
plot(1:size(y_data,2),y_data)

% 渲染画面
% rander_tanks(x_ext, u_ext, params, dt, 10, 'sampling.gif');

%% 在线控制

% 对象初始化
tanks.reset(state_ini);
Y_lb = C*zeros(4,1);
Y_ub = C*[params.TankHeight_1;params.TankHeight_2;...
        params.TankHeight_3;params.TankHeight_4];
U_lb = zeros(2,1);
U_ub = [params.PumpMaxFlow_1;params.PumpMaxFlow_2];

% 控制器初始化
settings = struct( ...
    'verbose', false, ...
    'eps_abs', 1e-8, ...
    'eps_rel', 1e-8, ...
    'max_iter', 50000);

controller = ctr.deepc(N,T_ini);
controller.setParams(Q,R,lambda_u,lambda_y,lambda_a);
controller.init(u_data,y_data);
controller.samplify(62,true);  % 手动指定保留奇异值
controller.setup('osqp',settings);
controller.setLogger(false);

% 数据初始化
tspan = 0:dt:T;
x = zeros(4,length(T));
    x(:,1) = state_ini;
u = zeros(2,length(tspan));
y = zeros(size(C,1),length(tspan));
    y(:,1) = C*state_ini;

u_rec = zeros(2,round(T/Ts));
y_rec = zeros(size(C,1),round(T/Ts));

%% 控制循环
k = 0;  % 数据样本计数器
for i = 1:length(tspan)-1
    t = tspan(i);
    if mod((i-1),round(Ts/dt)) == 0
        k = k+1;
        if k >= T_ini+1 % 此前初始样本数不足
            % 进入控制周期：更新
            u_ini = u_rec(:,k-T_ini:k-1);
            y_ini = y_rec(:,k-T_ini:k-1);
            controller.update(u_ini(:),y_ini(:),...
                              kron(ones(N,1),ur),kron(ones(N,1),yr),...
                              kron(ones(N,1),Y_lb),kron(ones(N,1),Y_ub),...
                              kron(ones(N,1),U_lb),kron(ones(N,1),U_ub));
            controller.solve;
            u(:,i) = controller.u_series(:,1);
            
            % 离散信息
            u_rec(:,k) = u(:,i);
            y_rec(:,k) = y(:,i);
        end
    else
        % 非控制周期：保持
        u(:,i) = u(:,i-1);
    end
    % 仿真推进
    u_act = clip(u(:,i),U_lb,U_ub);
    tanks.send(u_act(1),u_act(2));
    tanks.step(dt,t);
    x(:,i+1) = [tanks.x1;tanks.x2;tanks.x3;tanks.x4];
    y(:,i+1) = C*x(:,i+1);
end

figure
subplot(2,1,1)
plot(tspan,y(1,:))
hold on
plot(tspan,yr(1).*ones(2,length(tspan)),'--r');
subplot(2,1,2)
plot(tspan,y(2,:))
hold on
plot(tspan,yr(2).*ones(2,length(tspan)),'--r');

figure
plot(tspan,u)
