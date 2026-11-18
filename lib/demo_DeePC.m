clear,clc,close all

%% 仿真设置
rng(123);

%% 参数设置

% 被控对象模型
A = [0.921, 0.000, 0.041, 0.000;...
                 0.000, 0.918, 0.000, 0.033;...
                 0.000, 0.000, 0.924, 0.000;...
                 0.000, 0.000, 0.000, 0.937];
B = [0.017, 0.001;...
     0.001, 0.023;...
     0.000, 0.061;...
     0.072, 0.000];

C = [1,0,0,0; 0,1,0,0];

state_ini = [0;0;0;0];          % 初始状态

% 控制器参数
N = 30;             % 预测步数
T_ini = 10;         % 初始长度
Q = diag([35,35]);    % 输出权重
% R = diag([1e-4,1e-4]);    % 输入权重
R = diag([10,10]);    % 输入权重
lambda_u = 1e6;
lambda_y = 1e4;
lambda_a = 0;

Y_lb = [-2;-2];
Y_ub = [2;2];
U_lb = [-2;-2];
U_ub = [2;2];

% 离线收集数据
S = 400;           % 数据样本数
amplitude = 3;


% 在线控制
T = 200;    % 步数

ur = [0.22694;1.40062];
% ur = [0;0];
yr = [0.65;0.5];

%% 离线收集数据

% 信号初始化
u_ext = (2*rand(2,S)-1)*amplitude;
x_ext = zeros(4,S);
    x_ext(:,1) = state_ini;
y_ext = zeros(2,S);
    y_ext(:,1) = C*state_ini;

for k=1:S-1
    x_ext(:,k+1) = A*x_ext(:,k) + B*u_ext(:,k);
    y_ext(:,k+1) = C*x_ext(:,k+1);
end

u_data = u_ext;
y_data = y_ext;

% 加噪
% u_data = u_data + 0.002*(2*rand(size(y_data))-1);
% y_data = y_data + 0.002*(2*rand(size(y_data))-1);

% 预览数据
figure
plot(1:size(u_data,2),u_data)
figure
plot(1:size(y_data,2),y_data)

%% 在线控制

% 控制器初始化
settings = struct( ...
    'verbose', false, ...
    'eps_abs', 1e-8, ...
    'eps_rel', 1e-8, ...
    'max_iter', 50000);

controller = ctr.deepc(N,T_ini);
controller.setParams(Q,R,lambda_u,lambda_y,lambda_a);
controller.init(u_data,y_data);
controller.samplify(84,true);  % 手动指定保留奇异值
controller.setup('osqp',settings);
controller.setLogger(false);

% 数据初始化
x = zeros(4,T);
    x(:,1) = state_ini;
u = zeros(2,T);
y = zeros(2,T);
    y(:,1) = C*state_ini;

%% 控制循环
for k = 1:T-1
    if k >= T_ini+1 % 此前初始样本数不足
        % 进入控制周期：更新
        u_ini = u(:,k-T_ini:k-1);
        y_ini = y(:,k-T_ini:k-1);
        controller.update(u_ini(:),y_ini(:),...
                          kron(ones(N,1),ur),kron(ones(N,1),yr),...
                          kron(ones(N,1),Y_lb),kron(ones(N,1),Y_ub),...
                          kron(ones(N,1),U_lb),kron(ones(N,1),U_ub));
        controller.solve;
        u(:,k) = controller.u_series(:,1);
    end
    % 仿真推进
    u_act = clip(u(:,k),U_lb,U_ub);
    x(:,k+1) = A*x(:,k) + B*u_act;
    y(:,k+1) = C*x(:,k+1);
end

figure
plot(1:T,y(1,:),'b'); 
hold on
plot(1:T,yr(1).*ones(2,T),'--b');
plot(1:T,y(2,:),'r');
plot(1:T,yr(2).*ones(2,T),'--r');

figure
plot(1:T,u)
