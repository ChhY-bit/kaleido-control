clear,clc%,close all

%% 参数设置
% ========== 仿真参数 ==========
dt = 1e-4;
T = 20;
ifDelay = false;

% ========== 系统参数 ==========
Ac = [0,1;-1,-2];
Bc = [0;1];

% ========== 控制参数 ==========
Ts = 0.05;

N = 50;
Q = diag([2,1]);
R = 1;

W = diag([1,1]);
epsilon = 0;

warm_start = true;

% ========== 约束参数 =========
x_min = 1*[-0.5;-0.2];
x_max = 1*[1;0.2];
u_min = -0.2;
u_max = 0.2;

% ========== 任务参数 ==========
x_ini = [0.5;0];

% ========== 求解器 ==========
settings = struct( ...
    'verbose', false, ...
    'eps_abs', 1e-7, ...
    'eps_rel', 1e-7, ...
    'max_iter', 20000);

%% 初始化
% ---------- 系统初始化 ----------
sys = c2d(ss(Ac,Bc,[],[]),Ts);
A = sys.A;
B = sys.B;

tspan = 0:dt:T;
[n,m] = size(B);
dynamic_fun = @(x,u,t) Ac*x + Bc*u;

% ---------- 约束初始化 ----------
X_lb = kron(ones(N,1),x_min);
X_ub = kron(ones(N,1),x_max);
U_lb = kron(ones(N,1),u_min);
U_ub = kron(ones(N,1),u_max);


designer = ctr.terminal_designer(A,B,Q,R);
designer.setParams(W,epsilon,[x_min,x_max],[u_min,u_max]);
designer.design('lqr');

% controller = ctr.mpc_condensed;
controller = ctr.mpc_sparse;
controller.setProblem(A,B,N,Q,R);
% controller.setTerminal(designer.P,designer.Kt,designer.Et,designer.xt_min,designer.xt_max);
controller.init;
controller.setWarmStart(warm_start);
controller.setup('osqp',settings);

%% 控制循环
% ---------- 数据预分配 ----------
x = zeros(n,length(tspan));
x(:,1) = x_ini;
u = zeros(m,length(tspan));
solving_time = zeros(1,round(T/Ts));
run_time = zeros(1,round(T/Ts));
cost_value = zeros(1,round(T/Ts));
k = 0;

% ---------- 仿真循环 ----------
for i = 1:length(tspan)-1
    t = tspan(i);
    if mod((i-1),round(Ts/dt)) == 0 % ---------- 控制时刻 ----------
        tic;
        k = k+1;
        controller.update(x(:,i),X_lb,X_ub,U_lb,U_ub);
        controller.solve;
        u_temp = controller.U(:,1);
        solving_time(k) = controller.solving_time;
        cost_value(k) = controller.J;
        run_time(k) = toc;
        if ~ifDelay
            u(:,i) = u_temp;
        else
            % 求解结果延迟作用
            if controller.solving_time > Ts %#ok<*UNRCH>
                error("Solving time too long!");
            end
            % delayed_time = i + ceil(controller.solving_time/dt);    % 求解时间延迟
            delayed_time = i + ceil(run_time(k)/dt);    % 运行时间延迟
        end
        
    else                            % ---------- 保持时刻 ----------
        if ifDelay && i == delayed_time
            u(:,i) = u_temp;    % 执行求解结果
        else
            u(:,i) = u(:,i-1);  % 保持
        end
    end
    x(:,i+1) = ctr.update_rk4(dynamic_fun,x(:,i),u(:,i),t,dt);
end

%% 结果绘制
% 状态
figure(1);
subplot(2,1,1)
plot(tspan,x(1,:))
xlabel('t / s')
ylabel('x_1')

subplot(2,1,2)
plot(tspan,x(2,:))
xlabel('t / s')
ylabel('x_2')

% 控制
figure(2);

plot(tspan,u)
xlabel('t / s')
ylabel('u')

% 求解
figure(3)
bar(1:k,run_time)
hold on
bar(1:k,solving_time)
fprintf("[solving] warm: %d, avg: %.3f ms, max: %.3f ms.\n",warm_start,1000*mean(solving_time),1000*max(solving_time))
fprintf("[running] warm: %d, avg: %.3f ms, max: %.3f ms.\n",warm_start,1000*mean(run_time),1000*max(run_time))

% 代价
figure(4)
bar(1:k,cost_value)
