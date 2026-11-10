clear,clc,%close all

% save_dir = 'C:\Users\BIT_robot\Desktop\delay_test';
% save_name = 'no_delay.mat';
%% 参数设置
% ========== 仿真参数 ==========
dt = 1e-4;
T = 10;
ifDelay = 0;

% ========== 系统参数 ==========
Ac = [1,0;0,2];
Bc = [1;1];

% ========== 控制参数 ==========
Ts = 0.05;

N = 40;
Q = diag([2,1]);
R = 1;

W = diag([1,1]);
epsilon = 0;

warm_start = true;
logger_on = false;

% ========== 约束参数 =========
x_min = [-2;-1];
x_max = [2;1];
u_min = -1.5;
u_max = 1.5;

% ========== 任务参数 ==========
x_ini = [0.5;0];

% ========== 求解器 ==========
settings = struct( ...
    'verbose', false, ...
    'eps_abs', 1e-8, ...
    'eps_rel', 1e-8, ...
    'max_iter', 50000);

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
designer.design('lqr',[],false);

% controller = ctr.mpc_condensed;
controller = ctr.mpc_sparse;
controller.setProblem(A,B,N,Q,R);
controller.setTerminal(designer.P,designer.Kt,designer.Et,designer.xt_min,designer.xt_max);
controller.init;
controller.setWarmStart(warm_start);
controller.setLogger(logger_on);
controller.setup('osqp',settings);

%% 控制循环
% ---------- 数据预分配 ----------
x = zeros(n,length(tspan));
x(:,1) = x_ini;
u = zeros(m,length(tspan));
delay_flag = zeros(1,length(tspan));    % 记录延迟时刻
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
        % ----- 前向预测补偿延迟 -----
        % sys_p = d2d(sys,500*dt);
        % Ap = sys_p.A; Bp = sys_p.B;
        % controller.update(Ap*x(:,i)+Bp*u(:,i),X_lb,X_ub,U_lb,U_ub);
        % ----------------------------
        controller.solve;
        u_temp = controller.U(:,1);
        solving_time(k) = controller.solving_time;
        cost_value(k) = controller.J;
        run_time(k) = toc;
        if ifDelay
            % 模拟求解结果延迟作用
            % delay_steps = ceil(controller.solving_time/dt);
            delay_steps = 0*ceil(run_time(k)/dt)+500; %#ok<*UNRCH>
            if i+delay_steps <= length(tspan)
                u(:,i+delay_steps) = u_temp;
            end
            if i > 1 && ~delay_flag(i)
                u(:,i) = u(:,i-1);
            end
            delay_flag(i+delay_steps) = 1;
            
            % if delay_steps >= round(Ts/dt) %
            %     error("Solving time too long!");
            % end
            % delayed_time = i + delay_steps;
            % if delay_steps == 0
            %     u(:,i) = u_temp;       % 当步执行
            % elseif i > 1
            %     u(:,i) = u(:,i-1);     % 等待期间保持旧输入
            % else
            %     u(:,i) = zeros(m,1);   % 首次等待期间的初始输入
            % end
        else
            % delayed_time = i;
            u(:,i) = u_temp;
        end
        
    else                            % ---------- 保持时刻 ----------
        % if ifDelay && i == delayed_time
        %     u(:,i) = u_temp;    % 执行求解结果
        % else
        if ~delay_flag(i)
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
hold on

subplot(2,1,2)
plot(tspan,x(2,:))
xlabel('t / s')
ylabel('x_2')
hold on

% 控制
figure(2);

plot(tspan,u)
xlabel('t / s')
ylabel('u')
hold on

% 求解
figure(3)
clf
bar(1:k,[run_time(1:k).',solving_time(1:k).'],'grouped')
xlabel('Control step')
ylabel('Time / s')
legend('Running time','Solving time','Location','best')
grid on
fprintf("[solving] warm: %d, avg: %.3f ms, max: %.3f ms, freq_min = %.0f Hz.\n",...
        warm_start,1000*mean(solving_time),1000*max(solving_time),1/max(solving_time))
fprintf("[running] warm: %d, avg: %.3f ms, max: %.3f ms, freq_min = %.0f Hz.\n",...
        warm_start,1000*mean(run_time),1000*max(run_time),1/max(run_time))

% 代价
figure(4)
bar(1:k,cost_value)
hold on

%% 保存数据
% if ~isfolder(save_dir)
%     mkdir(save_dir);
% end
% save_path = fullfile(save_dir,save_name);
% save(save_path,'x','u','tspan');
% fprintf('[save] Data saved to: %s\n',save_path);
