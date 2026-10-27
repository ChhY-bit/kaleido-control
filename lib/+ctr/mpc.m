classdef mpc < handle
%MPC 标准线性模型预测控制

    properties(SetAccess=private)
        % ---------- 基本属性 ----------
        A   % [nxn] 系统矩阵
        B   % [nxm] 输入矩阵
        C   % [pxn] 输出矩阵
        D   % [pxm] 前馈矩阵

        Q   % [pxp] 跟踪权重矩阵
        R   % [mxm] 控制权重矩阵
        S   % [mxm)] 增量权重矩阵

        N   % 预测步数

        % ---------- 终端要素 ----------
        P       % [pxp] 终端权重矩阵
        K       % [mxn] 终端反馈增益
        yt_max  % 终端约束上界
        yt_min  % 终端约束下界
        
        % ---------- 约束界限 ----------
        % （可能在线更新）
        Y_lb   % [pxN] 输出约束下界
        Y_ub   % [pxN] 输出约束上界
        U_lb   % [mxN] 控制约束下界
        U_ub   % [mxN] 控制约束上界
        Delta_lb  % [mx(N-1)] 增量约束下界
        Delta_ub  % [mx(N-1)] 增量约束上界

        % ---------- 设置标志 ----------
        ifWarmStart
        ifLogger
        ifTiming
    end

    properties(Access = private)
        % ---------- 在线存储 ----------
        % 在mpc.setup()中被在线计算
        G_lb
        G_ub
        f
        c
        J       % 代价函数值
        x_ini   % 预测起点
        r       % [px(>=N)] 参考信号

        % ---------- 离线存储 ----------
        % 在 mpc.init()中被离线计算
        n,m,p   % 维数
        QQ, RR, SS  % 增广权重
        T_x
        T_u
        S_Delta
        Gamma_u
        C_t % [pxn]
        H   % [(mN)x(mN)]
        F   % [((2m+p)N-m+p)x(mN)]
    end


    properties(Access = public)
        U               % 求解结果
        solving_time    % 求解时间
        U_f             % 可行解，移位补终端
    end

    methods(Access=public)
        function obj = mpc(A,B,C,D,N,Q,R,S)
            obj.A = A; obj.B = B; obj.C = C; obj.D = D;
            obj.N = N; obj.Q = Q; obj.R = R; obj.S = S;
            obj.init();
        end

        function setTerminal(obj,P,K,yt_min,yt_max)
            obj.P = P;
            obj.K = K;
            obj.yt_max = yt_max;
            obj.yt_min = yt_min;
        end

        function setup(obj,x_ini,Y_b,U_b,Delta_b)
            if ~isempty(Y_b)
                obj.Y_lb = Y_b(:,1);
                obj.Y_ub = Y_b(:,2);
            else
                obj.Y_lb = [];
                obj.Y_ub = [];
            end
            if ~isempty(U_b)
                obj.U_lb = U_b(:,1);
                obj.U_ub = U_b(:,2);
            else
                obj.U_lb = [];
                obj.U_ub = [];
            end
            if ~isempty(Delta_b)
                obj.Delta_lb = Delta_b(:,1);
                obj.Delta_ub = Delta_b(:,2);
            else
                obj.Delta_lb = [];
                obj.Delta_ub = [];
            end
            obj.x_ini = x_ini;
        end

        function U = solve(obj,solver,options)
            Phi = obj.T_x*obj.x_ini - obj.r(:,1:obj.N);
            Psi = obj.C_t*obj.A^obj.N*obj.x_ini-obj.r(:,obj.N+1);
            obj.f = (Phi'*obj.P*obj.C_t*obj.Gamma_u + Phi'*obj.QQ*obj.T_u)';
            obj.c = Psi'*obj.P*Psi + Phi'*obj.QQ*Phi;

            switch solver
                case 'osqp'
                    sol = osqp;
                    sol.setup(obj.H,obj.f,obj.F,obj.G_lb,obj.G_ub);

                case 'default'
                otherwise
                    error("No solver named '%s'!",solver);
            end
        end

        function setLogger(obj,bool)
        end

        function setWarmStart(obj,bool)
        end

        function setTiming(obj,bool)
            % 记录的是求解时间，而非整个控制决策用时
            % 后者还包括状态观测、变量传入等耗时
        end
    end

    methods(Access=private)
        function init(obj)
        % INIT 进行初始化离线计算
            % ----- 系统维度 -----
            [obj.n,obj.m] = size(obj.B);
            [obj.p,~] = size(obj.C);
            
            % ----- 考虑终端情形 -----
            if isempty(obj.K)
                obj.C_t = obj.C;
                obj.P = zeros(obj.p);
                % obj.yt_max = inf.*ones(obj.p,1);  %（考虑在求解器中直接处理空情况）
                % obj.yt_min = -inf.*ones(obj.p,1);
            else
                obj.C_t = obj.C+obj.D*obj.K;
            end

            % ----- 权重矩阵 -----
            obj.QQ = kron(eye(obj.N),obj.Q);
            obj.RR = kron(eye(obj.N),obj.R);
            obj.SS = kron(eye(obj.N-1),obj.S);

            % ----- 计算 T_u -----
            obj.T_u = zeros(obj.p*obj.N,obj.m*obj.N);
            for i = 1:obj.N
                for j = 1:i
                    if j==i
                        obj.T_u((i-1)*obj.p+1:i*obj.p,(j-1)*obj.m+1:j*obj.m) = obj.D;
                        continue
                    else
                        obj.T_u((i-1)*obj.p+1:i*obj.p,(j-1)*obj.m+1:j*obj.m) = obj.C*obj.A^(i-1-j)*obj.B;
                    end
                end
            end

            % ----- 计算 T_x -----
            obj.T_x = zeros(obj.p*obj.N,obj.n);
            for i = 1:obj.N
                obj.T_x((i-1)*obj.p+1:i*obj.p,:) = obj.C*obj.A^(i-1);
            end

            % ----- 计算 S_Delta -----
            obj.S_Delta = zeros(obj.m*(obj.N-1),obj.m*obj.N);
            for i = 1:obj.N-1
                obj.S_Delta((i-1)*obj.m+1:i*obj.m,(i-1)*obj.m+1:i*obj.m) = -eye(obj.m);
                obj.S_Delta((i-1)*obj.m+1:i*obj.m,i*obj.m+1:(i+1)*obj.m) = eye(obj.m);
            end

            % ----- 计算 Gamma_u -----
            obj.Gamma_u = zeros(obj.n,obj.m*obj.N);
            for j = 1:obj.N
                obj.Gamma_u(:,(j-1)*obj.m+1:j*obj.m) = obj.A^(obj.N-j)*obj.B;
            end

            % % ----- 构建 H 和 F -----
            obj.H = obj.Gamma_u'*obj.C_t'*obj.P*obj.C_t*obj.Gamma_u + ...
                    obj.T_u'*obj.QQ*obj.T_u + obj.RR + ...
                    obj.S_Delta'*obj.SS*obj.S_Delta;
            obj.F = [eye(obj.m*obj.N);...
                     obj.S_Delta;...
                     obj.T_u;
                     obj.C_t*obj.Gamma_u];
        end

    end
end