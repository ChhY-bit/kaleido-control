classdef terminal_designer < handle
    %UNTITLED2 此处显示有关此类的摘要
    %   此处显示详细说明

    properties
        % ---------- 基本属性 ----------
        A   % [nxn] 系统矩阵
        B   % [nxm] 输入矩阵

        Q   % [nxn] 状态权重矩阵
        R   % [mxm] 控制权重矩阵

        % ---------- 外部约束 ----------
        x_ub
        x_lb
        u_ub
        u_lb

        % ---------- 终端要素 ----------
        P       % [nxn] 终端权重矩阵
        Kt      % [mxn] 终端反馈增益，约定 u_t = Kt*x；可为空
        xt_max  % [ntx1] Et*x 的终端约束上界
        xt_min  % [ntx1] Et*x 的终端约束下界
        Et      % [ntxn] 终端约束系数矩阵；为空时不启用终端约束

        % ---------- 中间变量 ----------
        M       % Lyapunov 解
        alpha   % 终端集大小
        At      % 终端闭环矩阵
        % ---------- 设计参数 ----------
        W       % [nxn] Lyapunov 正定矩阵
        epsilon % 标量，终端代价权重裕度
    end

    properties(Access=private)
        % ---------- 设置标志 ----------
        ifLogger = true
        ifPlot = true
    end

    methods
        function obj = terminal_designer(A,B,Q,R)
            obj.A = A;
            obj.B = B;
            obj.Q = Q;
            obj.R = R;
            obj.logInfo("Terinal set designer is created!");
        end

        function setParams(obj,W,epsilon,xb,ub)
            obj.W = W;
            obj.epsilon = epsilon;
            obj.x_lb = xb(:,1);
            obj.x_ub = xb(:,2);
            obj.u_lb = ub(:,1);
            obj.u_ub = ub(:,2);
            obj.logInfo("Terinal set designer is ready!");
        end

        function design(obj,Kt_type,params,ifplot)
            n = size(obj.A,1);
            % 确定增益
            switch lower(Kt_type)
                case 'lqr'
                    obj.logInfo("LQR is used for terminal control...");
                    obj.Kt = -dlqr(obj.A,obj.B,obj.Q,obj.R);
                case 'place'
                    obj.logInfo("Pole-placement is used for terminal control...");
                    obj.Kt = -place(obj.A,obj.B,params);
                otherwise
                    logError(sprintf("No method named '%s'!",Kt_type));     
            end
            obj.At = obj.A + obj.B*obj.Kt;

            % 终端代价权重
            obj.P = dlyap(obj.At',obj.Q+obj.Kt'*obj.R*obj.Kt+obj.epsilon*eye(n));

            % 确定形状
            obj.logInfo("Terminal rigion generated...");
            obj.M = dlyap(obj.At',obj.W);

            % 确定大小
            obj.alpha = obj.computeAlpha;
            if isinf(obj.alpha)
                obj.Et = [];
                obj.xt_min = [];
                obj.xt_max = [];
                obj.logWarning( ...
                    "Terminal region is unbounded; no terminal constraint is imposed.");
                return
            end
            if obj.ifPlot && n>=2 && n<=3 && ifplot
                fig = figure();
                ctr.plots.plotEllipsoid(obj.M,obj.alpha,[],fig);
            end

            % 终端约束线性化
            obj.logInfo("Terminal rigion linearized...");
            [obj.xt_min,obj.xt_max,obj.Et] = ctr.terminalSet_linear(obj.At,obj.M,obj.alpha);
            if obj.ifPlot && n>=2 && n<=3 && ifplot
                ctr.plots.plotPolytope(obj.Et,obj.xt_min,obj.xt_max,[],fig);
            end

            obj.logInfo("Terminal ingredient designing done!");
        end

        function setLogger(obj,bool)
            obj.ifLogger = bool;
        end

        function setPlot(obj,bool)
            obj.ifPlot = bool;
        end
    end

    methods(Access=private)
        function alpha = computeAlpha(obj)
            n = size(obj.A,1);
            G = [eye(n);-eye(n);obj.Kt;-obj.Kt];
            h = [obj.x_ub;-obj.x_lb;obj.u_ub;-obj.u_lb];
            len = size(G,1);
            alpha = zeros(len,1);

            if any(isnan(h)) || any(h < 0)
                obj.logError('State and input bounds must contain the origin.');
            end

            for i = 1:len
                if any(G(i,:))
                    alpha(i) = h(i)^2/(G(i,:)*(obj.M\G(i,:)'));
                else
                    alpha(i) = inf;
                end
            end
            alpha = min(alpha);
        end

        function logInfo(obj,str)
            % LOGINFO 为已格式化的消息添加时间戳、类前缀和对齐的 INFO 标签
            if ~obj.ifLogger
                return;
            end
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            fprintf('[%s] [DSN] %-9s %s\n', timestamp, '(INFO)', char(str));
        end

        function logWarning(obj,str)
            % LOGWARNING 为已格式化的消息添加时间戳、类前缀和对齐的 WARNING 标签
            if ~obj.ifLogger
                return;
            end
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            warning('[%s] [DSN] %-9s %s', timestamp, '(WARNING)', char(str));
        end
    end

    methods(Static,Access=private)
        function logError(str)
            % LOGERROR 为已格式化的消息添加时间戳、类前缀和对齐的 ERROR 标签
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            error('[%s] [DSN] %-9s %s', timestamp, '(ERROR)', char(str));
        end
    end
end
