classdef deepc < handle
    properties(Access=private)
        % ========== 数据信息 ==========
        data_u      % 输入数据
        data_y      % 输出数据
        Hankel_u    % 输入数据构成的 Hankel 矩阵
        Hankel_y    % 输出数据构成的 Hankel 矩阵
        m,p         % 输入输出维度
        S           % 数据样本数目

        % ========== 参数信息 ==========
        N           % 预测步数
        T_ini       % 初始数据

        Q               % 跟踪代价
        R               % 控制代价
        lambda_u = 1e6    % 输入松弛项正则化系数
        lambda_y = 1e4    % 输出松弛项正则化系数
        lambda_a = 1e2    % alpha项松弛正则化系数

        ratio=10;   % 拐点阈值比例

        % ========== 中间变量 ==========
        QQ      % 增广输出权重矩阵
        RR      % 增广输入权重矩阵

        H       % 二次项系数矩阵
        f       % 一次项系数矩阵
        c       % 常数项矩阵（实时更新）
        F       % 约束系数矩阵（实时更新）
        G_lb    % 约束下边界
        G_ub    % 约束上边界
        
        solver
        solver_name = ''        % 最近一次求解使用的求解器
        solver_settings = [];
        solver_info = struct()  % 最近一次求解器返回的统计信息
        solving_time

        % ========== 状态标志 ==========
        ifWarmStart = true  % 热启动（上一时刻解的移位+保持）
        ifLogger = true     % 打印日志
    end

    properties(SetAccess=private)
        % ========== 优化求解 ==========
        alpha       % 表征系数
        u_series    % [mxN] 预测输入序列
        y_series    % [pxN] 预测输出序列
        slack_u     % 关于输入数据的松弛变量
        slack_y     % 关于输出数据的松弛变量
        U           % [alpha;u;y;slack_u;slack_y] 原始解
        U_f         % 完整决策向量的热启动初值
        J           % 代价函数

    end

    methods(Access=public)
        function obj = deepc(N,T_ini)
            validateattributes(N,{'numeric'},{'scalar','real','finite','integer','positive'},mfilename,'N');
            validateattributes(T_ini,{'numeric'},{'scalar','real','finite','integer','positive'},mfilename,'T_ini');
            obj.N = double(N);
            obj.T_ini = double(T_ini);
        end

        function setParams(obj,Q,R,lambda_u,lambda_y,lambda_a)
            Q = obj.checkWeight(Q,'Q');
            R = obj.checkWeight(R,'R');
            validateattributes(lambda_u,{'numeric'},{'scalar','real','finite','positive'},mfilename,'lambda_u');
            validateattributes(lambda_y,{'numeric'},{'scalar','real','finite','positive'},mfilename,'lambda_y');
            validateattributes(lambda_a,{'numeric'},{'scalar','real','finite','nonnegative'},mfilename,'lambda_a');
            if ~isempty(obj.Hankel_u) && ...
                    (~isequal(size(Q),[obj.p,obj.p]) || ~isequal(size(R),[obj.m,obj.m]))
                obj.logError('Q and R must match the output and input dimensions.');
            end
            obj.Q = Q;
            obj.R = R;
            obj.lambda_u = double(lambda_u);
            obj.lambda_y = double(lambda_y);
            obj.lambda_a = double(lambda_a);
            if ~isempty(obj.Hankel_u)
                obj.buildCost();
                obj.resetSolver();
            end
        end

        function init(obj,u_d,y_d)
            % 理论保证另需输入数据满足 (N+T_ini+n) 阶 PE，且 T_ini >= 系统 lag。
            validateattributes(u_d,{'numeric'},{'2d','nonempty','real','finite'},mfilename,'u_d');
            validateattributes(y_d,{'numeric'},{'2d','nonempty','real','finite'},mfilename,'y_d');
            if size(u_d,2) ~= size(y_d,2)
                obj.logError('u_d and y_d must contain the same number of time-aligned samples.');
            end
            if isempty(obj.Q) || isempty(obj.R)
                obj.logError('Call .setParams() before .init().');
            end
            if ~isequal(size(obj.Q),[size(y_d,1),size(y_d,1)]) || ...
                    ~isequal(size(obj.R),[size(u_d,1),size(u_d,1)])
                obj.logError('Q and R must match the output and input dimensions.');
            end
            Hu = ctr.deepc.HankelMat(obj.N+obj.T_ini,u_d);
            Hy = ctr.deepc.HankelMat(obj.N+obj.T_ini,y_d);
            obj.data_u = u_d;
            obj.data_y = y_d;
            obj.Hankel_u = Hu;
            obj.Hankel_y = Hy;
            
            [obj.m,obj.S] = size(u_d);
            [obj.p,~] = size(y_d);

            obj.buildCost();
            obj.buildF();
            obj.resetSolver();
        end

        function r = samplify(obj,ratio,ifplot)
            if isempty(obj.Hankel_u)
                obj.logError('Call .init() before .samplify().');
            end
            if nargin < 2
                ratio = obj.ratio;
            end
            if nargin < 3
                ifplot = false;
            end
            validateattributes(ratio,{'numeric'},{'scalar','real','finite','>',1},mfilename,'ratio');
            validateattributes(ifplot,{'logical','numeric'},{'scalar','binary'},mfilename,'ifplot');
            [W,Sigma,~] = svd([obj.Hankel_u;obj.Hankel_y],'econ');
            sigma = diag(Sigma)+eps;
            r = length(sigma);
            for i = 1:length(sigma)-1
                if sigma(i)/sigma(i+1) >= ratio
                    r = i;
                end
            end
            Hkl = W(:,1:r)*Sigma(1:r,1:r);
            obj.Hankel_u = Hkl(1:obj.m*(obj.N+obj.T_ini),1:r);
            obj.Hankel_y = Hkl(obj.m*(obj.N+obj.T_ini)+1:(obj.m+obj.p)*(obj.N+obj.T_ini),1:r);
            obj.buildCost();
            obj.buildF();
            obj.resetSolver();

            if ifplot
                ctr.plots.plotSVD(sigma);
            end
        end

        function setup(obj,solver_name,settings)
            % 需要在 init() 之后使用
            if isempty(obj.H) || isempty(obj.F)
                obj.logError('Call .init() before .setup().');
            end
            solver_name = char(solver_name);
            if ~strcmp(solver_name,'osqp')
                obj.logError(sprintf('No solver named ''%s''!', solver_name));
            end
            try sol = osqp;
            catch 
                obj.logError('OSQP not installed!');
            end
            lb0 = -inf(size(obj.F,1),1);
            ub0 =  inf(size(obj.F,1),1);
            % OSQP 求解原代价的一半：0.5*U'*H*U + f*U。
            sol.setup(obj.H,[],obj.F,lb0,ub0,settings);

            % 初始化成功后再替换，失败时保留原有求解器及其配置。
            obj.clearResults();
            obj.U_f = [];
            obj.solver = sol;
            obj.solver_name = solver_name;
            obj.solver_settings = settings;
        end

        function update(obj,u_ini,y_ini,u_r,y_r,Y_lb,Y_ub,U_lb,U_ub)
            if isempty(obj.H) || isempty(obj.F)
                obj.logError('Call .init() before .update().');
            end
            % 接受按时间堆叠的行/列向量，内部统一为列向量。
            u_ini = obj.checkVector(u_ini,obj.m*obj.T_ini,'u_ini',false);
            y_ini = obj.checkVector(y_ini,obj.p*obj.T_ini,'y_ini',false);
            u_r = obj.checkVector(u_r,obj.m*obj.N,'u_r',false);
            y_r = obj.checkVector(y_r,obj.p*obj.N,'y_r',false);
            Y_lb = obj.checkVector(Y_lb,obj.p*obj.N,'Y_lb',true);
            Y_ub = obj.checkVector(Y_ub,obj.p*obj.N,'Y_ub',true);
            U_lb = obj.checkVector(U_lb,obj.m*obj.N,'U_lb',true);
            U_ub = obj.checkVector(U_ub,obj.m*obj.N,'U_ub',true);
            if any(Y_lb > Y_ub) || any(U_lb > U_ub)
                obj.logError('Every lower bound must be less than or equal to its upper bound.');
            end
            obj.clearResults();
            obj.f = [zeros(1,size(obj.Hankel_u,2)),...
                     -u_r'*obj.RR, -y_r'*obj.QQ,...
                     zeros(1,obj.m*obj.T_ini),...
                     zeros(1,obj.p*obj.T_ini)];
            obj.c = u_r'*obj.RR*u_r + y_r'*obj.QQ*y_r;
            obj.G_lb = [u_ini;zeros(obj.m*obj.N,1);...
                        y_ini;zeros(obj.p*obj.N,1);...
                        U_lb;Y_lb];
            obj.G_ub = [u_ini;zeros(obj.m*obj.N,1);...
                        y_ini;zeros(obj.p*obj.N,1);...
                        U_ub;Y_ub];
            if ~isempty(obj.U_f)
                alpha_guess = obj.U_f(1:size(obj.Hankel_u,2));
                % 历史窗口已更新，用当前历史信息恢复松弛变量初值。
                obj.U_f(end-(obj.m+obj.p)*obj.T_ini+1:end) = ...
                    [obj.Hankel_u(1:obj.m*obj.T_ini,:)*alpha_guess-u_ini;...
                     obj.Hankel_y(1:obj.p*obj.T_ini,:)*alpha_guess-y_ini];
            end
        end

        function solve(obj)
            if isempty(obj.solver) || isempty(obj.f) ...
                    || isempty(obj.G_lb) || isempty(obj.G_ub)
                obj.logError('Call .setup() and .update() before .solve().');
            end
            obj.clearResults();

            switch obj.solver_name
                case 'osqp'
                    obj.solver.update('q',obj.f(:),'l',obj.G_lb,'u',obj.G_ub);
                    if obj.ifWarmStart
                        if ~isempty(obj.U_f)
                            obj.solver.warm_start('x',obj.U_f);
                        end
                    else
                        obj.solver.warm_start( ...
                            'x',zeros(size(obj.H,1),1), ...
                            'y',zeros(size(obj.F,1),1));
                    end
                    result = obj.solver.solve;
                    obj.solver_info = result.info;
                    if ~strcmp(result.info.status,'solved')
                        obj.logError(sprintf('OSQP solving not successfully: %s',result.info.status));
                    end
                    obj.U = result.x;
                    obj.J = 2*result.info.obj_val + obj.c;
                    obj.solving_time = result.info.solve_time;
                otherwise
                    obj.logError(sprintf('No solver named ''%s''!',obj.solver_name));
            end

            n_alpha = size(obj.Hankel_u,2);
            obj.alpha = obj.U(1:n_alpha);
            obj.u_series = reshape(obj.U(n_alpha+(1:obj.m*obj.N)),obj.m,obj.N);
            obj.y_series = reshape(obj.U(n_alpha+obj.m*obj.N+(1:obj.p*obj.N)),obj.p,obj.N);
            obj.slack_u = obj.U(n_alpha+(obj.m+obj.p)*obj.N+(1:obj.m*obj.T_ini));
            obj.slack_y = obj.U(n_alpha+(obj.m+obj.p)*obj.N+obj.m*obj.T_ini+(1:obj.p*obj.T_ini));

            % alpha 是数据基的坐标，不按时间移位；输入输出序列移位并保持末项。
            % 热启动仅提供初值，不要求它满足下一轮的新约束。
            u_guess = [obj.u_series(:,2:end),obj.u_series(:,end)];
            y_guess = [obj.y_series(:,2:end),obj.y_series(:,end)];
            obj.U_f = [obj.alpha;u_guess(:);y_guess(:);obj.slack_u;obj.slack_y];

            obj.logInfo(sprintf('Solving done, took %.1f ms.',1000*obj.solving_time));
        end

        function setLogger(obj,bool)
            obj.ifLogger = bool;
        end

        function setWarmStart(obj,bool)
            obj.ifWarmStart = bool;
        end
    end

    methods(Static,Access=public)
        function [H,r,row] = HankelMat(order,signal)
            % signal 的每列对应一个采样时刻；单通道信号应传入行向量。
            validateattributes(signal, {'numeric'}, ...
                {'2d', 'nonempty', 'real', 'finite'}, mfilename, 'signal', 2);
            [m,T] = size(signal);   % [输入维数,数据长度]
            validateattributes(order, {'numeric'}, ...
                {'scalar', 'real', 'finite', 'integer', 'positive', '<=', T}, ...
                mfilename, 'order', 1);
            % 统一为 double，避免整数类型参与维度计算时发生饱和或类型冲突。
            order = double(order);
            window_length = T-order+1;
            % 构造 Hankel 矩阵不要求满行秩；持续激励条件应另行检查输入数据。
            H = zeros(m*order,window_length);
            for i = 1:window_length
                for j = 1:order
                    H((j-1)*m+1:j*m,i) = signal(:,j+i-1);
                end
            end
            if nargout > 1
                r = rank(H);
            end
            row = m*order;
        end
    end

    methods(Access=private)
        function W = checkWeight(obj,W,name)
            validateattributes(W,{'numeric'},{'2d','nonempty','real','finite','square'},mfilename,name);
            W = double(W);
            tol = 1e-12*norm(W,'fro');
            if norm(W-W','fro') > tol
                obj.logError(sprintf('%s must be symmetric.',name));
            end
            W = (W+W')/2;
            if min(eig(W)) < -tol
                obj.logError(sprintf('%s must be positive semidefinite.',name));
            end
        end

        function value = checkVector(~,value,count,name,allowInf)
            validateattributes(value,{'numeric'},{'vector','real','numel',count},mfilename,name);
            if any(isnan(value(:))) || (~allowInf && any(isinf(value(:))))
                error('deepc:InvalidVector','%s contains invalid nonfinite values.',name);
            end
            value = double(value(:));
        end

        function buildCost(obj)
            obj.RR = kron(eye(obj.N),obj.R);
            obj.QQ = kron(eye(obj.N),obj.Q);
            obj.H = blkdiag(obj.lambda_a*eye(size(obj.Hankel_u,2)),...
                            obj.RR,obj.QQ,obj.lambda_u*eye(obj.m*obj.T_ini),...
                            obj.lambda_y*eye(obj.p*obj.T_ini));
        end

        function clearResults(obj)
            obj.alpha = [];
            obj.u_series = [];
            obj.y_series = [];
            obj.slack_u = [];
            obj.slack_y = [];
            obj.U = [];
            obj.J = [];
            obj.solving_time = [];
            obj.solver_info = struct();
        end

        function resetSolver(obj)
            % 数据基变化后，旧求解器和热启动向量不再适用。
            obj.solver = [];
            obj.solver_name = '';
            obj.solver_settings = [];
            obj.U_f = [];
            obj.f = [];
            obj.c = [];
            obj.G_lb = [];
            obj.G_ub = [];
            obj.clearResults();
        end

        function logInfo(obj,str)
            if ~obj.ifLogger
                return;
            end
            timestamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
            fprintf('[%s] [DeePC] %-9s %s\n',timestamp,'(INFO)',char(str));
        end

        function logError(~,str)
            timestamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
            error('[%s] [DeePC] %-9s %s',timestamp,'(ERROR)',char(str));
        end

        function buildF(obj)
            n_alpha = size(obj.Hankel_u,2);
            obj.F = zeros((obj.m+obj.p)*(2*obj.N+obj.T_ini),...
                          n_alpha+(obj.m+obj.p)*(obj.N+obj.T_ini));

            % 分块行索引
            r1 = 1:obj.m*obj.T_ini;
            r2 = obj.m*obj.T_ini+(1:obj.m*obj.N);
            r3 = obj.m*(obj.T_ini+obj.N)+(1:obj.p*obj.T_ini);
            r4 = (obj.m+obj.p)*obj.T_ini+obj.m*obj.N+(1:obj.p*obj.N);
            r5 = (obj.m+obj.p)*(obj.T_ini+obj.N)+(1:obj.m*obj.N);
            r6 = (obj.m+obj.p)*(obj.T_ini+obj.N)+obj.m*obj.N+(1:obj.p*obj.N);

            % 分块列索引：alpha、u、y、sigma_u、sigma_y
            cg = 1:n_alpha;
            cu = n_alpha+(1:obj.m*obj.N);
            cy = n_alpha+obj.m*obj.N+(1:obj.p*obj.N);
            csu = n_alpha+(obj.m+obj.p)*obj.N+(1:obj.m*obj.T_ini);
            csy = n_alpha+(obj.m+obj.p)*obj.N+obj.m*obj.T_ini+(1:obj.p*obj.T_ini);

            % Hankel 块
            obj.F(r1,cg) = obj.Hankel_u(1:obj.m*obj.T_ini,:);
            obj.F(r2,cg) = obj.Hankel_u(obj.m*obj.T_ini+1:end,:);
            obj.F(r3,cg) = obj.Hankel_y(1:obj.p*obj.T_ini,:);
            obj.F(r4,cg) = obj.Hankel_y(obj.p*obj.T_ini+1:end,:);

            % 单位矩阵块
            obj.F(r1,csu) = -eye(obj.m*obj.T_ini);
            obj.F(r2,cu)  = -eye(obj.m*obj.N);
            obj.F(r3,csy) = -eye(obj.p*obj.T_ini);
            obj.F(r4,cy)  = -eye(obj.p*obj.N);
            obj.F(r5,cu)  =  eye(obj.m*obj.N);
            obj.F(r6,cy)  =  eye(obj.p*obj.N);
        end
    end
end
