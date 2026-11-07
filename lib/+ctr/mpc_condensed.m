classdef mpc_condensed < handle
%MPC_CONDENSED 与 ctr.mpc_sparse 使用相同控制接口的凝聚形式 MPC。
% n、m 分别为状态、输入维数；nt 为启用的终端约束行数。
% 启用终端约束时 nt = size(Et,1)，未启用时 nt = 0。

    properties(Access=private)
        % ---------- 基本属性 ----------
        A   % [nxn] 系统矩阵
        B   % [nxm] 输入矩阵

        Q   % [nxn] 状态权重矩阵
        R   % [mxm] 控制权重矩阵

        N   % 预测步数

        % ---------- 终端要素 ----------
        P       % [nxn] 终端权重矩阵
        Kt      % [mxn] 终端反馈增益，约定 u_t = Kt*x；可为空
        xt_max  % [ntx1] Et*x 的终端约束上界
        xt_min  % [ntx1] Et*x 的终端约束下界
        Et      % [ntxn] 终端约束系数矩阵；为空时不启用终端约束
        
        % ---------- 约束界限 ----------
        % （可能在线更新）
        X_lb   % [(nN)x1] 第 0 至 N-1 步状态约束下界，按步堆叠
        X_ub   % [(nN)x1] 第 0 至 N-1 步状态约束上界，按步堆叠
        U_lb   % [(mN)x1] 第 0 至 N-1 步控制约束下界，按步堆叠
        U_ub   % [(mN)x1] 第 0 至 N-1 步控制约束上界，按步堆叠
     
        % ---------- 设置标志 ----------
        hasTerminalSet = false
        ifWarmStart = true
        ifLogger = true

        % ---------- 在线存储 ----------
        % 在 update() 中被在线计算
        G_lb    % [((m+n)N+nt)x1] 增广约束下界
        G_ub    % [((m+n)N+nt)x1] 增广约束上界
        f       % [(mN)x1] 代价中的线性项系数，J = U'*H*U + 2*f'*U + c
        c       % 标量，代价中的常数项
        x_ini   % [nx1] 预测起点
        solver
        solver_name = '' % 最近一次求解使用的求解器
        solver_settings = [];
        solver_info = struct() % 最近一次求解器返回的统计信息

        % ---------- 离线存储 ----------
        % 在 mpc.init()中被离线计算
        n,m   % 维数
        QQ      % [(nN)x(nN)] 增广状态权重
        RR      % [(mN)x(mN)] 增广控制权重
        T_x     % [(nN)xn] 第 0 至 N-1 步的初始状态映射
        T_u     % [(nN)x(mN)] 第 0 至 N-1 步的控制序列映射
        Gamma_u % [nx(mN)] 第 N 步终端状态的控制序列映射
        H   % [(mN)x(mN)]
        F   % [((m+n)N+nt)x(mN)] 增广约束系数矩阵
    end


    properties(SetAccess = private, GetAccess = public)
        U               % [mxN] 第 0 至 N-1 步的最优控制序列
        solving_time    % 标量，求解时间，单位为秒
        U_f             % [mxN] 移位补项的热启动初值；可行性需额外条件保证
        X               % [nxN] 第 1 至 N 步的预测状态，不含预测起点
        J               % 标量，完整代价函数值，在 solve() 中更新
    end

    methods(Access=public)
        function obj = mpc_condensed
        end
        
        % ---------- 设置 MPC 问题 ----------
        function setProblem(obj,A,B,N,Q,R)
            % (A,B) 为离散状态方程 x(k+1) = A*x(k) + B*u(k) 的矩阵。
            obj.A = A; obj.B = B;
            obj.N = N; obj.Q = Q; obj.R = R;
            obj.resetCache();

            % 日志信息
            obj.logInfo('Successfully create MPC controller! Waiting for initialize...');
        end

        % ---------- 设置终端要素 ----------
        function setTerminal(obj,P,Kt,Et,xt_min,xt_max)
            % 终端要素由专门的终端设计器给出
            % setTerminal() 应当在 init() 之前使用
            % 若不使用，则终端要素留空，init() 中需另行处理
            obj.P = P;
            obj.Kt = Kt;
            obj.xt_max = xt_max;
            obj.xt_min = xt_min;
            obj.Et = Et;
            obj.resetCache();

            % 日志信息
            obj.logInfo('Successfully set MPC terminal ingredients! Waiting for initialize...');
        end

        % ---------- 初始化 MPC 控制器 ----------
        function init(obj)
        % INIT 进行初始化离线计算
            obj.resetCache();
            % ----- 系统维度 -----
            [obj.n,obj.m] = size(obj.B);
            
            % ----- 终端情形 -----
            P_local = obj.getTerminalWeight();
            obj.hasTerminalSet = ~isempty(obj.Et);

            % ----- 权重矩阵 -----
            obj.QQ = kron(eye(obj.N),obj.Q);
            obj.RR = kron(eye(obj.N),obj.R);

            % ----- 计算 T_u -----
            obj.T_u = zeros(obj.n*obj.N,obj.m*obj.N);
            for i = 2:obj.N
                for j = 1:i-1
                    obj.T_u((i-1)*obj.n+1:i*obj.n,(j-1)*obj.m+1:j*obj.m) = obj.A^(i-1-j)*obj.B;
                end
            end

            % ----- 计算 T_x -----
            obj.T_x = zeros(obj.n*obj.N,obj.n);
            for i = 1:obj.N
                obj.T_x((i-1)*obj.n+1:i*obj.n,:) = obj.A^(i-1);
            end

            % ----- 计算 Gamma_u -----
            obj.Gamma_u = zeros(obj.n,obj.m*obj.N);
            for j = 1:obj.N
                obj.Gamma_u(:,(j-1)*obj.m+1:j*obj.m) = obj.A^(obj.N-j)*obj.B;
            end

            % ----- 构建 H 和 F -----
            obj.H = obj.Gamma_u'*P_local*obj.Gamma_u + ...
                    obj.T_u'*obj.QQ*obj.T_u + obj.RR;
            if obj.hasTerminalSet
                obj.F = [eye(obj.m*obj.N);obj.T_u;obj.Et*obj.Gamma_u];
            else
                obj.F = [eye(obj.m*obj.N);obj.T_u];
            end

            % 清空旧有信息
            obj.U_f = [];

            % 日志信息
            obj.logInfo('MPC controller initialized! Waiting for setup...');
        end

        % ---------- 设置求解器 ----------
        function setup(obj,solver_name,settings)
        % 需要在 init() 之后使用。
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
            sol.setup(obj.H,[],obj.F,lb0,ub0,settings);

            % 初始化成功后再替换，失败时保留原有求解器及其配置。
            obj.clearResults();
            obj.solver = sol;
            obj.solver_name = solver_name;
            obj.solver_settings = settings;
        end

        % ---------- 在线更新 ----------
        function update(obj,x_ini,X_lb,X_ub,U_lb,U_ub)
        % UPDATE 进行在线更新计算
        % x_ini 为 nx1；X_lb/X_ub 为 (nN)x1；U_lb/U_ub 为 (mN)x1 列向量。
        % 更新预测起点
            obj.clearResults();
            obj.x_ini = x_ini;
            
            % 更新约束边界
            obj.X_lb = X_lb;
            obj.X_ub = X_ub;
            obj.U_lb = U_lb;
            obj.U_ub = U_ub;

            % 在线计算更新
            P_local = obj.getTerminalWeight();
            Phi = obj.T_x*obj.x_ini;
            psi = obj.A^obj.N*obj.x_ini;
            obj.f = (psi'*P_local*obj.Gamma_u + Phi'*obj.QQ*obj.T_u)';
            obj.c = psi'*P_local*psi + Phi'*obj.QQ*Phi;

            if obj.hasTerminalSet
                obj.G_lb = [obj.U_lb;obj.X_lb-Phi;obj.xt_min-obj.Et*psi];
                obj.G_ub = [obj.U_ub;obj.X_ub-Phi;obj.xt_max-obj.Et*psi];
            else
                obj.G_lb = [obj.U_lb;obj.X_lb-Phi];
                obj.G_ub = [obj.U_ub;obj.X_ub-Phi];
            end
            obj.logInfo('MPC controller updated! Waiting for solving...');
        end

        % ---------- 数值求解 ----------
        function solve(obj)
            if isempty(obj.solver) || isempty(obj.f) ...
                    || isempty(obj.G_lb) || isempty(obj.G_ub)
                obj.logError('Call .setup() and .update() before .solve().');
            end
            obj.clearResults();
            switch obj.solver_name
                case 'osqp' % ---------- 使用 OSQP 求解器 ----------
                    obj.solver.update('q',obj.f,'l',obj.G_lb,'u',obj.G_ub);
                    
                    % 开始求解：
                    if obj.ifWarmStart
                        if ~isempty(obj.U_f)
                            obj.solver.warm_start('x',obj.U_f(:));
                        end
                    else
                        % 复用求解器时，同时清零原始、对偶初值以实现冷启动。
                        obj.solver.warm_start( ...
                            'x', zeros(size(obj.H,1),1), ...
                            'y', zeros(size(obj.F,1),1));
                    end
                    result = obj.solver.solve;
                    obj.solver_info = result.info;
                    if ~strcmp(result.info.status, 'solved')
                        obj.logError(sprintf('OSQP solving not successfully: %s', result.info.status));
                        % obj.logWarning(sprintf('OSQP solving not successfully: %s', result.info.status));
                    end
                    u = result.x;
                    obj.J = 2*result.info.obj_val + obj.c;
                    obj.solving_time = result.info.solve_time;
                        % 记录的是求解时间，而非整个控制决策用时
                        % 后者还包括状态观测、变量传入等耗时
                % case 'SOLVER_NAME' % ---------- 使用 SOLVER_NAME 求解器 ----------
                otherwise
                    obj.logError(sprintf('No solver named ''%s''!', obj.solver_name));
            end

            % 返回前的处理
            x = reshape(obj.T_x*obj.x_ini + obj.T_u*u,obj.n,obj.N);
            xt = obj.A^obj.N*obj.x_ini + obj.Gamma_u*u;
            obj.X = [x(:,2:end),xt];
            obj.U = reshape(u,obj.m,obj.N);
            if ~isempty(obj.Kt)
                obj.U_f = [obj.U(:,2:end),obj.Kt*xt];
            else
                obj.U_f = [obj.U(:,2:end),obj.U(:,end)];
            end

            % 日志信息
            obj.logInfo(sprintf('Soving done, took %.1f ms.',1000*obj.solving_time));
        end

        function setLogger(obj,bool)
            obj.ifLogger = bool;
        end

        function setWarmStart(obj,bool)
            obj.ifWarmStart = bool;
        end
        
        function setWarmStartGuess(obj,U)
        % 接受 mxN 输入序列（或其列展开），在 init() 后使用。
            if isempty(U)
                obj.U_f = [];
            elseif numel(U) == obj.m*obj.N
                obj.U_f = reshape(U,obj.m,obj.N);
            else
                obj.logError('Warm-start guess must be an input sequence.');
            end
        end

        function info = viewMPC(obj,verbose)
        % VIEWMPC 查看 MPC 配置、约束及最近一次求解结果，不受 ifLogger 控制。
        % obj.viewMPC() 显示摘要；obj.viewMPC(true) 同时显示完整矩阵。
        % info = obj.viewMPC() 另外返回分组的原始数据结构体。
            if nargin < 2
                verbose = false;
            end
            validateattributes(verbose,{'logical','numeric'},{'scalar','binary'},mfilename,'verbose');

            n = []; m = []; decisionCount = [];
            if ~isempty(obj.B)
                [n,m] = size(obj.B);
            end
            if ~isempty(m) && ~isempty(obj.N)
                decisionCount = m*obj.N;
            end
            configured = ~isempty(obj.A) && ~isempty(obj.B) && ~isempty(obj.N) ...
                         && ~isempty(obj.Q) && ~isempty(obj.R);
            initialized = ~isempty(obj.H) && ~isempty(obj.F);
            solverReady = initialized && ~isempty(obj.solver);
            ready = solverReady && ~isempty(obj.f) && ~isempty(obj.G_lb) && ~isempty(obj.G_ub);
            solved = ~isempty(obj.U) && ~isempty(obj.J);
            stages = {'Not configured','Awaiting init','Awaiting setup', ...
                      'Awaiting update','Ready to solve','Solved'};
            if solved
                stage = stages{6};
            elseif ready
                stage = stages{5};
            elseif solverReady
                stage = stages{4};
            elseif initialized
                stage = stages{3};
            elseif configured
                stage = stages{2};
            else
                stage = stages{1};
            end
            firstControl = [];
            if solved
                firstControl = obj.U(:,1);
            end
            constraintCount = [];
            if initialized
                constraintCount = size(obj.F,1);
            end

            info.state = struct('stage',stage,'configured',configured, ...
                'initialized',initialized,'ready',ready,'solved',solved);
            info.model = struct('n',n,'m',m,'N',obj.N, ...
                'A',obj.A,'B',obj.B); %#ok<*PROPLC>
            P_local = obj.getTerminalWeight();
            info.weights = struct('Q',obj.Q,'R',obj.R,'P',P_local);
            info.terminal = struct('enabled',~isempty(obj.Et),'rows',size(obj.Et,1), ...
                'Et',obj.Et,'lower',obj.xt_min,'upper',obj.xt_max,'Kt',obj.Kt);
            info.bounds = struct('X_lb',obj.X_lb,'X_ub',obj.X_ub, ...
                'U_lb',obj.U_lb,'U_ub',obj.U_ub);
            info.qp = struct('variables',decisionCount,'constraint_rows',constraintCount, ...
                'QQ',obj.QQ,'RR',obj.RR,'T_x',obj.T_x,'T_u',obj.T_u, ...
                'Gamma_u',obj.Gamma_u,'H',obj.H,'F',obj.F, ...
                'f',obj.f,'c',obj.c,'G_lb',obj.G_lb,'G_ub',obj.G_ub);
            info.options = struct('logger',obj.ifLogger,'warm_start',obj.ifWarmStart, ...
                'warm_start_available',~isempty(obj.U_f));
            info.result = struct('solver',obj.solver_name,'statistics',obj.solver_info, ...
                'x_ini',obj.x_ini,'cost',obj.J,'solve_time_s',obj.solving_time, ...
                'solve_time_ms',1000*obj.solving_time,'first_control',firstControl, ...
                'U',obj.U,'X',obj.X,'U_f',obj.U_f);

            fprintf('\n=========== Condensed MPC overview ===========\n');
            obj.printViewSection('State',{'Stage',stage});
            obj.printViewSection('Model',{'n (states)',n;'m (inputs)',m; ...
                'N (horizon)',obj.N;'A',obj.A;'B',obj.B});
            obj.printViewSection('Weights',{'Q',obj.Q;'R',obj.R;'P',P_local});
            obj.printViewSection('Terminal',{'Constraint enabled',~isempty(obj.Et); ...
                'Constraint rows',size(obj.Et,1);'Et',obj.Et; ...
                'xt_min',obj.xt_min;'xt_max',obj.xt_max;'Kt (u_t = Kt*x)',obj.Kt});
            obj.printViewSection('Bounds (stacked, steps 0:N-1)', ...
                {'X_lb',obj.boundSummary(obj.X_lb);'X_ub',obj.boundSummary(obj.X_ub); ...
                 'U_lb',obj.boundSummary(obj.U_lb);'U_ub',obj.boundSummary(obj.U_ub)});
            obj.printViewSection('QP',{'Decision variables',decisionCount; ...
                'Constraint rows',constraintCount;'H size',mat2str(size(obj.H)); ...
                'F size',mat2str(size(obj.F));'f size',mat2str(size(obj.f))});
            obj.printViewSection('Options',{'Logger',obj.ifLogger;'Warm start',obj.ifWarmStart; ...
                'Warm start available',~isempty(obj.U_f)});
            obj.printViewSection('Latest solve',{'Solver',obj.solver_name; ...
                'Initial state',obj.x_ini;'Full cost J',obj.J; ...
                'Solve time (ms)',1000*obj.solving_time;'First control',firstControl; ...
                'U size (steps 0:N-1)',mat2str(size(obj.U)); ...
                'X size (steps 1:N)',mat2str(size(obj.X)); ...
                'U_f size',mat2str(size(obj.U_f))});
            if ~isempty(fieldnames(obj.solver_info))
                names = {'status','iter','setup_time','solve_time','update_time','polish_time','run_time'};
                rows = cell(0,2);
                for i = 1:numel(names)
                    if isfield(obj.solver_info,names{i})
                        rows(end+1,:) = {names{i},obj.solver_info.(names{i})}; %#ok<AGROW>
                    end
                end
                obj.printViewSection('Solver statistics (times in seconds)',rows);
            end
            if verbose
                groups = {'model','weights','terminal','bounds','qp','result'};
                fprintf('\n--- Full data ---\n');
                for i = 1:numel(groups)
                    data = info.(groups{i});
                    names = fieldnames(data);
                    for j = 1:numel(names)
                        fprintf('\n%s.%s =\n',groups{i},names{j});
                        disp(data.(names{j}));
                    end
                end
            end
            fprintf('==============================================\n');
        end

    end

    methods(Access=private)
        function P = getTerminalWeight(obj)
            % 默认零权重按当前维数生成，不写回用户配置。
            P = obj.P;
            if isempty(P) && ~isempty(obj.n)
                P = zeros(obj.n);
            end
        end

        function clearResults(obj)
            obj.U = []; obj.X = []; obj.J = []; obj.solving_time = [];
            obj.solver_info = struct();
        end

        function resetCache(obj)
            % 保留配置，清除依赖配置的计算数据和历史结果。
            names = {'n','m','QQ','RR','T_x','T_u','Gamma_u','H','F', ...
                'G_lb','G_ub','f','c','x_ini','X_lb','X_ub','U_lb','U_ub','U_f'};
            for i = 1:numel(names)
                obj.(names{i}) = [];
            end
            obj.solver = [];
            obj.solver_name = '';
            obj.solver_settings = [];
            obj.hasTerminalSet = false;
            obj.clearResults();
        end

        function logInfo(obj,str)
            % LOGINFO 为已格式化的消息添加时间戳和对齐的 INFO 标签
            if ~obj.ifLogger
                return;
            end
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            fprintf('[%s] [MPC] %-9s %s\n', timestamp, '(INFO)', char(str));
        end

        function logWarning(obj,str)
            % LOGWARNING 为已格式化的消息添加时间戳和对齐的 WARNING 标签
            if ~obj.ifLogger
                return;
            end
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            warning('[%s] [MPC] %-9s %s', timestamp, '(WARNING)', char(str));
        end
    end

    methods(Static,Access=private)
        function printViewSection(title,rows)
            fprintf('\n--- %s ---\n',title);
            for i = 1:size(rows,1)
                value = rows{i,2};
                if isempty(value)
                    summary = 'Not set / no data';
                elseif ischar(value) || isstring(value)
                    summary = char(value);
                elseif islogical(value) && isscalar(value)
                    labels = {'Off','On'};
                    summary = labels{double(value)+1};
                elseif numel(value) <= 9
                    summary = mat2str(full(value),6);
                else
                    summary = sprintf('%d x %d matrix',size(value,1),size(value,2));
                end
                fprintf('  %-29s : %s\n',rows{i,1},summary);
            end
        end

        function summary = boundSummary(value)
            if isempty(value)
                summary = 'Not set';
                return;
            end
            summary = sprintf('%d x %d, min=%.6g, max=%.6g, Inf=%d, NaN=%d', ...
                size(value,1),size(value,2),min(value(:)),max(value(:)), ...
                nnz(isinf(value)),nnz(isnan(value)));
        end

        function logError(str)
            % LOGERROR 为已格式化的消息添加时间戳和对齐的 ERROR 标签
            timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
            error('[%s] [MPC] %-9s %s', timestamp, '(ERROR)', char(str));
        end
    end
end
