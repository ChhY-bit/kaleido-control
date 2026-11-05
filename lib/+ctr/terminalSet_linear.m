function [xt_min,xt_max,Et] = terminalSet_linear(At,M,alpha,beta)
%TERMINALSET_LINEAR 从收缩椭球构造线性闭环系统的多面体终端集。
%   [xt_min,xt_max,Et] = ctr.terminalSet_linear(At,M,alpha)
%   [xt_min,xt_max,Et] = ctr.terminalSet_linear(At,M,alpha,beta)
%   输入、输出顺序与 ctr.terminalSet 相同，返回集合
%       P = {x : xt_min <= Et*x <= xt_max}。
%
%   At    - n x n 有限实数闭环矩阵。
%   M     - n x n 对称正定 Lyapunov 矩阵，须满足
%           At'*M*At < M；仅 At 为 Schur 稳定还不够。
%   alpha - 非负有限实数，原椭球为 E = {x : x'*M*x <= alpha}。
%   beta  - 可选的有限实数标量，须满足 beta_min <= beta < 1，
%           beta_min = norm(R*At/R,2)，其中 M = R'*R。
%           缺省时采用最小可行值 beta_min；增大 beta 可提高几何逼近
%           精度，通常也增加约束数量和离线计算量。
%
%   对 alpha > 0，构造 beta*E 包含于 P、P 包含于 E 的凸多面体。
%   系统满足 At*E 包含于 beta_min*E，且 beta_min <= beta，因此 At*P
%   包含于 P。该方法保守地内逼近 E，不保证最大不变集。
%   alpha = 0 时返回原点。支持奇异 At；一维直接返回区间。
%
%   实现 Alessio et al., Automatica 43 (2007), 2096-2103 的
%   Algorithm 1：单纯形仿射插值、QP 误差检查、最长边二分、子水平集
%   提取和凸包。归一化为单位球后利用坐标反射对称性，仅细分第一卦限。
%   n >= 2 且 alpha > 0 时需要 Optimization Toolbox 的 quadprog。
%   适合低维离线计算；维数大或 beta 接近 1 时可能触发资源上限并报错，
%   不会返回未经认证的近似集。输出 xt_min、xt_max 均为列向量。
%
%   示例：
%       [lb,ub,E] = ctr.terminalSet_linear([.7 .1;0 .6],[2 .3;.3 1],1);
%       ctr.plots.plotPolytope(E,lb,ub);

    narginchk(3,4);
    validateattributes(At,{'numeric'}, ...
        {'real','finite','2d','nonempty','square'},mfilename,'At');
    n = size(At,1);
    validateattributes(M,{'numeric'}, ...
        {'real','finite','size',[n n]},mfilename,'M');
    validateattributes(alpha,{'numeric'}, ...
        {'real','finite','scalar','nonnegative'},mfilename,'alpha');
    At = full(double(At));
    M = full(double(M));
    alpha = double(alpha);
    tol = 1e-10;
    if norm(M-M','fro') > tol*norm(M,'fro')
        error('terminalSet_linear:NotSymmetric','M 必须为对称矩阵。');
    end
    M = M/2+M'/2;
    [R,flag] = chol(M);
    if flag ~= 0
        error('terminalSet_linear:NotPositiveDefinite','M 必须为正定矩阵。');
    end
    beta_min = norm(R*At/R,2);
    if ~isfinite(beta_min) || beta_min >= 1
        error('terminalSet_linear:NotContractive', ...
            'M 不能证明一步严格收缩：norm(R*At/R,2) = %.16g，须小于 1。',beta_min);
    end
    if nargin < 4
        beta = beta_min;
    else
        validBeta = isnumeric(beta) && isreal(beta) && isscalar(beta) && isfinite(beta);
        if validBeta
            % 先统一精度，再检查边界，避免 single 比较时舍入 beta_min。
            beta = full(double(beta));
            validBeta = beta >= beta_min && beta < 1;
        end
        if ~validBeta
            error('terminalSet_linear:InvalidBeta', ...
                'beta 必须为有限实数标量，且处于区间 [%.17g, 1)。',beta_min);
        end
    end
    if alpha == 0
        Et = eye(n);
        xt_min = zeros(n,1);
        xt_max = zeros(n,1);
        return
    elseif n == 1
        Et = R;
        xt_max = sqrt(alpha);
        xt_min = -xt_max;
        return
    end
    delta = (1-beta)*(1+beta);
    if delta <= 100*tol
        error('terminalSet_linear:InsufficientMargin', ...
            ['beta = %.17g 的收缩裕度过小，无法在当前数值容差下可靠构造多面体。' ...
             'beta 的理论可行区间为 [%.17g, 1)，请减小 beta 或调整 M。'], ...
            beta,beta_min);
    end
    if exist('quadprog','file') ~= 2
        error('terminalSet_linear:MissingSolver', ...
            '此函数需要 Optimization Toolbox 的 quadprog。');
    end

    maxSimplices = 100000;
    maxVertices = 200000;
    if n > floor(log2(maxVertices))
        error('terminalSet_linear:ResourceLimit','维数过大，坐标反射超出顶点上限。');
    end
    options = optimoptions('quadprog','Display','off', ...
        'Algorithm','interior-point-convex','ConstraintTolerance',tol, ...
        'OptimalityTolerance',tol,'MaxIterations',200);

    % y = R*x/sqrt(alpha)，q(y) = y'*y。初始正卦限单纯形的所有反射
    % 组成 sqrt(n)*单位 l1 球，包含单位 l2 球，无须额外的 Delaunay 剖分。
    pending = {[zeros(1,n);sqrt(n)*eye(n)]};
    pieces = cell(0,1);
    vertexCount = 0;
    processed = 0;
    % 保留少量严格裕度，供凸包转换和浮点舍入认证使用。
    gap = max(10*tol,0.01*delta);
    while ~isempty(pending)
        processed = processed+1;
        if processed > maxSimplices
            error('terminalSet_linear:ResourceLimit', ...
                '单纯形数超过 %d；beta = %.6g，请使用更强收缩的终端椭球。', ...
                maxSimplices,beta);
        end
        V = pending{end};
        pending(end) = [];
        values = sum(V.^2,2);

        % 在重心坐标 lambda 中解论文式 (2)：
        % min ||V'*lambda||^2 - values'*lambda + delta，
        % s.t. lambda >= 0, sum(lambda) = 1。
        H = 2*(V*V');
        H = H/2+H'/2;
        [weights,~,exitflag] = quadprog(H,-values,[],[], ...
            ones(1,n+1),1,zeros(n+1,1),[],[],options);
        if exitflag <= 0 || isempty(weights) || any(~isfinite(weights))
            error('terminalSet_linear:QPFailure', ...
                '单纯形误差 QP 求解失败，quadprog exitflag = %d。',exitflag);
        end
        point = V'*weights;
        % 凸二次函数的切平面给出最优值的下界；避免把近似求解得到的
        % 目标值（最优值上界）误当作整个单纯形上的误差证书。
        lowerBound = delta-point'*point+min(2*V*point-values);
        if lowerBound > gap
            clipped = clipSimplex(V,values);
            if ~isempty(clipped)
                vertexCount = vertexCount+size(clipped,1);
                if vertexCount > maxVertices
                    error('terminalSet_linear:ResourceLimit','子水平集顶点数超出上限。');
                end
                pieces{end+1,1} = clipped; %#ok<AGROW>
            end
        else
            % 将最长边中点作为新顶点，得到覆盖原单纯形的两个子单纯形。
            longest = -Inf;
            for i = 1:n
                for j = i+1:n+1
                    distance = sum((V(i,:)-V(j,:)).^2);
                    if distance > longest
                        longest = distance;
                        first = i;
                        second = j;
                    end
                end
            end
            midpoint = (V(first,:)+V(second,:))/2;
            if isequal(midpoint,V(first,:)) || isequal(midpoint,V(second,:))
                error('terminalSet_linear:NumericalFailure','单纯形二分已达到浮点精度极限。');
            end
            left = V;
            right = V;
            left(first,:) = midpoint;
            right(second,:) = midpoint;
            pending(end+1:end+2) = {left,right};
        end
    end

    positive = uniquetol(vertcat(pieces{:}),tol,'ByRows',true,'DataScale',1);
    % 先取正卦限部分的凸包，去掉内部网格点，再反射，控制几何计算规模。
    facets = convhulln(positive,{'Qt','Qx'});
    positive = positive(unique(facets(:)),:);
    reflectionCount = 2^n;
    if reflectionCount*size(positive,1) > maxVertices
        error('terminalSet_linear:ResourceLimit','反射后的顶点数超出上限。');
    end
    signs = 2*(dec2bin(0:reflectionCount-1,n)-'0')-1;
    vertices = zeros(reflectionCount*size(positive,1),n);
    for i = 1:reflectionCount
        rows = (i-1)*size(positive,1)+(1:size(positive,1));
        vertices(rows,:) = positive.*signs(i,:);
    end
    vertices = uniquetol(vertices,tol,'ByRows',true,'DataScale',1);
    A = hullHalfspaces(vertices);

    % 中心对称性允许用 |A*y| <= 1 表示：每对相反法向只保留一个。
    [~,pivot] = max(abs(A),[],2);
    orientation = A(sub2ind(size(A),(1:size(A,1))',pivot));
    A = A.*sign(orientation);
    A = uniquetol(A,tol,'ByRows',true,'DataScale',1);

    % 对最终半空间表示再认证，而不只检查构造时的采样顶点：
    % P 的极集为 conv([A;-A])，极集各面的法向即 P 的顶点。
    % 这也检查合并共面约束时的数值误差，必要时将 P 微量缩小。
    actualVertices = hullHalfspaces([A;-A]);
    radius = max(vecnorm(actualVertices,2,2));
    if ~isfinite(radius) || radius <= 0
        error('terminalSet_linear:NumericalFailure','无法认证多面体的外包半径。');
    end
    A = A*max(1,radius/(1-10*tol));
    if beta*max(vecnorm(A,2,2)) >= 1-10*tol
        error('terminalSet_linear:NumericalFailure', ...
            '数值处理后的多面体无法可靠包含 beta*E，未返回终端约束。');
    end

    Et = A*R;
    if any(~isfinite(Et(:)))
        error('terminalSet_linear:NumericalFailure','终端约束矩阵超出浮点数范围。');
    end
    xt_max = sqrt(alpha)*ones(size(Et,1),1);
    xt_min = -xt_max;
end

function points = clipSimplex(V,values)
% 仿射插值 <= 1 的截集顶点：保留的原顶点及平面与单纯形边的交点。
    points = V(values <= 1,:);
    for i = 1:size(V,1)-1
        for j = i+1:size(V,1)
            if (values(i) < 1 && values(j) > 1) || ...
                    (values(j) < 1 && values(i) > 1)
                t = (1-values(i))/(values(j)-values(i));
                points(end+1,:) = (1-t)*V(i,:)+t*V(j,:); %#ok<AGROW>
            end
        end
    end
end

function A = hullHalfspaces(vertices)
% 原点位于内部的凸包转为 A*x <= 1；不使用会扰动输入点的 Qhull QJ。
    facets = convhulln(vertices,{'Qt','Qx'});
    n = size(vertices,2);
    A = zeros(size(facets,1),n);
    keep = false(size(facets,1),1);
    for i = 1:size(facets,1)
        face = vertices(facets(i,:),:);
        if rcond(face) <= 100*eps
            % Qt 对共面的非单纯形面做三角剖分时，可能同时输出零体积
            % 单纯形。这些低维片段不是独立支撑面，完整面由其余片段覆盖。
            continue
        end
        normal = face\ones(n,1);
        if any(~isfinite(normal))
            error('terminalSet_linear:NumericalFailure','无法计算凸包支撑平面。');
        end
        A(i,:) = normal';
        keep(i) = true;
    end
    A = A(keep,:);
    if size(A,1) < n+1 || rank(A) < n
        error('terminalSet_linear:NumericalFailure','凸包缺少可靠的支撑平面。');
    end
end
