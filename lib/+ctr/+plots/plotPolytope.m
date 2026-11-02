function [h,V] = plotPolytope(A,lb,ub,options,fig)
%PLOTPOLYTOPE 绘制二维或三维区域 lb <= A*x <= ub。
%   h = ctr.plots.plotPolytope(A,lb,ub)
%   h = ctr.plots.plotPolytope(A,lb,ub,options,fig)
%   [h,V] = ctr.plots.plotPolytope(...) 同时返回顶点，每行对应一个顶点。
%
%   A       - 有限实数矩阵，列数只能为 2 或 3。
%   lb,ub   - 每个约束的下、上界；可为行/列向量或标量。
%             可使用 -Inf 和 Inf 表示没有下界和上界。
%   options - 标量结构体，省略或 [] 使用默认值：
%       LineStyle = '-'            轮廓线型
%       LineWidth = 1              轮廓线宽
%       EdgeColor = [0 0 0]        轮廓颜色
%       EdgeAlpha = 1              轮廓透明度
%       FaceColor = [0 .447 .741]  填充颜色（'none' 表示不填充）
%       FaceAlpha = .25            填充透明度
%       Tolerance = 1e-9           几何计算相对容差
%       还可指定其他 patch 外观属性，例如 DisplayName、Marker。
%       不允许覆盖 Parent、Faces、Vertices、XData、YData、ZData。
%   fig     - figure 图窗句柄；省略或 [] 使用当前图窗。
%
%   h 为 patch 句柄。保留目标坐标区中的原有图形及 hold 状态。
%   空集返回 h = gobjects(0)、V = zeros(0,size(A,2))，不改变图窗。
%   无界区域报错 plotPolytope:UnboundedRegion，不进行任意裁剪。
%   有界退化区域仍可绘制；孤立点默认用 '.' 标记显示。
%   不依赖 Optimization Toolbox 或第三方多面体工具箱。
%   采用约束交点枚举法，适合约束数量适中的二维、三维问题。
%
%   示例（先将 lib 文件夹加入 MATLAB 路径）：
%       fig = figure(1);
%       options = struct('LineStyle','--','EdgeColor','r', ...
%                        'FaceColor',[1 .6 .2],'FaceAlpha',.3);
%       ctr.plots.plotPolytope(eye(2),[-1;-2],[1;2],options,fig);
%       ctr.plots.plotPolytope(eye(3),-1,1,[],figure(2));

    narginchk(3,5);
    validateattributes(A,{'numeric'},{'real','finite','2d'},mfilename,'A');
    [m,n] = size(A);
    if ~ismember(n,[2 3])
        error('plotPolytope:InvalidDimension','A 的列数必须为 2 或 3。');
    end
    A = full(double(A));
    lb = expandBound(lb,m,'lb');
    ub = expandBound(ub,m,'ub');
    if nargin < 4 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('plotPolytope:InvalidOptions','options 必须为标量结构体或 []。');
    end
    if nargin < 5
        fig = [];
    end
    if ~isempty(fig) && ~(isscalar(fig) && isgraphics(fig,'figure'))
        error('plotPolytope:InvalidFigure','fig 必须为有效的 figure 图窗句柄。');
    end

    style = struct('LineStyle','-','LineWidth',1,'EdgeColor',[0 0 0], ...
                   'EdgeAlpha',1,'FaceColor',[0 .447 .741],'FaceAlpha',.25);
    tol = 1e-9;
    names = fieldnames(options);
    reserved = {'Parent','Faces','Vertices','XData','YData','ZData'};
    for i = 1:numel(names)
        name = names{i};
        if strcmpi(name,'Tolerance')
            tol = options.(name);
        elseif any(strcmpi(name,reserved))
            error('plotPolytope:ReservedOption','options 不允许覆盖 %s。',name);
        else
            % 属性名不区分大小写；覆盖默认值，避免重复传入同一属性。
            existing = fieldnames(style);
            match = find(strcmpi(name,existing),1);
            if ~isempty(match)
                name = existing{match};
            end
            style.(name) = options.(names{i});
        end
    end
    validateattributes(tol,{'numeric'}, ...
        {'real','finite','scalar','positive','<',1},mfilename,'Tolerance');

    h = gobjects(0);
    V = zeros(0,n);
    if any(lb > ub) || any(lb == Inf) || any(ub == -Inf)
        return
    end

    % 将双侧约束转为 B*x <= b，去除无限界和恒真零行。
    upper = isfinite(ub);
    lower = isfinite(lb);
    B = [A(upper,:);-A(lower,:)];
    b = [ub(upper);-lb(lower)];
    lengths = vecnorm(B,2,2);
    if any(b(lengths == 0) < 0)
        return
    end
    nonzero = lengths > 0;
    B = B(nonzero,:)./lengths(nonzero);
    b = b(nonzero)./lengths(nonzero);

    % 每个顶点至少具有 n 个线性无关的活动约束。
    if size(B,1) >= n
        sets = nchoosek(1:size(B,1),n);
        for i = 1:size(sets,1)
            active = B(sets(i,:),:);
            if rcond(active) <= tol
                continue
            end
            point = active\b(sets(i,:));
            if isFeasible(B,b,point,tol)
                V(end+1,:) = point'; %#ok<AGROW>
            end
        end
    end

    % 无顶点不代表空集：半空间、平面、直线等仍可能可行。
    % 非空闭凸集的最近原点投影具有至多 n 个独立活动约束。
    feasible = ~isempty(V) || isFeasible(B,b,zeros(n,1),tol);
    if ~feasible
        for k = 1:min(n-1,size(B,1))
            sets = nchoosek(1:size(B,1),k);
            for i = 1:size(sets,1)
                active = B(sets(i,:),:);
                if rank(active,tol) < k
                    continue
                end
                point = lsqminnorm(active,b(sets(i,:)),tol);
                if isFeasible(B,b,point,tol)
                    feasible = true;
                    break
                end
            end
            if feasible
                break
            end
        end
    end
    if ~feasible
        return
    end
    if isUnbounded(B,n,tol)
        error('plotPolytope:UnboundedRegion', ...
              '区域无界，无法绘制完整轮廓；请增加有限的边界约束。');
    end
    if isempty(V)
        error('plotPolytope:NumericalFailure', ...
              '未能可靠计算有界区域的顶点，请检查约束条件或调整 Tolerance。');
    end
    V = uniquetol(V,tol,'ByRows',true,'DataScale',1);

    centered = V-mean(V,1);
    [~,S,basis] = svd(centered,0);
    dimension = sum(diag(S) > tol*max(1,norm(centered,'fro')));
    if dimension == 0
        faces = 1;
        if ~any(strcmpi(fieldnames(style),'Marker'))
            style.Marker = '.';
        end
    elseif dimension == 1
        [~,order] = sort(centered*basis(:,1));
        faces = [order(1),order(end),order(1)];
    elseif dimension == 2
        faces = orderFace(V,1:size(V,1));
    else
        % 按原始支撑平面构造多边形面，避免凸包三角剖分的对角线。
        faceList = {};
        for i = 1:size(B,1)
            residual = V*B(i,:)'-b(i);
            indices = find(abs(residual) <= tol*(1+abs(b(i))+abs(V)*abs(B(i,:)')));
            if numel(indices) < 3 || ...
                    rank(V(indices,:)-V(indices(1),:),tol) < 2
                continue
            end
            duplicate = cellfun(@(f) isequal(sort(f(:)),sort(indices)),faceList);
            if ~any(duplicate)
                face = orderFace(V,indices);
                normal = cross(V(face(2),:)-V(face(1),:), ...
                               V(face(3),:)-V(face(1),:));
                if normal*B(i,:)' < 0
                    face = fliplr(face);
                end
                faceList{end+1} = face; %#ok<AGROW>
            end
        end
        faces = nan(numel(faceList),max(cellfun(@numel,faceList)));
        for i = 1:numel(faceList)
            faces(i,1:numel(faceList{i})) = faceList{i};
        end
    end

    if isempty(fig)
        fig = gcf;
    end
    ax = get(fig,'CurrentAxes');
    if isempty(ax)
        ax = axes('Parent',fig);
    end
    style.Parent = ax;
    style.Faces = faces;
    style.Vertices = V;
    h = patch(style);
    axis(ax,'equal');
    if n == 3
        view(ax,3);
    end
end

function bound = expandBound(bound,m,name)
    validateattributes(bound,{'numeric'}, ...
        {'real','nonnan'},mfilename,name);
    if ~(isvector(bound) || isempty(bound)) || ...
            ~(isscalar(bound) || numel(bound) == m)
        error('plotPolytope:InvalidBounds', ...
              '%s 必须为标量或含有 size(A,1) 个元素的向量。',name);
    end
    bound = full(double(bound(:)));
    if isscalar(bound)
        bound = repmat(bound,m,1);
    end
end

function feasible = isFeasible(B,b,point,tol)
    feasible = all(isfinite(point)) && ...
        all(B*point <= b+tol*(1+abs(b)+abs(B)*abs(point)));
end

function unbounded = isUnbounded(B,n,tol)
    % 非平凡衰退方向 d 满足 B*d <= 0。
    if rank(B,tol) < n
        unbounded = true;
        return
    end
    sets = nchoosek(1:size(B,1),n-1);
    unbounded = false;
    for i = 1:size(sets,1)
        active = B(sets(i,:),:);
        if n == 2
            direction = [-active(2);active(1)];
        else
            direction = cross(active(1,:),active(2,:))';
        end
        length = norm(direction);
        if length <= tol
            continue
        end
        direction = direction/length;
        if all(B*direction <= tol) || all(-B*direction <= tol)
            unbounded = true;
            return
        end
    end
end

function face = orderFace(V,indices)
    points = V(indices,:)-mean(V(indices,:),1);
    [~,~,basis] = svd(points,0);
    plane = points*basis(:,1:2);
    order = convhull(plane(:,1),plane(:,2));
    face = reshape(indices(order(1:end-1)),1,[]);
end

