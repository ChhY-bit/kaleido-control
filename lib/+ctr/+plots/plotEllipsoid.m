function [h,V] = plotEllipsoid(P,alpha,options,fig)
%PLOTELLIPSOID 绘制二维椭圆或三维椭球区域 x'*P*x <= alpha。
%   h = ctr.plots.plotEllipsoid(P,alpha)
%   h = ctr.plots.plotEllipsoid(P,alpha,options,fig)
%   [h,V] = ctr.plots.plotEllipsoid(...) 返回边界采样点，每行对应一个点。
%
%   P       - 2x2 或 3x3 有限实数对称正定形状矩阵。
%             半轴长度为 sqrt(alpha/lambda_i)，lambda_i 为 P 的特征值。
%   alpha   - 有限实数标量。alpha > 0 绘制区域，alpha = 0 绘制原点，
%             alpha < 0 返回空句柄和空采样点矩阵，不改变图窗。
%   options - 标量结构体，省略或 [] 使用默认值（同 plotPolytope）：
%       LineStyle = '-'            轮廓/网格线型
%       LineWidth = 1              轮廓/网格线宽
%       EdgeColor = [0 0 0]        轮廓/网格颜色（'none' 表示不显示）
%       EdgeAlpha = 1              轮廓/网格透明度
%       FaceColor = [0 .447 .741]  填充颜色（'none' 表示不填充）
%       FaceAlpha = .25            填充透明度
%       Resolution = 128 (二维) 或 32 (三维)
%                                  正整数，至少为 4；三维为球面分段数
%       Tolerance = 1e-9           检查 P 对称性的相对容差
%       还可指定其他 patch 外观属性，例如 DisplayName、Marker。
%       不允许覆盖 Parent、Faces、Vertices、XData、YData、ZData。
%   fig     - figure 图窗句柄；省略或 [] 使用当前图窗。
%
%   h 为 patch 句柄。二维用闭合多边形近似椭圆，三维用面片近似椭球。
%   保留目标坐标区中的原有图形及 hold 状态，并设置等比例坐标轴。
%   三维 EdgeColor 控制曲面网格颜色；设置为 'none' 可显示光滑填充。
%   alpha = 0 时，原点默认用 '.' 标记显示。
%   不依赖 Optimization Toolbox 或第三方几何工具箱。
%
%   示例（先将 lib 文件夹加入 MATLAB 路径）：
%       options = struct('LineStyle','--','EdgeColor','r', ...
%                        'FaceColor',[1 .6 .2],'FaceAlpha',.3);
%       ctr.plots.plotEllipsoid([2 1;1 3],1,options,figure(1));
%       ctr.plots.plotEllipsoid(diag([1 2 4]),2,[],figure(2));

    narginchk(2,4);
    validateattributes(P,{'numeric'},{'real','finite','2d'},mfilename,'P');
    n = size(P,1);
    if ~ismember(n,[2 3]) || size(P,2) ~= n
        error('plotEllipsoid:InvalidDimension','P 必须为 2x2 或 3x3 方阵。');
    end
    P = full(double(P));
    validateattributes(alpha,{'numeric'}, ...
        {'real','finite','scalar'},mfilename,'alpha');
    alpha = double(alpha);
    if nargin < 3 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('plotEllipsoid:InvalidOptions','options 必须为标量结构体或 []。');
    end
    if nargin < 4
        fig = [];
    end
    if ~isempty(fig) && ~(isscalar(fig) && isgraphics(fig,'figure'))
        error('plotEllipsoid:InvalidFigure','fig 必须为有效的 figure 图窗句柄。');
    end

    style = struct('LineStyle','-','LineWidth',1,'EdgeColor',[0 0 0], ...
                   'EdgeAlpha',1,'FaceColor',[0 .447 .741],'FaceAlpha',.25);
    tol = 1e-9;
    if n == 2
        resolution = 128;
    else
        resolution = 32;
    end
    names = fieldnames(options);
    reserved = {'Parent','Faces','Vertices','XData','YData','ZData'};
    for i = 1:numel(names)
        name = names{i};
        if strcmpi(name,'Tolerance')
            tol = options.(name);
        elseif strcmpi(name,'Resolution')
            resolution = options.(name);
        elseif any(strcmpi(name,reserved))
            error('plotEllipsoid:ReservedOption','options 不允许覆盖 %s。',name);
        else
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
    validateattributes(resolution,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',4},mfilename,'Resolution');
    resolution = double(resolution);

    if norm(P-P','fro') > tol*norm(P,'fro')
        error('plotEllipsoid:NotSymmetric','P 必须为对称矩阵。');
    end
    % 消除容差范围内的非对称舍入误差；二次型仅取决于对称部分。
    P = P/2+P'/2;
    [R,flag] = chol(P);
    if flag ~= 0
        error('plotEllipsoid:NotPositiveDefinite', ...
              'P 必须为正定矩阵，才能确定有界椭圆或椭球。');
    end

    h = gobjects(0);
    V = zeros(0,n);
    if alpha < 0
        return
    elseif alpha == 0
        V = zeros(1,n);
        faces = 1;
        if ~any(strcmpi(fieldnames(style),'Marker'))
            style.Marker = '.';
        end
    else
        if n == 2
            theta = (0:resolution-1)'*(2*pi/resolution);
            unit = [cos(theta),sin(theta)];
            faces = 1:resolution;
        else
            % 复用 sphere 的经纬参数；四边形面避免添加三角剖分对角线。
            [X,Y,Z] = sphere(resolution);
            unit = [X(:),Y(:),Z(:)];
            indices = reshape(1:numel(X),size(X));
            first = indices(1:end-1,1:end-1);
            second = indices(1:end-1,2:end);
            third = indices(2:end,2:end);
            fourth = indices(2:end,1:end-1);
            faces = [first(:),second(:),third(:),fourth(:)];
        end
        % P = R'*R，将单位圆/球的边界通过 R\ 映射为目标二次型边界。
        V = (R\(sqrt(alpha)*unit'))';
        if any(~isfinite(V(:)))
            error('plotEllipsoid:NumericalFailure', ...
                  '边界坐标超出浮点数范围，请调整 P 或 alpha 的尺度。');
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
