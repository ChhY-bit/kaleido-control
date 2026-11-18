function h = plotSVD(sigma,position,fig)
%PLOTSVD 在对数纵轴上绘制 Hankel 矩阵的奇异值并标记截断位置。
%   h = ctr.plots.plotSVD(sigma,position) 新建图窗并返回奇异值曲线句柄。
%   h = ctr.plots.plotSVD(sigma,position,fig) 使用指定图窗。
%   sigma 为已计算的奇异值向量，本函数不执行 SVD。
%   position 为最后一个保留的奇异值索引，取值为 1 到 numel(sigma)。
%   截断点以红色标记，并用红色虚线指示其横、纵坐标及索引。

    narginchk(2,3);
    validateattributes(sigma,{'numeric'}, ...
        {'vector','nonempty','real','finite','nonnegative'},mfilename,'sigma');
    validateattributes(position,{'numeric'}, ...
        {'scalar','real','finite','integer','>=',1,'<=',numel(sigma)}, ...
        mfilename,'position');
    if sigma(position) == 0
        error('plotSVD:ZeroCutoffValue', ...
            '截断位置的奇异值必须大于零，才能在对数纵轴上标记。');
    end
    if nargin < 3 || isempty(fig)
        fig = figure;
    elseif ~(isscalar(fig) && isgraphics(fig,'figure'))
        error('plotSVD:InvalidFigure','fig 必须为有效的 figure 图窗句柄。');
    end
    ax = get(fig,'CurrentAxes');
    if isempty(ax)
        ax = axes('Parent',fig);
    end
    h = semilogy(ax,sigma,'o-');
    grid(ax,'on');
    xlabel(ax,'Singular Value Index');
    ylabel(ax,'Singular Value');
    title(ax,'The SVD Result of Hankel Matrices');

    xline(ax,position,'r--',sprintf('index = %d',position), ...
        'LabelOrientation','horizontal', ...
        'LabelVerticalAlignment','bottom','HandleVisibility','off');
    yline(ax,sigma(position),'r--','HandleVisibility','off');
    line('Parent',ax,'XData',position,'YData',sigma(position), ...
        'LineStyle','none','Marker','o','MarkerSize',8, ...
        'Color','r','MarkerFaceColor','r','HandleVisibility','off');
end
