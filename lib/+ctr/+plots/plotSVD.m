function h = plotSVD(sigma,fig)
%PLOTSVD 在对数纵轴上绘制 Hankel 矩阵的奇异值。
%   h = ctr.plots.plotSVD(sigma) 新建图窗并返回曲线句柄。
%   h = ctr.plots.plotSVD(sigma,fig) 使用指定图窗。
%   sigma 为已计算的奇异值向量，本函数不执行 SVD。

    validateattributes(sigma,{'numeric'}, ...
        {'vector','nonempty','real','finite','nonnegative'},mfilename,'sigma');
    if nargin < 2 || isempty(fig)
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
end
