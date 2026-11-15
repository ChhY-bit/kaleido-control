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
    window_lenth = T-order+1;
    % 构造 Hankel 矩阵不要求满行秩；持续激励条件应另行检查输入数据。
    H = zeros(m*order,window_lenth);
    for i = 1:window_lenth
        for j = 1:order
            H((j-1)*m+1:j*m,i) = signal(:,j+i-1);
        end
    end
    if nargout > 1
        r = rank(H);
    end
    row = m*order;
end
