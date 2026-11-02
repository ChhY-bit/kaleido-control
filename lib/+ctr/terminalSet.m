function [xt_min,xt_max,Et] = terminalSet(At,M,alpha)
    % At    -   Schur 稳定的非奇异闭环系统矩阵
    % M     -   终端集形状矩阵（Lyapunov解）
    % alpha -   终端集大小因子
    n = size(At,1);
    r = sqrt(alpha*diag(inv(M)));
    G = M;
    k = 0;
    Et = eye(n);
    while r'*abs(G)*r > alpha
        G = At'*G*At;
        k = k+1;
        Et = [Et;(At\eye(n))^k]; %#ok<AGROW>
        xt_min = kron(ones(k+1,1),-r);
        xt_max = kron(ones(k+1,1),r);
    end
    % xt_min = kron(ones(k+1,1),-r);
    % xt_max = kron(ones(k+1,1),r);
end