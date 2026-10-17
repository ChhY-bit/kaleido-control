classdef mpc < handle
%MPC 标准线性模型预测控制

    properties(SetAccess=private)
        control_lb
        control_ub
        state_lb
        state_ub
        output_lb
        output_ub

        Q
        R
        N
    end

    
end