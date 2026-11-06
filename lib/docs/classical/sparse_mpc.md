- **Dense MPC**
    - **Advantages**
        - Simplicity of model expression
        - Easy to send into QP solver interface
        - High efficiency of solving due to less dimension
    - **Shortages**
        - Numerical pathological problem
        - Sensitivity to predictive horizon with open-loop unstable system
    - **Ralated reference:** <u> Horizon-Independent Preconditioner Design for Linear Predictive Control (IEEE Transactions on Automatic Control,2023)</u>

- **Sparse MPC**
    - Transist single-variable into multi-variable decision.
    - Introduce prediction by equality constraints of dynamics.

Consider the classic QP-MPC optimization problem:

$$
\begin{aligned}
    \min_{u(\cdot|k)}\quad & J(k) = \left\| x(k+N|k) \right\|_P^2 + \sum_{i=0}^{N-1} \left(\|x(k+i|k)\|_Q^2 + \left\|u(k+i|k)\right\|_R^2\right). \\
    \text{s.t.}\quad & x(k|k) = x(k), \\
    & x(k+i+1|k) = Ax(k+i|k) + Bu(k+i|k), \quad i=0,1,\cdots N-1\\
    & x(k+i|k) \in\mathcal{X}, \quad i=0,1,\cdots N-1\\
    & u(k+i|k)\in\mathcal{U}, \quad i=0,1,\cdots N-1\\
    & x(k+N|k)\in\mathcal{X}_t.
\end{aligned}
$$

Take the immediate variable $x$ as decision variable. In other words, define a new vector $$V(k):=\left\{\begin{bmatrix} x(k|k) \\ x(k+1|k) \\ \vdots \\ x(k+N|k) \end{bmatrix}^\top, \begin{bmatrix} u(k|k) \\ u(k+1|k) \\ \vdots \\ u(k+N-1|k) \end{bmatrix}^\top\right\}^\top. $$

Then the recursive equations can be obtained according to the dynamics:

$$
\begin{matrix}
x_k & & & & & & & & & = & x(k) \\
Ax_{k} & -x_{k+1} & & & & +Bu_k & & & & = & 0 \\
& Ax_{k+1} & -x_{k+2} & & & & +Bu_{k+1} & & & = & 0 \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
& & & Ax_{k+N-1} & -x_{k+N} & & & & +Bu_{k+N-1} & = & 0 \\
\end{matrix}
$$

Rewrite the dynamics with matrices:

$$
\underbrace{\begin{bmatrix} 
I & 0 & 0 & \cdots & 0 & 0 & 0 & 0 & 0 \\
A & -I & 0 & \cdots & 0 & B & 0 & 0 & 0 \\
0 & A & -I & \cdots & 0 & 0 & B & 0 & 0 \\
\vdots & \vdots & \ddots & \ddots & \vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & A & -I & 0 & 0 & \cdots & B
\end{bmatrix}}_{\mathcal{S}\in \mathbb{R}^{n(N+1)\times [n(N+1)+mN]}} V(k) = 
\underbrace{\begin{bmatrix}
    x(k) \\ 0 \\ 0 \\ \vdots \\ 0
\end{bmatrix}}_{\mathcal{O}(k)}
$$