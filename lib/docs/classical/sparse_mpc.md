> 2026-11-07 实验记录：
> 使用 Condensed-MPC，在开环不稳定（系数1、2）、大预测时域（N=100）时，控制量、最优代价抖动，甚至无法控制。
> 推测原因是A^N叠加后数值爆炸，造成数值病态
> 改用 Sparse-MPC 后，大N值的病态结果消失
> Sparse-MPC 存在新的问题：
>   - 约束稍严就可能造成初始不可行（已验证：收紧求解容差可改善；推测是贴边运行时的过度越界问题。容差如果过紧，可能导致 Condensed 形式因数值病态问题而极大损失求解精度，若只接受“solved”结果则会报错）；
>   - 求解结果似乎容易有毛刺（待排查）。

- **Condensed MPC**
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

Rewrite the dynamics and constraints with matrices:

$$
\underbrace{\begin{bmatrix} 
-I_n & 0 & 0 & \cdots & 0 & 0 & 0 & 0 & 0 \\
A & -I_n & 0 & \cdots & 0 & B & 0 & 0 & 0 \\
0 & A & -I_n & \cdots & 0 & 0 & B & 0 & 0 \\
\vdots & \vdots & \ddots & \ddots & \vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & A & -I_n & 0 & 0 & \cdots & B
\end{bmatrix}}_{\mathcal{S}\in \mathbb{R}^{n(N+1)\times [n(N+1)+mN]}} V(k) = 
\underbrace{\begin{bmatrix}
    -x(k) \\ 0 \\ 0 \\ \vdots \\ 0
\end{bmatrix}}_{\mathcal{O}(k)\in\mathbb{R}^{n(N+1)}},
$$

$$\begin{bmatrix}
    X_{lb}(k) \\ x_{t,\min} \\ U_{lb}(k)
\end{bmatrix} \le \mathcal{D} V(k) \le \begin{bmatrix}
    X_{ub}(k) \\ x_{t,\max} \\ U_{ub}(k)
\end{bmatrix},$$

where $$ \mathcal{D} = \begin{bmatrix} I_{nN} & & \\ & E_t & \\ & & I_{mN} \end{bmatrix}.$$

Together, 

$$
\underbrace{\begin{bmatrix}
    \mathcal{O}(k) \\ X_{lb}(k) \\ x_{t,\min} \\ U_{lb}(k)
\end{bmatrix}}_{
    \mathcal{G}_{lb}(k)} \le \underbrace{\begin{bmatrix}
\mathcal{S} \\ \mathcal{D}
\end{bmatrix}}_{\mathcal{F}} V(k) \le \underbrace{\begin{bmatrix}
    \mathcal{O}(k) \\ X_{ub}(k) \\ x_{t,\max} \\ U_{ub}(k)
\end{bmatrix}}_{\mathcal{G}_{ub}(k)}.
$$

Rewrite the cost function:

$$
J(k) = \left\| V(k) \right\|_{\mathcal{W}}^2,
$$

where 

$$
\mathcal{W}:=\begin{bmatrix}
I_{N} \otimes Q & & \\
& P & \\
& & I_N \otimes R \\
\end{bmatrix}.
$$

In conclusion, the Sparse-MPC problem is

$$
    \min_V \; V^\top \mathcal{W} V, \quad \text{s.t.} \; \mathcal{G}_{lb} \le \mathcal{F}V \le \mathcal{G}_{ub}.
$$