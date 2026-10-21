对于离散 LTI 系统

$$
\begin{aligned}
    x_{k+1} &= Ax_k + Bu_k, \\
    y_k &= Cx_k + Du_k,
\end{aligned}
$$

控制参考信号为 $r(k),\;k = 0,1,2,\cdots$

记 $x(k+i|k)$ 为 $k$ 时刻对此后第 $i$ 个时间步的状态预测，同理有 $y(k+i|k),u(k+i|k).$

定义控制增量 $\Delta u(k+i|k) = u(k+i|k)-u(k+i-1|k)$

$$ \begin{aligned}
\min \quad & J(k) = \sum_{i=0}^{N-1} \left( \left\|r(k+i) - y(k+i|k)\right\|_Q^2 + \left\|u(k+i|k)\right\|^2_R \right) \\
&\quad + \left\| r(k+N)-y(k+N|k) \right\|_P^2 + \sum_{i=1}^{N-1}\left\|\Delta u(k+i|k)\right\|_S^2. \\ 
\text{s.t.} \quad & x(k|k) = x_k,\\
& x(k+i+1|k) = Ax(k+i|k) + Bu(k+i|k), \quad i = 0,1,2,\cdots,N-1\\
& y(k+i|k) = Cx(k+i|k) + Du(k+i|k), \quad i = 0,1,2,\cdots,N-1\\
& u(k+i|k) \in \mathcal{U}, \quad i = 0,1,2,\cdots,N-1\\
& y(k+i|k) \in \mathcal{Y}, \quad i = 0,1,2,\cdots,N-1\\
& \Delta u(k+i|k)\in \mathcal{D}, \quad i = 1,2,\cdots,N-1\\
& y(k+N|k) \in \mathcal{Y}_f.
\end{aligned}
$$

$$\begin{aligned}
    x(k+1|k) &= Ax(k|k) + Bu(k|k), \\
    x(k+2|k) &= Ax(k+1|k) + Bu(k+1|k) \\
    &= A^2x(k|k) + ABu(k|k) + Bu(k+1|k), \\
    x(k+3|k) &= Ax(k+2|k) + Bu(k+2|k) \\
    &= A^3x(k|k) + A^2Bu(k|k) + ABu(k+1|k) + Bu(k+2|k), \\
    \vdots \\
    x(k+i|k) &= A^ix(k|k) + \sum_{m=0}^{i-1}A^{i-m-1}Bu(k+m|k), \\
    \vdots
\end{aligned}
$$

整理为 

$$
X(k|k) = \mathcal{S}_x x(k|k) + \mathcal{S}_u U(k)
$$

其中

$$
\begin{aligned}
X(k|k) &:= \begin{bmatrix}
    x(k+1|k) \\ x(k+2|k) \\ x(k+3|k) \\ \vdots \\ x(k+N|k)
\end{bmatrix}, \quad
U(k) := \begin{bmatrix}
    u(k|k) \\ u(k+1|k) \\ u(k+2|k) \\ \vdots \\ u(k+N-1|k)
\end{bmatrix}, \\[2ex]
\mathcal{S}_x &:= \begin{bmatrix}
    A \\ A^2 \\ A^3 \\ \vdots \\ A^{N}
\end{bmatrix}, \quad
\mathcal{S}_u := \begin{bmatrix}
    B & 0 & 0 & \cdots & 0 \\
    AB & B & 0 & \cdots & 0 \\
    A^2B & AB & B & \cdots & 0 \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    A^{N-1}B & A^{N-2}B & A^{N-3}B & \cdots & B
\end{bmatrix}.
\end{aligned}
$$

考虑到 $y(k+i|k) = Cx(k+i|k) + Du(k+i|k)$，因此有

$$
Y(k|k) = \mathcal{T}_x x(k|k) + \mathcal{T}_u U(k)
$$

其中

$$
\begin{aligned}
Y(k|k) &:= \begin{bmatrix}
    y(k|k) \\ y(k+1|k) \\ y(k+2|k) \\ \vdots \\ y(k+N-1|k)
\end{bmatrix}, \quad
\mathcal{T}_x := \begin{bmatrix}
    C \\ CA \\ CA^2 \\ \vdots \\ CA^{N-1}
\end{bmatrix}, \\[2ex]
\mathcal{T}_u &:= \begin{bmatrix}
    D & 0 & 0 & \cdots & 0 & 0 \\
    CB & D & 0 & \cdots & 0 & 0 \\
    CAB & CB & D & \cdots & 0 & 0 \\
    \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
    CA^{N-2}B & CA^{N-3}B & CA^{N-4}B & \cdots & CB & D
\end{bmatrix}.
\end{aligned}
$$

特别地，针对终端时刻，有 

$$
\begin{aligned}
y(k+N|k) &= CA^{N}x(k|k) \\
&\quad + \left[\begin{array}{cccc:c}
    CA^{N-1}B & CA^{N-2}B & \cdots & CB & D
\end{array}\right] \left[\begin{array}{c}
    U(k) \\ \hdashline
    u(k+N|k)
\end{array}\right].
\end{aligned}
$$

针对控制增量，定义 $$\Delta U(k):=\begin{bmatrix}
    \Delta u(k+1|k) \\ \Delta u(k+2|k) \\ \cdots \\ \Delta u(k+N-1|k)
\end{bmatrix}$$ 那么有

$$
\Delta U(k) = \underbrace{\begin{bmatrix}
    -I & I & 0 & \cdots & 0 & 0 \\
    0 & -I & I & \cdots & 0 & 0 \\
    \vdots & \vdots & \ddots & \ddots & \vdots & \vdots \\
    0 & 0 & 0 & \cdots & -I & I \\
\end{bmatrix}}_{\mathcal{S}_\Delta \in \mathbb{R}^{m(N-1)\times mN}}
U(k),\quad I\in\mathbb{R}^{m\times m}.
$$

则原代价函数可以重写为
