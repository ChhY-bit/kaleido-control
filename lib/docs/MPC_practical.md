
## 1 问题定义

对于离散 LTI 系统

$$ x(k+1) = Ax(k) + Bu(k),$$

$x(k)\in\mathbb{R}^{n},u(k)\in\mathbb{R}^{m},y(k)\in\mathbb{R}^{p},k\in\mathbb{Z}_{\ge 0}.$

同时给定期望轨迹序列 $x_r(k),u_r(k).$

## 2 线性二次模型预测控制

记 $x(k+i|k)$ 为 $k$ 时刻对此后第 $i$ 个时间步的状态预测，同理有 $u(k+i|k).$

定义控制增量 $\Delta u(k+i|k) = u(k+i|k)-u(k+i-1|k)$

$$ \begin{aligned}
\min_{u(\cdot|k)} \quad & J(k) = \sum_{i=1}^{N} \left\|x(k+i|k)-x_r(k+i)\right\|_Q^2 + \sum_{i=0}^{N-1} \left\|u(k+i|k)-u_r(k+i)\right\|^2_R + \\
&\qquad \quad \left\| x(k+N|k) - x_r(k+N) \right\|_P^2 + \sum_{i=1}^{N-1}\left\|\Delta u(k+i|k)\right\|_S^2. \\
\text{s.t.} \quad & x(k|k) = x(k),\\
& x(k+i+1|k) = Ax(k+i|k) + Bu(k+i|k), \quad i = 0,1,2,\cdots,N-1\\
& u(k+i|k) \in \mathcal{U}, \quad i = 0,1,2,\cdots,N-1\\
& x(k+i|k) \in \mathcal{X}, \quad i = 1,2,3,\cdots,N\\
& \Delta u(k+i|k)\in \mathcal{D}, \quad i = 1,2,\cdots,N-1\\
& x(k+N|k)-x_r(k+N) \in \mathcal{E}_t.
\end{aligned}
$$F

## 3 优化求解

### 3.1 递推预测

$$\begin{aligned}
    x(k+1|k) &= Ax(k|k) + Bu(k|k), \\
    x(k+2|k) &= Ax(k+1|k) + Bu(k+1|k) \\
    &= A^2x(k|k) + ABu(k|k) + Bu(k+1|k), \\
    x(k+3|k) &= Ax(k+2|k) + Bu(k+2|k) \\
    &= A^3x(k|k) + A^2Bu(k|k) + ABu(k+1|k) + Bu(k+2|k), \\
    \vdots \\
    x(k+N|k) &= A^Nx(k|k) + \sum_{m=0}^{N-1}A^{N-m-1}Bu(k+m|k).
\end{aligned}
$$

整理为 

$$
X(k) = \mathcal{T}_x x(k|k) + \mathcal{T}_u U(k)
$$

其中

$$
\begin{aligned}
X(k) &:= \begin{bmatrix}
    x(k+1|k) \\ x(k+2|k) \\ x(k+3|k) \\ \vdots \\ x(k+N|k)
\end{bmatrix}, \quad
U(k) := \begin{bmatrix}
    u(k|k) \\ u(k+1|k) \\ u(k+2|k) \\ \vdots \\ u(k+N-1|k)
\end{bmatrix}, \\[2ex]
\mathcal{T}_x &:= \begin{bmatrix}
    A \\ A^2 \\ A^3 \\ \vdots \\ A^{N}
\end{bmatrix}, \quad
\mathcal{T}_u := \begin{bmatrix}
    B & 0 & 0 & \cdots & 0 \\
    AB & B & 0 & \cdots & 0 \\
    A^2B & AB & B & \cdots & 0 \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    A^{N-1}B & A^{N-2}B & A^{N-3}B & \cdots & B
\end{bmatrix}.
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
\end{bmatrix}}_{\mathcal{T}_\Delta \in \mathbb{R}^{m(N-1)\times mN}}
U(k),\quad I\in\mathbb{R}^{m\times m}.
$$

特别地，针对终端时刻，有 

$$
x(k+N|k) = A^Nx(k|k) + \Gamma_u U(k).
$$

其中 $$\Gamma_u = \begin{bmatrix} O_{n\times n(N-1)} & I_n \end{bmatrix} \mathcal{T}_u = \begin{bmatrix}A^{N-1}B & A^{N-2}B & \cdots & B \end{bmatrix}.$$

### 3.3 重写代价函数

定义权重增广矩阵：

$$\mathcal{Q}:=I_N\otimes Q,\quad \mathcal{R} := I_N\otimes R,\quad \mathcal{S}:=I_{N-1}\otimes S.$$

> 若考虑时变权重，则将不同时刻权重构成分块对角矩阵即可，例如 $\text{blkdiag}[Q(k+1),Q(k+2),\cdots,Q(k+N)].$

取参考向量：

$$
X_r(k) := \begin{bmatrix} x_r(k+1) \\ x_r(k+2) \\ \vdots \\ x_r(k+N) \end{bmatrix},\quad
U_r(k) := \begin{bmatrix} u_r(k) \\ u_r(k+1) \\ \vdots \\ u_r(k+N-1) \end{bmatrix}.
$$

重写为标准二次规划问题：

$$ \begin{aligned}
J(k) &= \| x(k+N|k) - x_r(k+N) \|^2_P + \left\|X(k)-X_r(k)\right\|^2_\mathcal{Q} + 
\left\| U(k) - U_r(k)\right\|^2_\mathcal{R} + 
\left\| \Delta U(k) \right\|^2_\mathcal{S} \\
&= \left\| \Gamma_u U(k) + \underbrace{A^N x(k|k) - x_r(k+N)}_{\Psi(k)} \right\|^2_P + \left\|\mathcal{T}_u U(k)+\underbrace{\mathcal{T}_x x(k|k)-X_r(k)}_{\Phi(k)}\right\|^2_\mathcal{Q} + \\
& \quad \; \left\| U(k) - U_r(k)\right\|^2_\mathcal{R} + 
\left\| \mathcal{T}_\Delta U(k) \right\|^2_\mathcal{S} \\
& \quad \; \\
&= U^\top(k) \underbrace{\left( \Gamma_u^\top P \Gamma_u + \mathcal{T}_u^\top \mathcal{Q} \mathcal{T}_u + \mathcal{R} + \mathcal{T}_\Delta^\top \mathcal{S}\mathcal{T}_\Delta \right)}_{H} U(k) + \\
& \quad \; 2\underbrace{\left( \Psi^\top(k) P\Gamma_u + \Phi^\top(k) \mathcal{Q} \mathcal{T}_u - U_r^\top(k) \mathcal{R} \right)}_{f^\top(k)} U(k) + \\ 
& \quad \; \underbrace{\Psi^\top(k) P \Psi(k) + \Phi^\top(k)\mathcal{Q}\Phi(k) + U^\top_r(k)\mathcal{R} U_r(k)}_{c(k)}
\end{aligned}
$$

### 3.4 重写约束条件

采用箱形约束（假设终端约束也已经化为箱形约束）：

$$
x_{\min}(k+i) \le x(k+i|k) \le x_{\max}(k+i), \quad i=1,2,\cdots,N \\
u_{\min}(k+i) \le u(k+i|k) \le u_{\max}(k+i), \quad i=0,1,\cdots,N-1 \\
\Delta u_{\min}(k+i) \le \Delta u(k+i|k) \le \Delta u_{\max}(k+i), \quad i=1,2,\cdots,N-1\\
e_{t,\min} \le x(k+N|k) - x_r(k+N) \le e_{t,\max}.\\
$$

定义

$$
X_{lb}(k):=
\begin{bmatrix}
x_{\min}(k+1) \\ x_{\min}(k+2) \\ \vdots \\ x_{\min}(k+N)
\end{bmatrix},\quad
X_{ub}(k):=
\begin{bmatrix}
x_{\max}(k+1) \\ x_{\max}(k+2) \\ \vdots \\ x_{\max}(k+N)
\end{bmatrix},
$$

$$
\Delta_{lb}(k):=
\begin{bmatrix}
\Delta u_{\min}(k+1) \\ \Delta u_{\min}(k+2) \\ \vdots \\ \Delta u_{\min}(k+N-1)
\end{bmatrix},\quad
\Delta_{ub}(k):=
\begin{bmatrix}
\Delta u_{\max}(k+1) \\ \Delta u_{\max}(k+2) \\ \vdots \\ \Delta u_{\max}(k+N-1)
\end{bmatrix},
$$

$$
U_{lb}(k):=
\begin{bmatrix}
u_{\min}(k) \\ u_{\min}(k+1) \\ \vdots \\ u_{\min}(k+N-1)
\end{bmatrix},\quad
U_{ub}(k):=
\begin{bmatrix}
u_{\max}(k) \\ u_{\max}(k+1) \\ \vdots \\ u_{\max}(k+N-1)
\end{bmatrix}.
$$

> 若约束界限不时变，则退化为 $ X_{lb} := \mathbf{1}_N \otimes x_{\min} $，余者同理。

重写为

$$
X_{lb}(k) \le \mathcal{T}_x x(k|k) + \mathcal{T}_u U(k) \le X_{ub}(k), \\
U_{lb}(k) \le U(k) \le U_{ub}(k), \\
\Delta_{lb} \le \mathcal{T}_\Delta U(k) \le \Delta_{ub}, \\
e_{t,\min} \le \Gamma_u U(k) + A^N x(k|k) - x_r(k+N) \le e_{t,\max}.\\
$$

整理为矩阵：

$$
\underbrace{\begin{bmatrix}
U_{lb}(k) \\ \Delta_{lb}(k) \\ X_{lb}(k)-\mathcal{T}_x x(k|k) \\ e_{t,\min} - \Psi(k)
\end{bmatrix}}_{G_{lb}(k)} \le 
\underbrace{\begin{bmatrix}
I_{mN} \\ \mathcal{T}_\Delta \\ \mathcal{T}_{u} \\ \Gamma_u
\end{bmatrix}}_{F} U(k) \le
\underbrace{\begin{bmatrix}
U_{ub}(k) \\ \Delta_{ub}(k) \\ X_{ub}(k)-\mathcal{T}_x x(k|k) \\ e_{t,\max} - \Psi(k)
\end{bmatrix}}_{G_{ub}(k)}.
$$

如果 QP 求解器要求的形式为 $F U \le G $，那么可以再通过引入负号的方式加以处理。

### 3.5 状态观测

每次采样时，需获知当前 $x(k|k) = x(k).$

如果 $x(k)$ 不能直接获知，但系统 $(A,C)$ 可检测，则可采用状态观测器

$$ \begin{cases}
\hat{x}(k+1) = A \hat{x}(k) + Bu(k) + L\left(y(k)-\hat{y}(k)\right), \\ \hat{y}(k) = C\hat{x}(k) + Du(k).
\end{cases}$$

## 4 算法总结

**给定信息：**
- 系统模型 $A,B,C,D$
- 设定权重 $Q,R,S$
- 预测步数 $N$
- 终端要素 $e_{t,\min/\max},P$
- 约束界限 $x_{\min / \max}(k),u_{\min / \max}(k),\Delta u_{\min / \max}(k)$

**离线部分：**
- 构建 $\mathcal{Q},\mathcal{R},\mathcal{S}.$
- 构建 $\mathcal{T}_x,\mathcal{T}_u,\mathcal{T}_\Delta,\Gamma_u.$
- 构建 $H,F$

**在线部分：**
- 获取 $x(k|k)$ （或 $\hat{x}(k|k)$）
- 截取 $X_{lb/ub}(k),U_{lb/ub}(k),\Delta_{lb/ub}(k).$
- 截取 $X_r(k),U_r(k),x_r(k+N).$
- 计算 $\Phi(k),\Psi(k),G_{lb/ub}(k).$
- 计算 $f^\top(k),c(k).$
- 求解 $\min_{U} U^\top H U + 2f^\top(k) U + c(k),\;\text{s.t.}G_{lb} \le FU \le G_{ub}.$
- 取 $U(k)$ 的第一组控制量执行。
- （可选）$U(k)$ 移位，补充终端控制量，作为可行解/热启动解