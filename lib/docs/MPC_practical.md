
## 1 问题定义

对于离散 LTI 系统

$$
\begin{aligned}
    x(k+1) &= Ax(k) + Bu(k), \\
    y(k) &= Cx(k) + Du(k),
\end{aligned}
$$

$x(k)\in\mathbb{R}^{n},u(k)\in\mathbb{R}^{m},y(k)\in\mathbb{R}^{p},k\in\mathbb{Z}_{\ge 0}.$

控制参考信号为 $r(k)\in\mathbb{R}^{p},\;R_{N}(k) := [r^\top(k),r^\top(k+1),\cdots,r^\top(k+N-1)]^\top.$

## 2 模型预测控制框架

记 $x(k+i|k)$ 为 $k$ 时刻对此后第 $i$ 个时间步的状态预测，同理有 $y(k+i|k),u(k+i|k).$

定义控制增量 $\Delta u(k+i|k) = u(k+i|k)-u(k+i-1|k)$

$$ \begin{aligned}
\min \quad & J(k) = \sum_{i=0}^{N-1} \left( \left\|r(k+i) - y(k+i|k)\right\|_Q^2 + \left\|u(k+i|k)\right\|^2_R \right) + \\
&\qquad \quad \left\| r(k+N)-y(k+N|k) \right\|_P^2 + \sum_{i=1}^{N-1}\left\|\Delta u(k+i|k)\right\|_S^2. \\
\text{s.t.} \quad & x(k|k) = x(k),\\
& x(k+i+1|k) = Ax(k+i|k) + Bu(k+i|k), \quad i = 0,1,2,\cdots,N-1\\
& y(k+i|k) = Cx(k+i|k) + Du(k+i|k), \quad i = 0,1,2,\cdots,N-1\\
& u(k+i|k) \in \mathcal{U}, \quad i = 0,1,2,\cdots,N-1\\
& y(k+i|k) \in \mathcal{Y}, \quad i = 0,1,2,\cdots,N-1\\
& \Delta u(k+i|k)\in \mathcal{D}, \quad i = 1,2,\cdots,N-1\\
& y(k+N|k) \in \mathcal{Y}_t.
\end{aligned}
$$

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
X(k+1|k) = \mathcal{S}_x x(k|k) + \mathcal{S}_u U(k)
$$

其中

$$
\begin{aligned}
X(k+1|k) &:= \begin{bmatrix}
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

### 3.2 终端情况

特别地，针对终端时刻，有 

$$
\begin{aligned}
y(k+N|k) &= CA^{N}x(k|k) + \sum_{m=0}^{N-1} CA^{N-1-m}B\, u(k+m|k) + D\, u(k+N|k) \\
&= CA^Nx(k|k) + C\Gamma_u U(k) + Du(k+N|k).
\end{aligned}
$$

其中 $$\Gamma_u = \begin{bmatrix}A^{N-1}B & A^{N-2}B & \cdots & B \end{bmatrix}.$$

并且终端控制量 $u(k+N|k)$ 一般设计为状态反馈 $$ u(k+N|k) = K\underbrace{\left( A^N x(k|k) + \Gamma_u U(k) \right)}_{x(k+N|k)}, $$ 其中 $K$ 可通过极点配置或LQR等方法获得。那么，终端输出应为

$$
y(k+N|k) = \underbrace{(C + DK)}_{C_t}(A^N x(k|k) + \Gamma_u U(k)).
$$

### 3.3 重写代价函数

定义权重增广矩阵：

$$\mathcal{Q}:=I_N\otimes Q,\quad \mathcal{R} := I_N\otimes R,\quad \mathcal{S}:=I_{N-1}\otimes S.$$

> 若考虑时变权重，则将不同时刻权重构成分块对角矩阵即可，例如 $\text{blkdiag}[Q(k),Q(k+1),\cdots,Q(k+N-1)].$

重写为标准二次规划问题：

$$ \begin{aligned}
J(k) &= \left\|Y(k|k)-R_N(k)\right\|^2_\mathcal{Q} + 
\left\| U(k)\right\|^2_\mathcal{R} + 
\left\| \Delta U(k) \right\|^2_\mathcal{S} + \| r(k+N)-y(k+N|k) \|^2_P \\
&= \left\|\mathcal{T}_u U(k)+\overbrace{\mathcal{T}_x x(k|k)-R_N(k)}^{\Phi(k)}\right\|^2_\mathcal{Q} + 
\left\| U(k)\right\|^2_\mathcal{R} + 
\left\| \mathcal{S}_\Delta U(k) \right\|^2_\mathcal{S} + \\
& \quad \; \left\| C_t\Gamma_u U(k) + \underbrace{C_tA^N x(k|k) - r(k+N)}_{\Psi(k)} \right\|^2_P \\
&= U^\top(k) \underbrace{\left( \Gamma_u^\top C_t^\top P C_t \Gamma_u + \mathcal{T}_u^\top \mathcal{Q} \mathcal{T}_u + \mathcal{R} + \mathcal{S}_\Delta^\top \mathcal{S}\mathcal{S}_\Delta \right)}_{H} U(k) + \\
& \quad \; 2\underbrace{\left( \Psi^\top(k) PC_t\Gamma_u + \Phi^\top(k) \mathcal{Q} \mathcal{T}_u \right)}_{f^\top(k)} U(k) + \underbrace{\Psi^\top(k) P \Psi(k) + \Phi^\top(k)\mathcal{Q}\Phi(k)}_{c(k)}
\end{aligned}
$$

### 3.4 重写约束条件

采用凸多面体约束（假设终端约束也已经化为多面体约束）：

$$
y_{\min}(k+i) \le y(k+i|k) \le y_{\max}(k+i), \quad i=0,1,\cdots,N-1 \\
u_{\min}(k+i) \le u(k+i|k) \le u_{\max}(k+i), \quad i=0,1,\cdots,N-1 \\
\Delta u_{\min}(k+i) \le \Delta u(k+i|k) \le \Delta u_{\max}(k+i), \quad i=1,2,\cdots,N-1\\
y_{t,\min} \le y(k+N|k) \le y_{t,\max}.\\
$$

定义

$$
Y_{lb}(k):=
\begin{bmatrix}
y_{\min}(k) \\ y_{\min}(k+1) \\ \vdots \\ y_{\min}(k+N-1)
\end{bmatrix},\quad
Y_{ub}(k):=
\begin{bmatrix}
y_{\max}(k) \\ y_{\max}(k+1) \\ \vdots \\ y_{\max}(k+N-1)
\end{bmatrix},
$$

> 若约束界限不时变，则退化为 $ Y_{lb} := \mathbf{1}_N \otimes y_{\min} $，余者同理。

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

重写为

$$
Y_{lb}(k) \le \mathcal{T}_x x(k|k) + \mathcal{T}_u U(k) \le Y_{ub}(k), \\
U_{lb}(k) \le U(k) \le U_{ub}(k), \\
\Delta_{lb} \le \mathcal{S}_\Delta U(k) \le \Delta_{ub}, \\
y_{t,\min} \le C_t\Gamma_u U(k) + C_t A^N x(k|k) \le y_{t,\max}.\\
$$

整理为矩阵：

$$
\underbrace{\begin{bmatrix}
U_{lb}(k) \\ \Delta_{lb}(k) \\ Y_{lb}(k)-\mathcal{T}_x x(k|k) \\ y_{t,\min} - C_t A^N x(k|k)
\end{bmatrix}}_{G_{lb}(k)} \le 
\underbrace{\begin{bmatrix}
I_{mN} \\ \mathcal{S}_\Delta \\ \mathcal{T}_{u} \\ C_t\Gamma_u
\end{bmatrix}}_{F} U(k) \le
\underbrace{\begin{bmatrix}
U_{ub}(k) \\ \Delta_{ub}(k) \\ Y_{ub}(k)-\mathcal{T}_x x(k|k) \\ y_{t,\max} - C_t A^N x(k|k)
\end{bmatrix}}_{G_{ub}(k)}.
$$

如果 QP 求解器要求的形式为 $F U \le G $，那么可以再通过引入负号的方式加以处理。

### 3.5 状态观测

每次采样时，需获知当前 $x(k|k) = x(k).$

若 $p=n$，$C$ 可逆且 $D=0$，则 $$x(k) = C^{-1}y(k).$$

一般情况下，采用状态观测器

$$ \begin{cases}
\hat{x}(k+1) = A \hat{x}(k) + Bu(k) + L\left(y(k)-\hat{y}(k)\right), \\ \hat{y}(k) = C\hat{x}(k) + Du(k).
\end{cases}$$


## 4 算法总结

**给定信息：**
- 系统模型 $A,B,C,D$
- 设定权重 $Q,R,S$
- 预测步数 $N$
- 终端要素 $K,y_{t,\min/\max},P$
- 控制目标 $r(k),\;k=0,1,2,\cdots$
- 约束界限 $y_{\min / \max}(k),u_{\min / \max}(k),\Delta u_{\min / \max}(k)$

**离线部分：**
- 计算 $C_t,\mathcal{Q},\mathcal{R},\mathcal{S}.$
- 构建 $\mathcal{T}_x,\mathcal{T}_u,\mathcal{S}_\Delta,\Gamma_u.$
- 构建 $H,F$

**在线部分：**
- 获取 $x(k|k)$ （或 $\hat{x}(k|k)$）
- 截取 $R_N(k),r(k+N),Y_{lb/ub}(k),U_{lb/ub}(k),\Delta_{lb/ub}(k).$
- 计算 $\Phi(k),\Psi(k),G_{lb/ub}(k).$
- 计算 $f^\top(k),c(k).$
- 求解 $\min_{U} U^\top H U + 2f^\top(k) U + c(k),\;\text{s.t.}G_{lb} \le FU \le G_{ub}.$
- 取 $U(k)$ 的第一组控制量执行。
- （可选）$U(k)$ 移位，补充终端控制量，作为热启动解