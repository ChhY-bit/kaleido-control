## 1 可行性问题

### 1.1 可行解的构造

优化问题：

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
$$

其中 $\Delta u(k+i|k) := u(k+i|k)-u(k+i-1|k).$

假设已有满足各类约束的最优解：

$$ U^*(k) = \begin{bmatrix}
    u^*(k|k) \\ u^*(k+1|k) \\ \vdots \\ u^*(K+N-1|k)
\end{bmatrix}, $$

尽管只有 $u^*(k|k)$ 被实际执行，但剩余项仍可用于构造下一时刻的可行解，即 $$\bar{u}(k+i|k+1) = \begin{cases} u^*(k+i|k) , & i=1,2,\cdots,N-1 \\ \kappa(\bar{x}(k+N|k+1)), & i = N \end{cases} $$ 其中 $\kappa(\cdot)$ 为终端控制器。

### 1.2 终端控制器、终端不变集的设计

#### 1.2.1 给定参考为系统轨迹的情形

如果给定参考 $(x_r,u_r)$ 是一条可行的系统轨迹，即满足 $$x_r(k+1) = A x_r(k) + Bu_r(k),$$ 那么可以定义 $x_e:=x-x_r,u_e:=u-u_r$ 并得到误差动态 $$x_e(k+1) = A x_e(k) + Bu_e(k),$$

进而可以证明，取终端控制器为 $$\kappa(\bar{x}(k+N|k+1)) = K_t\left[\bar{x}(k+N|k+1)-x_r(k+N)\right] + u_r(k+N),$$ 并且 $K_t$ 能使得闭环系统 $A_t:=A+BK_t$ 是 Schur 稳定的，即 Lyapunov 方程 $$A_t^\top M A_t - M = - N,\quad N^\top = N \succ 0$$ 有对称正定解 $M$ 时，以下集合为不变集： $$\mathcal{E}_t = \left\{ e\in \mathbb{R}^n | e^\top M e \le \alpha ,\;\alpha \ge 0 \right\}.$$ 

$K_t$ 一般由 LQR 或极点配置等方法获得。基于此，可行解为

$$\bar{U}(k+1) = \begin{bmatrix}
    u^*(k+1|k) \\ u^*(k+2|k) \\ \vdots \\ u^*(k+N-1|k) \\ K_t \bar{x}_e(k+N|k+1)
\end{bmatrix}. $$

在系统标称的条件下执行可行解，有 $\bar{x}_e(k+N|k+1) = x_e^*(k+N|k)$，因此只要选取合适的终端集半径 $\alpha$，使得 $$x_e^*(k+N|k) \in \mathcal{E}_t \subseteq \mathcal{X} - x_r(k+N), $$ 并且 $$K_t x_e^*(k+N|k) \in \mathcal{U},$$ 就能保证 $k+1$ 时刻的优化问题依然可行。

#### 1.2.2 给定参考为任意信号的情形

暂略，需要考虑鲁棒不变集。

## 2 稳定性问题

考虑 $k+1$ 时刻的执行可行解后的代价函数

$$ \begin{aligned}
\bar{J}(k+1) &= \sum_{i=1}^{N} \left\|\bar{x}(k+1+i|k+1)-x_r(k+1+i)\right\|_Q^2 + \\
& \quad \sum_{i=0}^{N-1} \left\|\bar{u}(k+1+i|k+1)-u_r(k+1+i)\right\|^2_R + \\
& \quad \left\| \bar{x}(k+1+N|k+1) - x_r(k+1+N) \right\|_P^2 + \\ 
& \quad \sum_{i=1}^{N-1}\left\|\Delta \bar{u}(k+1+i|k+1)\right\|_S^2.
\end{aligned}
$$

考虑到标称条件下，$$\bar{x}(k+i|k+1) = x^*(k+i|k),\quad i=1,2,\cdots,N$$ 那么作差可得

$$
\bar{J}(k+1) - J^*(k) = \left\|\bar{x}(k+1+N|k+1)-x_r(k+1+N)\right\|_Q^2 - \left\|x^*(k+1|k)-x_r(k+1)\right\|_Q^2 
$$