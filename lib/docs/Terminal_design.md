## 1 可行性问题

代价函数：

$$ \begin{aligned}
    J(k) =& \sum_{i=0}^{N-1}\left(\left\|y(k+i|k)\right\|_Q^2 + \left\|u(k+i|k)\right\|^2_{R} \right) + \\
    & \left\| y(k+N|k) \right\|^2_P + \sum_{i=1}^{N-1}\left\| \Delta u(k+i|k) \right\|^2_S.
\end{aligned} $$

其中 $\Delta u(k+i|k) := u(k+i|k)-u(k+i-1|k).$

假设已有满足各类约束的最优解：

$$ U^*(k) = \begin{bmatrix}
    u^*(k|k) \\ u^*(k+1|k) \\ \vdots \\ u^*(K+N-1|k)
\end{bmatrix}, $$

尽管只有 $u^*(k|k)$ 被实际执行，但剩余项仍可用于构造下一时刻的可行解，即 $$\bar{u}(k+i|k+1) = \begin{cases} u^*(k+i|k) , & i=1,2,\cdots,N-1 \\ \kappa(x^*(k+N|k)), & i = N
\end{cases}    $$ 或者写作

$$\bar{U}(k+1) = \begin{bmatrix}
    u^*(k+1|k) \\ u^*(k+2|k) \\ \vdots \\ u^*(k+N-1|k) \\ Kx^*(k+N|k)
\end{bmatrix}, $$

其中 $x^*(k+N|k) = A^N x(k|k) + \Gamma_u U^*(k)$；终端反馈增益 $K$ 可通过 LQR 或其他方法获得，只需保证闭环系统矩阵 $A_t:=A+BK$ 是 Schur 稳定的，即 Lyapunov 方程 $$A_t^\top M A_t - M = - N,\quad N^\top = N \succ 0$$ 有对称正定解 $M.$ 可以证明在相应的终端集

$$\mathcal{Y}_f = \left\{ y\in \mathbb{R}^p | y^\top M y \le \alpha ,\;\alpha \ge 0 \right\} \subseteq \mathcal{Y}$$ 

内实施终端控制时，该集合是不变集；其中 $\alpha$ 可用于调节集合范围的大小。

由于控制量受到集合 $\mathcal{U}$ 的约束，需确保终端控制 $Kx^*(k+N|k)\in\mathcal{U}$

在终端约束下，能保证 $$ y^*(k+N|k) = (C+DK)x^*(k+N|k) \in \mathcal{Y}_f,$$ 