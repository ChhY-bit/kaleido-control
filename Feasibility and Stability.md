# 引言

上一篇已经介绍了经典无约束/有约束 MPC 的算法，该算法足以应对大多数精度要求不高、约束条件不复杂的场景。然而，上一篇博客所介绍的方法忽略了优化问题的递归可行性和控制系统的闭环稳定性问题，进而直接扔掉了终端代价 $\left\| \cdot \right\|_P^2$ ——这在理论上是不完美的，在实践中是有隐患的。

因此，本篇在上篇经典 MPC 的基础上，进一步分析递归可行性和闭环稳定性问题，并探讨与终端代价有关的问题，以及实际应用时的工程处理方法。

# MPC 的递归可行性问题

## 递归可行性的概念

MPC 本质上是在线求解优化问题 $$\min_u J \quad \text{s.t.}\; \mathrm{constraints}.$$ 其中，只要解序列 $\left\{u_0,u_1,\cdots,u_{N-1}\right\}$ 能够使得整个优化问题的约束条件得到满足，那么就被称为该优化问题的一个**可行解**；所有可行解中，能使目标函数 $J$ 取得最优值的，被称为**最优解**。显然，可行解可以不只有一个，并且不一定是最优解。

基于上述定义，**递归可行性 (Recursive Feasibility)** 是指：若优化问题在时刻 $k$ 有可行解，则在执行最优控制量后的任意时刻，优化问题依然能有可行解。保证递归可行性的根本目的是：
- 使滚动时域优化在理论上始终有解，这是 MPC 长期运行的前提；
- 为闭环稳定性证明提供关键工具。

在工程实现上，递归可行性还带来两个附带好处：
- 当求解器因偶然因素（如数值异常返回 NaN）未能给出解时，可将上一步解移位后得到的可行解作为备用控制序列执行，避免程序报错停机；
- 该可行解也可作为数值求解器的初始迭代点进行 **热启动 (Warm Start)**，提高求解效率。

## 终端控制器、终端集与终端约束

保证 MPC 递归可行性的最主要困难在于：尽管任意时刻都可以构建出相应的优化问题并付诸求解，但不同时刻的优化问题之间是彼此孤立的，当前时刻是否具有可行解完全取决于优化问题本身，且能否在限定时间内找到可行解，还取决于求解器性能——这就为 MPC 的长期滚动运行带来的隐患——毕竟没有理论依据能保证从此刻开始往后任意时刻的优化问题都有可行解。

解决上述问题，使用最广泛的方法是 **终端约束(Terminal Constraint)** 法。其核心思想是：在优化问题中加入针对 **终端状态（预测时域末端时刻系统状态）** 的特定约束，该约束要求系统在经过时长为预测时域的 MPC 作用后，进入名为 **终端集 (Terminal Set)** 的区域，即 $$ x(t+T|t)\in\Omega,$$ 其中终端集 $\Omega$ 具有这样的性质：

> 若 $ x(t_0) \in \Omega$，则 $\forall t>t_0, x(t)\in \Omega.$

这被称为终端集的 **正不变性 (Positive Invariance)**，或者说终端集合是一个 **不变集 (Invariant Set)**。换言之，只要系统状态进入了终端集的范围，那么它将不再离开这个集合——这就为寻找优化问题的可行解带来了思路：只要能找到一个状态约束集的子集，并能保证它为不变集，那么自系统状态进入该不变集起，往后任意时刻的状态约束都不会被违反，即 $x(t)\in\Omega\subseteq\mathcal{X}.$

然而，终端集 $\Omega$ 的不变性并不是与生俱来的，因为系统的演化必须要符合状态方程 $\dot{x}(t) = Ax(t)+Bu(t)$，而终端集如果要具有不变性，就必须在终端集区域内执行适当的控制操作，通过控制量的作用使系统状态继续停留在终端集内。换言之，终端集的不变性需要在相应时刻执行一定的控制作用来保证，该控制作用被称为 **终端控制器 (Terminal Controller)**：

$$
u(\tau|t) = \kappa(x(\tau|t)),\quad x(\tau|t)\in \Omega.
$$

对于 LTI 系统，终端控制器可以选取为状态反馈的形式：

$$
\kappa(x) = K_t x, 
$$

> $K_t$ 作为反馈增益，通常可以使用 **极点配置、LQR** 等已经能确保稳定性的方法来获得。

这样，如果能选取合适的 $K_t$ ，使得闭环系统矩阵 $A_t := A + BK_t$ 是 Hurwitz 的 ，那么 Lyapunov 方程

$$
S_tA_t + A_t^\top S_t = -Q_t,\; Q_t\succ 0

$$

必有对称正定的解 $S_t$，并且由该解构成的 $V(x) = x^\top S_t x$ 是一个 Lyapunov 函数，即 $\forall x\neq 0,\;\dot{V}(x)<0.$

如果选取 $\alpha>0$，并假设某一时刻 $t_0$ 时的函数值为

$$
x^\top(t_0) S_t x(t_0) = \alpha, 
$$

那么由 $\dot{V}(x)<0$ 可知

$$
\forall t>t_0, \quad x^\top(t) S_t x(t) < \alpha.
$$

由此构成的集合

$$
\Omega = \left\{ x\in\mathbb{R}^n \big| x^\top S_t x \le \alpha \right\}

$$

就是一个不变集。
> **需要注意的是：** 由于终端控制器 $\kappa(x) = K_t x$ 的执行必须满足约束，导致相应系统状态 $x$ 的允许范围不能是任意的，因此必须在 $\kappa(x) \in \mathcal{U}$ 的条件下，选择合适的 $\alpha$ 值来确定$\Omega.$

当不变集确定以后，我们知道：只要系统状态落入 $\Omega$，并且从此开始施加终端控制器 $\kappa(\cdot)$，那么系统状态必将维持在 $\Omega \subseteq \mathcal{X}$ 之内，且终端控制量也在控制约束集 $\mathcal{U}$ 之内——换言之，此时的控制量 $u(\tau|t)$ 与状态量 $x(\tau|t)$ 均是不违反 MPC 优化问题约束的。

下面将说明，只要 MPC 优化问题

$$
\begin{aligned}
 \min_{u(\tau|t)} \quad & J(x(t),u(t)) = \left\| x(t+T|t) \right\|^2_P + \int_{t}^{t+T} {\left\| x(\tau|t) \right\|^2_Q + \left\| u(\tau|t) \right\|^2_R} \mathrm{d}\tau, \\
\text{s.t.} \quad & \dfrac{\mathrm{d}}{\mathrm{d}\tau}x(\tau|t) = Ax(\tau|t) + Bu(\tau|t), \quad \tau\in [t,t+T],\\
& x(t|t) = x(t), \\
& x(\tau|t) \in \mathcal{X}, \\
& u(\tau|t) \in \mathcal{U}, \\
& x(t+T|t)\in \Omega
\end{aligned}

$$

在当前时刻可行（有解），那么在滚动优化机制下此后任意时刻也都可行。

假设 $t$ 时刻已经求得 MPC 优化问题的控制序列

$$
U^* = \left[ u^*_0,u^*_1,\cdots,u^*_{N-1} \right]^\top,
$$

对应于连续时域中的分段常值特性：

$$
u^*(\tau|t) = u^*_{i}, \quad i = \lfloor \dfrac{\tau-t}{\delta} \rfloor ,\quad \tau \in [t,t+T).
$$

那么在 $t+\delta$ 时刻，其控制量 $u(\tau|t+\delta),\;\tau\in[t+\delta,t+\delta+T)$ 同样应当由 $N$ 段函数构成。对于这样的控制序列，我们可以令其前 $N-1$ 段值沿用自上一时刻的常值 $U^*$，而最后一段采用状态反馈的终端控制器 $\kappa(x)$，并将这样组成的控制量记作 $\bar{u}(\tau|t+\delta)$，即

$$
\bar{u}(\tau|t+\delta) = 
\begin{cases}
u^*(\tau|t), & \tau \in [t+\delta,t+T)\\
K_t x(\tau|t+\delta), & \tau \in [t+T,t+\delta+T)
\end{cases}

$$

这个解就是 MPC 优化问题的一个可行解，因为：
- $\tau \in [t+\delta,t+T)$ 时，控制量完全沿用自已经求出的最优解，对于理想系统，状态量和控制量都会沿着“预测”出的方向演变，期间的控制约束 $\mathcal{U}$ 和状态约束 $\mathcal{X}$ 必然是被满足的（否则上一时刻的 MPC 问题就根本无法求出解）；
- $\tau \in [t+T,t+\delta+T)$ 时，已经超出了预测时域的作用范围，没有现成的不违反约束的控制量可用，此时切换为精心设计好的终端控制器 $\bar{u} = \kappa(x) = K_t x$，同样可以保证 $\bar{u}(\tau|t+\delta)\in\mathcal{U}$；并且由于优化问题中已经靠终端约束 $x(t+T|t+\delta) = x(t+T|t)\in \Omega$ 保证了 $\tau = t+T$ 时刻系统状态进入了由终端控制器所保证的不变集，因此在后续 $\tau\in[t+T,+\infty)$ 时必然有 $x(\tau|t+\delta)\in\Omega\subseteq\mathcal{X}.$

以上构造表明，$t+\delta$ 时刻存在可行解。由归纳法可知，对任意后续时刻 $t+2\delta,t+3\delta,\dots$ 均能构造出可行候选，因此滚动优化过程中优化问题始终保持可行——即递归可行性成立。

具体而言，记 $t_k := t+k\delta,\;k = 0,1,2,\cdots$，那么对于任意 $t_k$，后一时刻 $t_{k+1}$ 时都可以构造出可行解：

$$
\bar{u}(\tau|t_{k+1}) =
\begin{cases}
u^*(\tau|t_k), & \tau \in [t_{k+1},t_k+T) \\
K_t x(\tau|t_{k+1}), & \tau \in [t_k+T,t_{k+1}+T)
\end{cases}
$$

其中 $ u^*(\tau|t_k)$ 为上一时刻 $t_k$ 时真正由求解器求解 MPC 优化问题所得的最优解，这样所得的可行解均是可行的（不超出$\mathcal{U}$，也不会使状态超出$\mathcal{X}$）。

递归可行性的工程价值正在于此：它保证我们手中**始终握有一个不违反约束的备用解**。上述移位候选解可由上一时刻的解和终端控制器构造，无需等待当前时刻的优化求解。因此，当求解器因超时或数值异常未能返回新的解时，可以执行该候选解，维持约束可行性。不过，可行性本身并不保证控制性能或闭环稳定性，后者还需要进一步的代价下降条件。相比之下，执行不可行的控制量则可能违反系统的安全或运行限制，造成更为严重的后果。

---

# MPC 的闭环稳定性问题

上述分析说明了MPC优化问题的迭代可行性在引入终端约束、终端控制器之后是可以被保证的。以下将讨论 MPC 的稳定性问题。

## 稳定性问题产生的原因

在传统控制（例如 PID）中，控制律是显式的代数表达式，将其代入系统方程即可直接分析闭环稳定性。然而，对于带有约束的 MPC 问题，控制量是数值求解得到的序列，它由数值求解器在每个采样时刻独立求解优化问题给出，而非显式的函数关系。因此，系统在 MPC 作用下是否仍然稳定，无法像传统方法那样直接获知。

更关键的是，MPC 的预测时域是有限的：每一时刻求解出的控制量仅保证在当前预测窗口 \([t_k,t_k+T)\) 内代价函数的最优性，却无法天然保证前后两拍求解所对应的代价函数继续递减。这意味着，即使每一拍都成功求解了优化问题，代价函数在整个滚动过程中仍可能出现振荡甚至发散——有限时域的最优性并不自动蕴含闭环稳定性。

如果能够做到对无穷远的未来进行预测和优化，得到整个无穷时域上的最优解：

$$
\min J_\infty (t_k) = \int_{t_k}^{\infty} {\left\| x(\tau|t_k) \right\|^2_Q + \left\| u(\tau|t_k) \right\|^2_R} \mathrm{d}\tau,
$$

其中$Q,R$均正定，那么确保整个控制过程的稳定性将变得比较容易（例如证明$J_\infty$ 的时间导数在$x\neq0$时为负）。但是，在实际运行过程中，特别是对于有约束条件的问题，几乎不可能做到无穷时域的求解，预测时域长度通常只能取为有限值。

因此，有必要针对有限时域 MPC 的稳定性进行专门的分析和保证。

## 处理稳定性：利用终端代价

在 $t_k$ 时刻，假设已经求得最优控制量 $u^*(\tau|t_k)$ 及其作用下的状态量 $x^*(\tau|t_k)$ ，对应的代价函数为：

$$
J^*(t_k) = \left\| x^*(t_k+T|t_k) \right\|^2_P + \int_{t_k}^{t_k+T} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau.
$$

对于此后一个控制周期内的任意时刻 $\sigma\in(t_k,t_{k+1}]$ ，由于系统状态受到终端集的约束，因而具备递归可行性，故在其预测时域区间 $[\sigma,\sigma+T)$ 内可取控制量为：

$$
\bar{u}(\tau|\sigma) = 
\begin{cases}
u^*(\tau|t_k), & \tau \in [\sigma,t_{k}+T) \\
K_t \bar{x}(\tau|\sigma), & \tau \in [t_{k}+T,\sigma+T)
\end{cases}
$$

该控制量由两部分构成：第一部分均取自控制时刻 $t_k$ 时的优化结果；第二部分超出了预测时域，因而取自终端控制器。这意味着，在理想条件下，系统将按照“预测”结果运行，即

$$
\bar{x}(\tau|\sigma) = x^*(\tau|t_k),\quad \tau \in [\sigma,t_k+T].
$$

因此，对应的代价函数为：

$$
\begin{aligned}
\bar{J}(\sigma) &= \left\| \bar{x}(\sigma+T|\sigma) \right\|^2_P + \int_{\sigma}^{\sigma+T} \left\| \bar{x}(\tau|\sigma) \right\|^2_Q + \left\| \bar{u}(\tau|\sigma) \right\|^2_R \mathrm{d}\tau \\
&= \left\| \bar{x}(\sigma+T|\sigma) \right\|^2_P + \int_{\sigma}^{t_{k}+T} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau \\
& \quad + \int_{t_{k}+T}^{\sigma+T} \left\| \bar{x}(\tau|\sigma) \right\|^2_Q + \left\| K_t \bar{x}(\tau|\sigma) \right\|^2_R \mathrm{d}\tau
\end{aligned}
$$

选取待定 Lyapunov 函数

$$
V(x) = x^\top P x,
$$

由于终端闭环系统矩阵 $A_t:=A+BK_t$ 是 Hurwitz 的，因此 Lyapunov 方程

$$
PA_t + A_t^\top P = -Q_p,\quad Q_p\succ 0
$$

必有唯一正定解$P$，则当 $x\in\Omega$ 时，对于 $\dot{x} = A_t x$，有

$$
\forall x\neq 0, \quad \dot{V}(x) = x^\top \left( PA_t + A_t^\top P \right)x = -x^\top Q_p x < 0.
$$

上式两边同时从 $t_1$ 到 $t_2$ 积分，可得

$$
V(x(t_2))-V(x(t_1)) = -\int_{t_1}^{t_2} \left\| x(\tau) \right\|^2_{Q_p} \mathrm{d}\tau,
$$

对于 $t_k$ 时刻最优解、$\sigma$ 时刻可行解作用下的系统状态 $x^*(\cdot|t_k),\;\bar{x}(\cdot|\sigma)$，令上式中的 $t_1 = t_k+T,\;t_2 = \sigma + T$，则有

$$
\begin{aligned}
\left\| \bar{x}(\sigma+T|\sigma) \right\|^2_{P} - \left\| x^*(t_k+T|t_k) \right\|^2_{P} &= \left\| \bar{x}(\sigma+T|\sigma) \right\|^2_{P} - \left\| \bar{x}(t_k+T|\sigma) \right\|^2_{P} \\
&= -\int_{t_k+T}^{\sigma+T} \left\| \bar{x}(\tau|\sigma) \right\|^2_{Q_p} \mathrm{d}\tau
\end{aligned}
$$

将不同时刻的代价 $J^*(t_k)$ 和 $\bar{J}(\sigma)$ 作差，得到

$$
\begin{aligned}
\bar{J}(\sigma) - J^*(t_k) &= 
\left\| \bar{x}(\sigma+T|\sigma) \right\|^2_P - \left\| x^*(t_k+T|t_k) \right\|^2_P \\
& \quad - \int_{t_k}^{\sigma} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau \\
& \quad + \int_{t_{k}+T}^{\sigma+T} \left\| \bar{x}(\tau|\sigma) \right\|^2_Q + \left\| K_t \bar{x}(\tau|\sigma) \right\|^2_R \mathrm{d}\tau \\
&= -\int_{t_k+T}^{\sigma+T} \left\| \bar{x}(\tau|\sigma) \right\|^2_{Q_p} \mathrm{d}\tau \\
& \quad - \int_{t_k}^{\sigma} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau \\
& \quad + \int_{t_{k}+T}^{\sigma+T} \left\| \bar{x}(\tau|\sigma) \right\|^2_{Q+K_t^\top R K_t} \mathrm{d}\tau \\
&= -\int_{t_k}^{\sigma} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau \\ 
& \quad + \int_{t_{k}+T}^{\sigma+T} \bar{x}(\tau|\sigma)^\top \left(Q+K_t^\top R K_t-Q_p\right)\bar{x}(\tau|\sigma) \mathrm{d}\tau
\end{aligned}.
$$

由于二次型的权重矩阵 $Q+K_t^\top R K_t-Q_p$ 尚未保证半负定，无法保证 $ \bar{J}(\sigma) - J^*(t_k) $ 一定为负。要使代价函数递减，可以令

$$
Q+K_t^\top R K_t \preceq Q_p,
$$

> 从直观意义上来看，取 $ Q+K_t^\top R K_t \preceq Q_p $ 能够保证在超出预测时域的 $[t_k+T,\sigma+T]$ 这段时间内，终端代价的减少量（即$Q_p$对应的项），不低于终端控制作用下的代价增量。

则由于 $\sigma > t_k$ 可知

$$
\bar{J}(\sigma) - J^*(t_k)  \le -\int_{t_k}^{\sigma} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau \le 0.
$$

由于 $\bar{J}(\sigma)$ 是在可行解 $\bar{u}(\tau|\sigma),\;\tau\in(\sigma,\sigma+T]$ 作用下得到的，不考虑有限维参数化带来的最优性损失，客观上必然存在某个不差于 $\bar{J}(\sigma)$ 代价 $J^*(\sigma)$，即

$$
J^*(\sigma) \le \bar{J}(\sigma),
$$

进而有

$$
J^*(\sigma) - J^*(t_k) \le -\int_{t_k}^{\sigma} \left\| x^*(\tau|t_k) \right\|^2_Q + \left\| u^*(\tau|t_k) \right\|^2_R \mathrm{d}\tau \le 0.
$$

上式对任意 $\sigma\in(t_k,t_{k+1}]$ 都成立。当 $x(t_k)\neq0$ 时，右端积分为正，故 $J^*(\sigma)<J^*(t_k)$；特别地，$J^*$ 在相邻采样时刻严格下降。结合最优价值函数的正定性、适当的上下界及采样间状态演化性质，可证明原点在相应可行域内渐近稳定。

## 无穷时域优化角度：二次型准无限时域 MPC

上述利用终端代价保证稳定性的 MPC 方法，称为**二次型准无限时域 MPC**。接下来阐述其与无穷时域优化问题的关系。

可以将无穷时域的代价拆分为：

$$
\begin{aligned}
J_\infty (t_k) &= \int_{t_k}^{t_k+T} \left\|x(\tau|t_k)\right\|_Q^2 + \left\| u(\tau|t_k)\right\|_R^2 \mathrm{d}\tau \\
& \quad + \underbrace{\int_{t_k+T}^{\infty} \left\|x(\tau|t_k)\right\|_Q^2 + \left\| u(\tau|t_k)\right\|_R^2 \mathrm{d}\tau}_{\text{Infinity Tail}},
\end{aligned}
$$

其中，有限时域部分 $[t_k,t_k+T]$ 相对而言是容易处理的（因为其长度有限），而无限余项 (Infinity Tail) 部分处理起来则是比较困难的。

而在引入终端要素（终端约束、终端控制器、终端代价）后，在终端控制器 $u(\tau|t_k) = K_t x(\tau|t_k),\;\tau \in [t_k+T,+\infty)$ 的作用下，无穷余项 (Infinity Tail) 将变为：

$$
\int_{t_k+T}^{\infty} \left\|x(\tau|t_k)\right\|_Q^2 + \left\| u(\tau|t_k)\right\|_R^2 \mathrm{d}\tau = \int_{t_k+T}^{\infty} \left\lVert x(\tau|t_k) \right\rVert ^2_{Q+K_t^\top R K_t} \mathrm{d}\tau.
$$

根据前一节证明可知，$$ -\left\lVert x(\tau|t_k) \right\rVert ^2_{Q+K_t^\top R K_t} \ge -\left\lVert x(\tau|t_k) \right\rVert ^2_{Q_p}, $$

两边在 $[t_k+T,+\infty)$ 上积分可得（其中 $V(\cdot):=\left\| \cdot \right\|_P^2$，满足 $\dfrac{\mathrm{d}V(\cdot)}{\mathrm{d}t} = -\left\|\cdot\right\|_{Q_p}^2$）

$$
-\left(\int_{t_k+T}^{\infty} \left\|x(\tau|t_k)\right\|_Q^2 + \left\| u(\tau|t_k)\right\|_R^2 \mathrm{d}\tau\right) \ge \lim_{\tau\to\infty}V\left(x(\tau|t_k)\right) - V \left(x(t_k+T|t_k)\right),
$$

由于终端闭环矩阵 $A_t$ 为 Hurwitz，终端控制器作用下的状态渐近收敛到原点，因此 $\lim_{\tau\to\infty}V(x(\tau|t_k))=0$，那么有

$$
\int_{t_k+T}^{\infty} \left\|x(\tau|t_k)\right\|_Q^2 + \left\| u(\tau|t_k)\right\|_R^2 \mathrm{d}\tau \le V \left(x(t_k+T|t_k)\right) = \left\lVert x(t_k+T|t_k) \right\rVert _P^2,
$$

换言之，这里选定的终端代价 $\left\|\cdot\right\|^2_P$ 正是无穷余项的一个上界。不难发现，若取 $Q+K_t^\top R K_t = Q_p$ 则终端代价恰好等于无穷余项。



---

> 由此可见，上篇中取 $P = 0$，系统也能稳定运行，是因为……，是……的一个特例。

# 完善的 LTI-MPC

