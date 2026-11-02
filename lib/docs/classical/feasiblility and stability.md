# Feasibility and Stability of LQ-MPC

Consider the standard LQ-MPC optimization problem:

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


## 1 Recursive Feasibility

### 1.1 Feasible Solution

Assume that the optimal control squence $$ U^*(k) = \begin{bmatrix}
u^*(k|k) \\ u^*(k+1|k) \\ u^*(k+2|k) \\ \vdots \\ u^*(k+N-1|k)
\end{bmatrix} $$ is available at time instant $k$, then we can construct feasible control input at time instant $k+1$ without solving the optimization by shifting the optimal control squence:

$$\bar{u}(k+i|k+1) = \begin{cases} u^*(k+i|k), & i=1,2,\cdots,N-1 \\ 
\kappa(\bar{x}(k+N|k+1)), & i=N \end{cases}$$

The $\kappa(\cdot)$ is called the **terminal controller**, which can be selected as a state feedback law $$\kappa(\bar{x}(k+N|k+1)) = K_t \bar{x}(k+N|k+1),$$ where $K_t$ is a constant gain matrix.

### 1.2 Terminal Invariant Set

The feedback gain $K_t$ usually determined by the LQR or Pole-Placement such that the *Lypapunov* equation holds:

$$A_t^\top M A_t - M = -W,\quad W^\top = W \succ 0$$

with the unique symetric positive definite solution $M$, where $A_t:=A+BK_t$ is *Schur* and nonsingular.

It can be proved that the set 

$$\mathcal{X}_t = \left\{x\in\mathbb{R}^{n} | x^\top M x \le \alpha , \;\alpha>0 \right\} \subseteq \mathcal{X} $$

is an invariant set after exicuting the feasible control input, with the nominal predictive state at time instant $k+1$ :

$$
\bar{x}(k+i|k+1) = \begin{cases} x^*(k+i|k), & i = 1,2,\cdots,N \\
A_t x^*(k+N|k), & i = N+1\end{cases} 
$$

Then if the control constraints are satisfied meanwhile, i.e. $K_t x^*(k+N|k) \in \mathcal{U},\forall x^*(k+N|k)\in\mathcal{X}_t$, or the linear mapping form $$\mathscr{K}(\mathcal{X}_t) = \left\{ K_t x \; \big| \; x\in\mathcal{X}_t \right\} \subseteq \mathcal{U}.$$


## 2 Closed-Loop Stability

### 2.1 Optimal Cost and Feasible Cost

The optimal cost at time instant $k$ after executing the optimal control sequence $U^*(k)$ is

$$J^*(k) = \left\| x^*(k+N|k) \right\|_P^2 + \sum_{i=0}^{N-1} \left(\|x^*(k+i|k)\|_Q^2 + \left\|u^*(k+i|k)\right\|_R^2\right).$$

The feasible cost at time instant $k+1$ after executing the feasible control sequence $\bar{U}(k+1)$ is

$$\begin{aligned}
\bar{J}(k+1) &= \left\| \bar{x}(k+1+N|k+1) \right\|_P^2 + \sum_{i=1}^{N} \left(\|\bar{x}(k+i|k+1)\|_Q^2 + \left\|\bar{u}(k+i|k+1)\right\|_R^2\right) \\
&= \left\| A_tx^*(k+N|k) \right\|_P^2 + \sum_{i=1}^{N-1} \left(\|x^*(k+i|k)\|_Q^2 + \left\|u^*(k+i|k)\right\|_R^2\right) + \\ 
& \quad \; \| x^*(k+N|k)\|_Q^2 + \left\|K_t x^*(k+N|k)\right\|_R^2\\
&= \| x^*(k+N|k) \|_{A_t^\top P A_t + Q + K_t^\top R K_t}^2 + \sum_{i=1}^{N-1} \left(\|x^*(k+i|k)\|_Q^2 + \left\|u^*(k+i|k)\right\|_R^2\right).
\end{aligned}$$

where $\|z\|_{M}^2$ denotes the quadratic form $z^\top M z$ rather than 2-norm.

### 2.2 Terminal Cost and Its Weight

Take the difference of the two cost functions, we have

$$
\bar{J}(k+1) - J^*(k) = \| x^*(k+N|k) \|_{A_t^\top P A_t - P + Q + K_t^\top R K_t}^2 - \left(\|x^*(k|k)\|_Q^2 + \left\|u^*(k|k)\right\|_R^2\right)
$$

If the wight matrix of terminal cost can be chosen properly so that the *Lyapunov* inequation holds:

$$A_t^\top P A_t - P \preceq -\left(Q + K_t^\top R K_t\right),$$

then the difference of the two cost functions will be 

$$\bar{J}(k+1) - J^*(k) \le - \left(\|x^*(k|k)\|_Q^2 + \left\|u^*(k|k)\right\|_R^2\right) \le 0$$

If resolving the receding optimization problem, it is obvious that the optimal cost will be more optimal than the feasible cost, i.e. $J^*(k+1)\le \bar{J}(k+1).$ Hence $$J^*(k+1) - J^*(k) \le - \left(\|x^*(k|k)\|_Q^2 + \left\|u^*(k|k)\right\|_R^2\right) \le 0.$$ The latter equal sign holds if and only if $x^*(k|k) = 0$ and $u^*(k|k) = 0$, which means the sysem is then can be proved asymptotically stable.

Therefore, the inequality above of $P$ holding is the sufficient (but not necessary) condition for the system to be asymptotically stable.

## 3 Algorithm for Parameter Tuning

**Given:**
- System model $A,B$
- Constraints $\mathcal{X}, \mathcal{U}$
- Symetric positive-definite matrix $W$
- Weight matrices $Q,R$
- Terminal wight margin $\varepsilon > 0$

**Calculate:**
- Determine $K_t$ and calculate $A_t$
- Solve for $M$
- Optimize the radius: 

$$\begin{aligned} \max \quad & \alpha, \\ \text{s.t.}\quad & \mathcal{X}_t \subseteq \mathcal{X}, \\ & \mathscr{K}(\mathcal{X}_t) \subseteq \mathcal{U}. \end{aligned}$$

- Determine $P$ by $$A_t^\top P A_t - P = -\left(Q + K_t^\top R K_t + \varepsilon I_n \right)$$ 

**Return:**
- Terminal set $\mathcal{X}_t$
- Terminal controller gain $K_t$
- Terminal cost weight $P$

## 4 Quadratic to Linear Gap

The quadratic constraints $$ x^\top(k+N|k) M x(k+N|k) \le \alpha $$ cannot be inputed into the LQ-MPC framework directly. If only the linear-quadratic solver available (such as OSQP), it is necessary to replace the quadratic constraint with a linear/polyhedral constraint $$x_{t.\min} \le E_t x(k+N|k) \le x_{t,\max}.$$

### Algorithm of Determining Linear Terminal Constraint

**Input:**
- Closed-loop system matrix $A_t$
- Shape of the terminal set $M$
- Radius of the terminal set $\alpha$

**Initialize:**
- $r = [r_1,r_2,\cdots,r_n]^\top, r_i = \sqrt{\alpha\cdot [M^{-1}]_{ii}}.$
- $\mathcal{S} = \{x\in\mathbb{R}^n | -r \le x \le r\}$
- $\mathscr{A}(S):=\{A_t x | x\in\mathcal{S} \}$

**Algorithm:**
$k \leftarrow 0,\;G\leftarrow M,\;E_t\leftarrow I_n$
**while** $ r^\top \mathrm{abs}(G) r > \alpha $ **do:**
$G \leftarrow A_t^\top G A_t$
$k \leftarrow k+1$
$ E_t \leftarrow \mathrm{col}(E_t,A_t^{-k})$
**end while**
$ x_{t,\max} \leftarrow \boldsymbol{1}_{k+1} \otimes r $
$ x_{t,\min} \leftarrow -\boldsymbol{1}_{k+1} \otimes r $


**Output:** $E_t,x_{t,\max},x_{t,\min}$

