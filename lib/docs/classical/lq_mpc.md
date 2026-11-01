# Linear Quadratic MPC

## 1 Problem Formulation

- Discrete system model: $$ x(k+1) = A x(k) + B u(k), $$ where $x(k)\in\mathbb{R}^{n},u(k)\in\mathbb{R}^{m},k\in\mathbb{Z}_{\ge 0}.$

- The regulation problem with linear-quadratic cost function: $$ J(k) = \left\| x(k+N|k) \right\|_P^2 + \sum_{i=0}^{N-1} \left(\|x(k+i|k)\|_Q^2 + \left\|u(k+i|k)\right\|_R^2\right), $$ where $P,Q,R$ are symetric positive definite matrices with proper dimension.

- Constraints: 
    - State constraint: $x(k+i|k)\in\mathcal{X}$
    - Control constraint: $u(k+i|k)\in\mathcal{U}$
    - Terminal constraint: $x(k+N|k)\in\mathcal{X}_t$

- The MPC optimization problem:

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

## 2 Optimal Solution 

### 2.1 Recursive Prediction

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

Represent them as 

$$
X(k) = \mathcal{T}_x x(k|k) + \mathcal{T}_u U(k),
$$

where

$$
\begin{aligned}
X(k) &:= \begin{bmatrix}
    x(k|k) \\ x(k+1|k) \\ x(k+2|k) \\ \vdots \\ x(k+N-1|k)
\end{bmatrix}, \quad
U(k) := \begin{bmatrix}
    u(k|k) \\ u(k+1|k) \\ u(k+2|k) \\ \vdots \\ u(k+N-1|k)
\end{bmatrix}, \\[2ex]
\mathcal{T}_x &:= \begin{bmatrix}
    I \\ A \\ A^2 \\ \vdots \\ A^{N-1}
\end{bmatrix}, \quad
\mathcal{T}_u := \begin{bmatrix}
    0 & 0 & 0 & \cdots & 0 & 0 \\
    B & 0 & 0 & \cdots & 0 & 0 \\
    AB & B & 0 & \cdots & 0 & 0 \\
    \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
    A^{N-2}B & A^{N-3}B & A^{N-4}B & \cdots & B & 0
\end{bmatrix}.
\end{aligned}
$$

Especially, for the terminal time instance, the following representation holds: 

$$
x(k+N|k) = A^Nx(k|k) + \Gamma_u U(k),
$$

where $$\Gamma_u = \begin{bmatrix}A^{N-1}B & A^{N-2}B & \cdots & B \end{bmatrix}.$$

### 2.2 Rewrite the Cost

Define the extended wight matrices:

$$\mathcal{Q}:=I_N\otimes Q,\quad \mathcal{R} := I_N\otimes R.$$

> Take block-diagonal matrices like $\text{blkdiag}[Q(k),Q(k+1),\cdots,Q(k+N-1)]$ if regarding time-variant weights.

Then the cost has a standard QP formualation:

$$ \begin{aligned}
J(k) &= \| x(k+N|k) \|_P^2 + \|X(k)\|_{\mathcal{Q}}^2 + \|U(k)\|_{\mathcal{R}}^2 \\
&= \|\underbrace{A^Nx(k|k)}_{\psi(k)} + \Gamma_u U(k)\|^2_P + \|\underbrace{\mathcal{T}_x x(k|k)}_{\Phi(k)} + \mathcal{T}_u U(k)\|_{\mathcal{Q}}^2 + \|U(k)\|_{\mathcal{R}}^2 \\
&= U^\top(k) \underbrace{\left( \Gamma_u^\top P \Gamma_u + \mathcal{T}_u^\top \mathcal{Q} \mathcal{T}_u + \mathcal{R} \right)}_{H} U(k) + \\
& \quad \; 2\underbrace{\left( \psi^\top(k) P\Gamma_u + \Phi^\top(k) \mathcal{Q} \mathcal{T}_u \right)}_{f^\top(k)} U(k) + \underbrace{\psi^\top(k) P \psi(k) + \Phi^\top(k)\mathcal{Q}\Phi(k)}_{c(k)}.
\end{aligned}
$$

### 2.3 Rewrite the Constraints

Assume the polyhedral constraints：

$$
x_{\min}(k+i) \le x(k+i|k) \le x_{\max}(k+i), \quad i=0,1,\cdots,N-1 \\
u_{\min}(k+i) \le u(k+i|k) \le u_{\max}(k+i), \quad i=0,1,\cdots,N-1 \\
x_{t,\min} \le E_t x(k+N|k) \le x_{t,\max}.\\
$$

Define

$$
X_{lb}(k):=
\begin{bmatrix}
x_{\min}(k) \\ x_{\min}(k+1) \\ \vdots \\ x_{\min}(k+N-1)
\end{bmatrix},\quad
X_{ub}(k):=
\begin{bmatrix}
x_{\max}(k) \\ x_{\max}(k+1) \\ \vdots \\ x_{\max}(k+N-1)
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

> If the constraint bound is constant，we can define $ X_{lb} := \mathbf{1}_N \otimes x_{\min} $, so do the others.

Then the constraints can be rewritten as

$$
X_{lb}(k) \le \mathcal{T}_x x(k|k) + \mathcal{T}_u U(k) \le X_{ub}(k), \\
U_{lb}(k) \le U(k) \le U_{ub}(k), \\
x_{t,\min} \le E_t\Gamma_u U(k) + E_tA^N x(k|k) \le x_{t,\max}.\\
$$

Transform into the standard QP form:

$$
\underbrace{\begin{bmatrix}
U_{lb}(k) \\ X_{lb}(k)-\Phi(k) \\ x_{t,\min} - E_t\psi(k)
\end{bmatrix}}_{G_{lb}(k)} \le 
\underbrace{\begin{bmatrix}
I_{mN} \\ \mathcal{T}_{u} \\ E_t\Gamma_u
\end{bmatrix}}_{F} U(k) \le
\underbrace{\begin{bmatrix}
U_{ub}(k) \\ X_{ub}(k)-\Phi(k) \\ x_{t,\max} - E_t\psi(k)
\end{bmatrix}}_{G_{ub}(k)}.
$$

如果 QP 求解器要求的形式为 $F U \le G $，那么可以再通过引入负号的方式加以处理。

## 3 State Observer

It is necessary to aquire the initial state $x(k|k) = x(k)$.

If $x(k)$ cannot be aquired directly, while system $(A,C)$ is detectable, an observer can be developed:

$$ \begin{cases}
\hat{x}(k+1) = A \hat{x}(k) + Bu(k) + L\left(y(k)-\hat{y}(k)\right), \\ \hat{y}(k) = C\hat{x}(k) + Du(k),
\end{cases}$$

so that $x(k)$ can be replaced by its estimation $\hat{x}(k)$.

## 4 The LQ-MPC Algorithm

**Given:**
- System model $A,B$
- Cost weight matrices $Q,R$
- Prediction horizon $N$
- Terminal ingredients $E_t,x_{t,\min/\max},P$
- Constraints $x_{\min / \max}(k),u_{\min / \max}(k)$

**Offline:**
- Construct $\mathcal{Q},\mathcal{R}$
- Construct $\mathcal{T}_x,\mathcal{T}_u,\Gamma_u$
- Construct $H,F$

**Online:**
- Aquire $x(k|k)$ (or $\hat{x}(k|k)$)
- Extract $X_{lb/ub}(k),U_{lb/ub}(k)$
- Calculate $\Phi(k),\psi(k),G_{lb/ub}(k).$
- Calculate $f^\top(k),c(k).$
- Solve $\min_{U} U^\top H U + 2f^\top(k) U + c(k),\;\text{s.t.}G_{lb} \le FU \le G_{ub}.$
- Take the first $m$ elements of $U(k)$ as the control input $u(k)$
- (Optional) Move $U(k)$ and use terminal controller as the feasible solution (under certain conditions) and warm-start the next optimization.

## 5 Beyond Zero-Stabilization 

When exsiting a reference trajectory $x_r(k)$ satisfying $$ x_r(k+1) = Ax_r(k) + Bu_r(k),$$ the LQ-MPC problem can be reformulated as:

$$ x_e(k+1) = Ax_e(k)+Bu_e(k), $$

where $ x_e(k):=x(k)-x_r(k),u_e(k):=u(k)-u_r(k).$ 

In MPC problem, the error vectors are defined as

$$ x_e(k+i|k) = x(k+i|k) - x_r(k+i),\; u_e(k+i|k) = u(k+i|k) - u_r(k+i).$$

Meanwhile, the form of  other ingredients of the MPC Optimizing problem can remain the same by only replacing $x,u$ with $x_e,u_e$, except for the constraints, which should be modified as: $$
x_e(k+i|k)\in\mathcal{X}-x_r(k+i), \\
u_e(k+i|k)\in\mathcal{U}-u_r(k+i), \\
x_e(k+N|k)\in\mathcal{E}_t.
$$
