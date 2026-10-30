# Linear Quadratic MPC

## 1 Problem Formulation

- Discrete system model: $$ x(k+1) = A x(k) + B u(k), $$ where $x(k)\in\mathbb{R}^{n},u(k)\in\mathbb{R}^{m},k\in\mathbb{Z}_{\ge 0}.$

- The regulation problem with linear-quadratic cost function: $$ J(k) = \| x(k+N|k) \|_P^2 + \sum_{i=0}^{N-1}\|x(k+1+i|k)\|_Q^2 + \|u(k+i|k)\|_R^2. $$
The LQ form above is suitable for: zero-stablization or trajectory tracking with $$ x_r(k) = Ax_r(k) + Bu_r(k)$$ and we can define $ x_e:=x-x_r$ such that $$ x_e(k+1) = Ax_e(k)+Bu(k). $$

- Constraints: 
    - State constraint: $x(k+i|k)\in\mathcal{X}$
    - Control constraint: $u(k+i|k)\in\mathcal{U}$
    - Terminal constraint: $x(k+N|k)\in\mathcal{X}_t$

- The MPC optimization problem:

$$
\begin{aligned}
    \min_{u(\cdot|k)}\quad & J(k) = \| x(k+N|k) \|_P^2 + \sum_{i=0}^{N-1}\|x(k+1+i|k)\|_Q^2 + \|u(k+i|k)\|_R^2. \\
    \text{s.t.}\quad & x(k|k) = x(k), \\
    & x(k+i+1|k) = Ax(k+i|k) + Bu(k+i|k), \quad i=0,1,\cdots N-1\\
    & x(k+i|k) \in\mathcal{X}, \quad i=1,2,\cdots N\\
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

Especially, for the terminal time instance, the following representation holds: 

$$
x(k+N|k) = A^Nx(k|k) + \Gamma_u U(k),
$$

where $$\Gamma_u = \begin{bmatrix} O_{n\times n(N-1)} & I_n \end{bmatrix} \mathcal{T}_u = \begin{bmatrix}A^{N-1}B & A^{N-2}B & \cdots & B \end{bmatrix}.$$

### 2.2 Rewrite the Cost

Define the extended wight matrices:

$$\mathcal{Q}:=I_N\otimes Q,\quad \mathcal{R} := I_N\otimes R,\quad \mathcal{S}:=I_{N-1}\otimes S.$$

> Take block-diagonal matrices like $\text{blkdiag}[Q(k+1),Q(k+2),\cdots,Q(k+N)]$ if regarding time-variant weights.

Then the cost has a standard QP formualation:

$$ \begin{aligned}
J(k) &= \| x(k+N|k) \|_P^2 + \|X(k)\|_{\mathcal{Q}}^2 + \|U(k)\|_{\mathcal{R}}^2 \\
&= \|\underbrace{A^Nx(k|k)}_{\Psi(k)} + \Gamma_u U(k)\|^2_P + \|\underbrace{\mathcal{T}_x x(k|k)}_{\Phi(k)s} + \mathcal{T}_u U(k)\|_{\mathcal{Q}}^2 + \|U(k)\|_{\mathcal{R}}^2 \\
&= 
\end{aligned}
$$