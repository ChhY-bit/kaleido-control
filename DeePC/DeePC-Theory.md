# Data-enalbed Predictive Control (DeePC)

## 1 Problem Formulation

### 1.1 The Actual System Model

Assume that the actual LTI system is

$$ \begin{aligned}
x(k+1) &= Ax(k) + Bu(k), \\
y(k) &= Cx(k) + Du(k),    
\end{aligned} $$

where $x\in\mathbb{R}^n,u\in\mathbb{R}^{m},y\in\mathbb{R}^p$ are state, input and output vector, respectively; $A,B,C,D$ are the matrices with proper dimensions. However, the matrices $(A,B,C,D)$ are unknown, leading to unavailability of the sysytem model. Therefore, the model-based method (e.g. MPC) are unable to reach full performance.

### 1.2 Available Data in Control

Model identification $(\hat{A},\hat{B},\hat{C},\hat{D})$ is commonly of high cost, while the input-output pairs $\{u_k^d,y_k^d\}_{k=0}^{S-1}$ with scale $S$ are much more easy to acquire:

$$\begin{aligned}
u^d = \mathrm{col} \left[ u_0^d,u_1^d,\cdots,u_{S-1}^d \right] \in \mathbb{R}^{mS},\\
y^d = \mathrm{col} \left[ y_0^d,y_1^d,\cdots,y_{S-1}^d \right] \in \mathbb{R}^{pS}.
\end{aligned}$$

It is required that $\{u_k^d,y_k^d\}_{k=0}^{S-1}$ is a trajetory of the unknown system, which actually satisfies $$ y^d_{k+1} = CA x^d_{k} + CB u^d_k + D u^d_{k+1}, \quad \forall k \in\mathbb{I}_{[0,S-2]}. $$ Then all of other pairs $\{\hat{u}_k,\hat{y}_k\}_{k=0}^{N-1}$ with horizon length $N$ belonging to this trajectory can be represented by $\{u_k^d,y_k^d\}_{k=0}^{S-1}$ instead of the actual model $(A,B,C,D)$, if the input sequence $u^d$ meet specific condition called **persistently exciting**.

In details, there exists the coordinates $\alpha \in \mathbb{R}^{S-N+1}$ such that the sequence

$$\begin{aligned}
\hat{u} = \mathrm{col} \left[ \hat{u}_0,\hat{u}_1,\cdots,\hat{u}_{N-1} \right] \in \mathbb{R}^{mN},\\
\hat{y} = \mathrm{col} \left[ \hat{y}_0,\hat{y}_1,\cdots,\hat{y}_{N-1} \right] \in \mathbb{R}^{pN},
\end{aligned}$$

can be represented by 

$$\begin{bmatrix}
    \mathcal{H}_N(u^d) \\ \mathcal{H}_N(y^d)
\end{bmatrix} \alpha = \begin{bmatrix}
    \hat{u} \\ \hat{y}
\end{bmatrix},$$

where $\mathcal{H}_N(\cdot)$ is the $N$-order *Hankel* matrix:

$$
\mathcal{H}_{N}(u^d) = \begin{bmatrix}
u^d_0 & u^d_1 & \cdots & u^d_{S-N} \\
u^d_1 & u^d_2 & \cdots & u^d_{S-N+1} \\
\vdots & \vdots & \ddots & \vdots \\
u^d_{N-1} & u^d_{N} & \cdots & u^d_{S-1} \\
\end{bmatrix}\in\mathbb{R}^{mN\times (S-N+1)},
$$

the similar as $\mathcal{H}_N(y^d)\in \mathbb{R}^{pN\times (S-N+1)}.$

### 1.3 DeePC Scheme

The core difference between DeePC and MPC is the way of prediction:
- **MPC:** Predict the next $N$ steps by known model $(A,B,C,D)$
- **DeePC:** Predict the next $N$ steps by using data-set $\{u_k^d,y_k^d\}_{k=0}^{S-1}$ and a section of initial trajectory $\{u_k,y_k\}_{k=-T}^{-1}$ with the length of $T$, whose vector forms are 

$$\begin{aligned}
u_{\text{ini}} = \mathrm{col} \left[ u_{(-T)},u_{(-T+1)},\cdots,u_{(-1)} \right] \in \mathbb{R}^{mT},\\
y_{\text{ini}} = \mathrm{col} \left[ y_{(-T)},y_{(-T+1)},\cdots,y_{(-1)} \right] \in \mathbb{R}^{pT},
\end{aligned}$$

Note that $\{u_k,y_k\}_{k=-T}^{-1}$ are all actual data during control process, neighter offline data-set nor virtual prediction.

Extend the representing interval to $L=T+N$, then the classic DeePC scheme at real time $t$ is

$$
\begin{aligned}
\min_{U(t)} \quad & J(t) = \sum_{k=0}^{N-1}\ell \left( \hat{u}_k(t),\hat{y}_k(t) \right) \\ 
\text{s.t.} \quad & \begin{bmatrix}
    \mathcal{H}_{L}(u^d) \\ \mathcal{H}_{L}(y^d)
\end{bmatrix} \alpha(t) = \begin{bmatrix}
    u_{\text{ini}}(t) \\ \hat{u}(t) \\ y_{\text{ini}}(t) \\ \hat{y}(t)
\end{bmatrix}, \\
& \hat{u}_k(t) \in \mathcal{U},\; \hat{y}_k(t) \in \mathcal{Y},\;k\in\mathbb{I}_{[0,N-1]},
\end{aligned}
$$

where $$U(t):=\begin{bmatrix} \alpha^\top(t) & \hat{u}^\top(t) & \hat{y}^\top(t) \end{bmatrix}^\top,\quad \alpha \in \mathbb{R}^{S-L+1}$$ is the decision variable vector, which indicated that the optimization problem is in a sparse form.

## 2 Usage of DeePC

### 2.1 Conditions before Solving

#### 2.1.1 Persistenly Excitation

If the length of the sequence to be represented is $L$ (actually $L = N+T$), it is sufficient that the historical input sequence $u^d$ is **persistenly exciting** of order $L+n$, which means

$$
\mathrm{rank} \left( \mathcal{H}_{L+n}(u^d) \right) = m(L+n).
$$

Practically, due to the noise when sampling data-set, the *Hankel* matrix can be of full rank too easily, which is not truely persistenly exciting in fact.

Alternatively, the Singular Value Decomposition (SVD) technique can be used to identify the actual effective order of the persistenly excitation.

We can firstly collect a large amount of data with the length of $K$ and decompose  them by

$$V = \mathrm{SVD}\left(\mathcal{H}^\top_{K}(u^d) \mathcal{H}_{K}(u^d)\right)$$

Then we can take the largest $L+n$ elements of $V$, where the $(L+n+1)$-th value is much less than the $L+n$-th.

Another challenge is to determine the system order $n$ (and its lag $l$, see <u> Data-Enabled Predictive Control: In the Shallows of the DeePC, ECC, 2019</u>) because of its unknowns. An acceptable method is using a sufficiently large upper bound $\bar{n} \ge n$ to replace $n$ and $l$.

#### 2.1.2 Dimension Requirement

- **Hankel matrix:** $S-(L+n)+1\ge m(L+n)$, i.e. $S \ge (m+1)(L+n)-1$, for the full rank requirement.
- **Representative window:** $S \ge (m+1)(T+N+n)-1$, and $T\ge l$ due to the obserbavility principle of linear system.

- In conclusion, the scale $S$ of the data-set should be  $$S \ge (m+1)(T+N+n)-1,$$ where
    - $m$ is the dimension of input;
    - $T$ is the length of initial window, and usually requring $T\ge n$;
    - $N$ is the predictive horizon;
    - $n$ is the dimension of system;
    - $L+n = T+N+n$ is usually estimated by SVD technique.

### 2.2 Calculation and Solution

#### 2.2.1 Rewrite the constraints

Decompose the *Hankel* matrices into

$$
\begin{bmatrix}
    \mathcal{H}_{N+T}(u^d) \\ \mathcal{H}_{N+T}(y^d)
\end{bmatrix} = \begin{bmatrix}
    U^d_{T} \\ U^d_{N} \\ Y^d_T \\ Y^d_N
\end{bmatrix},
$$

then we have 

$$ \begin{aligned}
\begin{bmatrix}
    \mathcal{H}_{N+T}(u^d) \\ \mathcal{H}_{N+T}(y^d)
\end{bmatrix} \alpha(t) = \begin{bmatrix}
    u_{\text{ini}}(t) \\ \hat{u}(t) \\ y_{\text{ini}}(t) \\ \hat{y}(t)
\end{bmatrix} & \Leftrightarrow \begin{bmatrix}
    U^d_{T}\alpha(t) \\ U^d_{N}\alpha(t) - \hat{u}(t) \\ Y^d_T\alpha(t) \\ Y^d_N\alpha(t)-\hat{y}(t) \end{bmatrix}= \begin{bmatrix} u_{\text{ini}}(t) \\ 0 \\ y_{\text{ini}}(t) \\ 0 \end{bmatrix} \\ & \Leftrightarrow 
    \begin{bmatrix}
        U_T^d & 0 & 0 \\
        U_N^d & -I & 0 \\
        Y_T^d & 0 & 0 \\
        Y_N^d & 0 & -I
    \end{bmatrix} \begin{bmatrix}
        \alpha(t) \\ \hat{u}(t) \\ \hat{y}(t) \\ 
    \end{bmatrix} = \begin{bmatrix} u_{\text{ini}}(t) \\ 0 \\ y_{\text{ini}}(t) \\ 0 \end{bmatrix}.
\end{aligned} $$

Introduce the box-constraints 

$$
y_{\min,k}(t) \le \hat{y}_k(t) \le y_{\max,k}(t), \quad k\in\mathbb{I}_{[0,N-1]}\\
u_{\min,k}(t) \le \hat{u}_k(t) \le u_{\max,k}(t), \quad k\in\mathbb{I}_{[0,N-1]}\\
$$

with the vectors $Y_{\max}(t) := \mathrm{col}_{k=0}^{N-1}(y_{\max,k}),\;Y_{\min}(t) := \mathrm{col}_{k=0}^{N-1}(y_{\min,k})$ and the similar as $U_{\min}(t),U_{\max}(t)$ into matrix form

$$
\underbrace{\begin{bmatrix} u_{\text{ini}}(t) \\ 0 \\ y_{\text{ini}}(t) \\ 0  \\ U_{\min}(t) \\ Y_{\min}(t) \end{bmatrix}}_{G_{lb}(t)} \le
\underbrace{\begin{bmatrix}
        U_T^d & 0 & 0 \\
        U_N^d & -I_{mN} & 0 \\
        Y_T^d & 0 & 0 \\
        Y_N^d & 0 & -I_{pN} \\
        0 & I_{mN} & 0 \\
        0 & 0 & I_{pN}
\end{bmatrix}}_{F} \underbrace{\begin{bmatrix}
        \alpha(t) \\ \hat{u}(t) \\ \hat{y}(t) \\ 
\end{bmatrix}}_{U(t)} \le \underbrace{\begin{bmatrix} u_{\text{ini}}(t) \\ 0 \\ y_{\text{ini}}(t) \\ 0  \\ U_{\max}(t) \\ Y_{\max}(t) \end{bmatrix}}_{G_{ub}(t)}.
$$

临时记：DeePC方法优化变量多，求解通常更耗时，存在对同一模型的冗余表示。

#### 2.2.2 Rewrite the Cost

Generally, the cost function is selected as

$$
\begin{aligned}
    J(t) &= \sum_{k=0}^{N-1}\ell \left( \hat{u}_k(t),\hat{y}_k(t) \right) \\
    &= \sum_{k=0}^{N-1} \left(\left\|\hat{y}_k(t) - y_k^r(t)\right\|^2_Q + \left\|\hat{u}_k(t) - u_k^r(t)\right\|^2_R \right) \\
    &= \left\| \hat{y}(t) - y^r(t) \right\|_{\mathcal{Q}}^2 + \left\| \hat{u}(t) - u^r(t) \right\|_{\mathcal{R}}^2 \\
    &= U^\top(t)\underbrace{\begin{bmatrix} 0 & & \\ & \mathcal{R} & \\ & & \mathcal{Q} \\ \end{bmatrix}}_{H} U(t) + 2 \underbrace{\begin{bmatrix} 0 & -(u^r(t))^\top\mathcal{R} &  -(y^r(t))^\top \mathcal{Q} \end{bmatrix}}_{f^\top(t)} U(t) + \\
    & \quad \; \underbrace{(u^r(t))^\top \mathcal{R} u^r(t)+(y^r(t))^\top \mathcal{Q} y^r(t)}_{c(t)} .
\end{aligned}
$$

where $u^r(t) = \mathrm{col}_{k=0}^{N-1}(u^r_k(t)),\;y^r(t) = \mathrm{col}_{k=0}^{N-1}(y^r_k(t))$ and $$\mathcal{R} = I_N \otimes R,\;\mathcal{Q} = I_N \otimes Q $$ with symmetric positive-definite weight matrix $Q\in\mathbb{R}^{p\times p},R\in\mathbb{R}^{m\times m}$; $\otimes$ is the kronecker-product.

## 3 Beyond Theory to Practice

### 3.1 Slack Variables

Due to the unavoidable noise in the actual system, the equation constraints cannot be always satisfied, leading to infeasibility of the optimization problem. Thus the slack variables are introduced to the optimization: 

$$
\begin{bmatrix}
    \mathcal{H}_{L}(u^d) \\ \mathcal{H}_{L}(y^d)
\end{bmatrix} \alpha(t) = \begin{bmatrix}
    u_{\text{ini}}(t) \\ \hat{u}(t) \\ y_{\text{ini}}(t) \\ \hat{y}(t)
\end{bmatrix} + \begin{bmatrix}\sigma_u(t) \\ 0 \\ \sigma_y(t) \\ 0\end{bmatrix}.
$$

Let

$$
\tilde{U}(t) = \begin{bmatrix} \alpha^\top(t) & \hat{u}^\top(t) & \hat{y}^\top(t) & \sigma_u^\top(t) & \sigma_y^\top(t) \end{bmatrix}^\top,
$$

then 

$$
G_{lb}(t) \le
\underbrace{\begin{bmatrix}
        U_T^d & 0 & 0 & -I_{mT} & 0 \\
        U_N^d & -I_{mN} & 0 & 0 & 0 \\
        Y_T^d & 0 & 0 & 0 & -I_{pT}\\
        Y_N^d & 0 & -I_{pN} & 0 & 0\\
        0 & I_{mN} & 0 & 0 & 0\\
        0 & 0 & I_{pN} & 0 & 0
\end{bmatrix}}_{\tilde{F}} \underbrace{\begin{bmatrix}
        \alpha(t) \\ \hat{u}(t) \\ \hat{y}(t) \\  \sigma_u(t) \\ \sigma_y(t)
\end{bmatrix}}_{\tilde{U}(t)} \le {G_{ub}(t)}.
$$

### 3.2 Regulation Terms

The slack variables should not be too large, or it can destroy the accuracy of data-based representation. Meanwhile, the decision variable $\alpha$ should also be limited into a reasonable range. Hence,

$$
\begin{aligned}
\tilde{J}(t) &= \left\| \hat{y}(t) - y^r(t) \right\|_{\mathcal{Q}}^2 + \left\| \hat{u}(t) - u^r(t) \right\|_{\mathcal{R}}^2 + \lambda_u \|\sigma_u(t)\|^2 + \lambda_y \|\sigma_y(t)\|^2 + \lambda_\alpha \|\alpha(t)\|^2 \\ 
& = \tilde{U}^\top(t)\underbrace{\mathrm{diag}(\Lambda_\alpha,\mathcal{R},\mathcal{Q},\Lambda_u,\Lambda_y)}_{\tilde{H}} \tilde{U}(t) + \\ 
& \quad \; 2 \underbrace{\left[ 0, -(u^r(t))^\top\mathcal{R}, -(y^r(t))^\top \mathcal{Q},0,0 \right]}_{\tilde{f}^\top(t)} \tilde{U}(t) + c(t).
\end{aligned}
$$

where $$\Lambda_\alpha := \lambda_\alpha I_{S-L+1},\quad \Lambda_u := \lambda_u I_{mT},\quad \Lambda_y := \lambda_y I_{pT},$$ with $\lambda_\alpha \ge 0$ being the regulation weight and $\lambda_u>0,\lambda_y>0$ being both sufficiently large constant.