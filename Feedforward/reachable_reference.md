# Tracking MPC for Reachale Reference

- Model: $$x_{t+1} = f(x_t,u_t),\;t\in\mathbb{N}.$$

- Input and State Constraints <u>(point-wise in time)</u>: $$(x_t,u_t)\in\mathcal{Z}.$$

- **Reachable** Reference: $r_t:=\mathrm{col}(x^r_t,u^r_t)$
    1. Dynamics: $$x^r_{t+1} = f(x^r_t,u^r_t),$$
    2. Constraints: $$r_t \in \mathcal{Z}_r \subseteq \mathrm{int}(\mathcal{Z}),$$

- Tracking MPC Scheme: $(\cdot)_{k|t}$ denotes the state/input at predictive time $t+k$. 
$$ \begin{aligned}
\min \quad & V_t = \sum_{k=0}^{N-1} \ell \left(x_{k|t},u_{k|t},r_{k|t}\right) + V_f\left( x_{N|t},r_{N|t} \right). \\
\mathrm{s.t.} \quad & x_{0|t} = x_t, \\
& x_{k+1|t} = f(x_{k|t},u_{k|t}), \\
& (x_{k|t},u_{k|t}) \in \mathcal{Z}, \quad k \in \mathbb{I}_{[0,N-1]} \\
& x_{N|t} \in \mathcal{X}_f(r_{N|t}). 
\end{aligned} $$ where $$ \ell(x_{k|t},u_{k|t},r_{k|t}) = \|x_{k|t}-x^r_{k|t}\|^2_Q + \|u_{k|t}-u^r_{k|t}\|^2_R $$ and $\mathcal{X}_f(r) \subset \mathbb{R}^n$ is the terminal set depend on corresponding reference, $V_f: \mathbb{R}^n \times \mathcal{Z}_r \to \mathbb{R}$ is terminal cost following certain condition.

- Conditions:
    - Exists feedback-feedforward controller: $$\kappa_f:\mathbb{R}^n \times \mathcal{Z}_r \to \mathbb{R}^m. $$
    - Satisfies constraints: $x^+ = f\left(x,\kappa_f(x,r)\right)$
        1. Terminal cost: $$V_f(x^+,r^+) - V_f(x,r) \le -\ell \left(x,k_f(x,r),r\right).$$
        2. Terminal set: $$(x,\kappa_f(x,r))\in \mathcal{Z}, \quad x^+ \in \mathcal{X}_f(r^+).$$
        3. Initial cost (ensured by terminal controller): $$V_t \le c_u \left\|x_t-x^r_t\right\|_Q^2.$$

- Results:
    - Feasible control sequences:
        $$u_{k|t+1} = u^*_{k+1|t},\;k\in\mathbb{I}_{[0,N-2]}$$ $$u_{N-1|t+1} = \kappa_f \left(x^*_{N|t},r_{N|t}\right).$$
    - Cost reducing:
        $$ V_{t+1} - V_t \le - \ell(x_{t},u_{t}) \le -\left\|x_t-x^r_t\right\|_Q^2,$$ and $\left\|x_t-x^r_t\right\|_Q^2 \le V_t \le c_v \left\|x_t-x^r_t\right\|_Q^2$ for some $c_v \ge c_u \ge 1.$

