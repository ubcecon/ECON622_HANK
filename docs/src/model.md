 # Model Description

 ## Households

```math
\begin{aligned}
 \max_{c_t, \ell_t, a_t, b_t, d_t} E_0 & \int_0^\infty
e^{−(ρ+ζ)t} u( c_t, \ell_t) dt \\
\text{s.t.} \\
\frac{db}{dt} & = (1 − τ_t)w_t z_t ℓ_t + r^b_t(b_t) b_t + T_t − d_t − \chi(d_t, a_t) − c_t, \\
\frac{da}{dt} & = r_t^a a_t + d_t \\
b_t & \geq −b_ ,\;\; a_t \geq 0
\end{aligned}
 ```
 where ``a_0``, ``b_0`` and ``d_0`` are given, and
 ```math
 r^b_t(b_t) = r_t^b + \kappa
 ```
 and
 ```\math
\chi(d_t, a_t) =  \chi_0 |d_t| - \chi_1 \left\vert \frac{d_t}{a_t}\right \vert^{\chi_2}
```


## Final Goods

```math
Y_t = \left( \int_0^1 y_{j,t}^{\frac{\epsilon-1}{\epsilon}} dj\right)^\frac{\epsilon}{\epsilon-1}
```

```math
y_{j,t}(p_{j,t}) = \left(\frac{p_{j,t}}{P_t}\right)^{-\epsilon}
```

```math
P_t = \left( \int_0^1 p_{j,t}^{1-\epsilon} dj \right)^{\frac{1}{1-\epsilon}}
```

## Intermediate Goods

```math
\max_{p_t} \int_0^\infty e^{\int_0^t r_s^a ds} \underbrace{\left[ \tilde{\Pi}(p_t) - \Theta\left(\frac{\dot{p}_t}{p_t} \right) \right]}_{\Pi_t} dt
```
where
```math
\Theta\left(\frac{\dot{p}_t}{p_t} \right) = \frac{\theta}{2}\left(\frac{\dot{p}_t}{p_t}\right)^2 Y_t
```
,
```math
\tilde{\Pi}(p_t) = \left( \frac{p_t}{P_t} - m_t \right) \left(\frac{p_t}{P_t} \right)^{-\epsilon} Y_t
```
and
```math
m_t = \left( \frac{r_t^k}{\alpha} \right)^\alpah \left( \frac{w_t}{1-\alpha} \right)^{1-\alpha}
```

## Illiquid Walth

```math
\dot{k}_t + q_t \dot{s}_t = (r_t^k - \delta) k_t + \Pi_t s_t + d_t
```

```math
\frac{\Pi_t + \dot{q}_t}{q_t} = r_t^k -\delta \equiv r_t^a
```

## Monetary Authority

```math
i_t = \bar{r}^b + \phi \pi_t \epsilon_t
```

```math
r_t^b = i_t - \pi_t
```

## Bond Market

## Government

```math
\dot{B}_t^g + G_t + T_t = \tau_t \int w_t z\ell_t(a, b, z) d\mu_t + r_t^b B_t^g
```

## Market Clearing

```math
\begin{aligned}
0 = & \int b d\mu_t + B_t^g \\
\int a d\mu_t = & K_t + q_t \\
N_t = & \int z \ell_t(a,b,z) d\mu_t \\
Y_t = & C_t + I_t + G_t + \Theta_t + \chi_t + \kappa \int\max\{-b, 0\} d\mu_t
\end{aligned}
```
