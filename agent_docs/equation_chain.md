# Equation Chain: Microprobe Motion Control

## Eq.1 — Brownian Dynamics

    dz[k+1] = dz[k] + a_x * (fd[k] + fT[k])

- dz[k]: position increment at step k
- a_x = Ts / gamma: mobility (um/pN)
- fd[k]: deterministic control force
- fT[k]: thermal (Brownian) force, white noise with variance sigma_fT^2

## Eq.7 — AR(1) Closed-Loop Error

With proportional control fd[k] = -(1 - lc) / a_x * dz[k]:

    e[k+1] = lc * e[k] - a_x * fT[k]

- lc (lambda_c): closed-loop pole, 0 < lc < 1
- e[k] = dz[k] - dz_ref[k]: tracking error

The error is a stable AR(1) process with bandwidth set by lc.

## Eq.9/10 — IIR Filter Decomposition

Separate dz into deterministic and stochastic components via IIR filtering:

- Deterministic: low-pass IIR extracts the trajectory-following part
- Stochastic residual: dz_r[k] = dz[k] - dz_det[k]

The stochastic residual variance is used for mobility estimation.

## Eq.11 — Variance via Impulse Response

    sigma^2_dxr = sum(h[k]^2, k=0..inf) * sigma_fT^2 * a_x^2

where h[k] is the impulse response of the stochastic transfer function.

## Eq.12 — Variance with C-Factor

    sigma^2_dxr = C(lc) * 4 * kb * T * a_x

where:

    C(lc) = 2 + 1 / (1 - lc^2)

This is the key result connecting steady-state variance to mobility.

## Eq.13 — Mobility Recovery

Rearranging Eq.12 to recover mobility from measured variance:

    a_xm = sigma^2_dxr / (4 * kb * T * C(lc))

This is the estimator equation: measure variance, divide by known constants to get a_x.

## Eq.17 — d-Step Delay Control Law

When the controller has a d-step delay (observer-based):

    fd[k] = (1/a_x) * (1 - lc) * dzm - (1 - lc) * sum(fd[k-i], i=1..d)

- dzm: observed/estimated position increment
- The sum accounts for previously applied forces during the delay

This modifies the closed-loop dynamics and changes the variance expression.

## C_obs Extension — Observer-Augmented Variance

With a 3-state pole-placement (PP) observer, the closed-loop system becomes a 4-state augmented system. The variance factor becomes:

    C_obs(lc, le) > C(lc)

where le (lambda_e) is the observer pole.

### Observer Gains (3-state PP)

    L1 = 1 - 3*le
    L2 = 1 - 3*le + 3*le^2
    L3 = (1 - le)^3

### Computation

C_obs is obtained by solving the Lyapunov equation on the 4-state augmented system (controller + observer states) and extracting the (1,1) element of the steady-state covariance.

### Deadbeat Property

When le = 0 (deadbeat observer):

    DeltaC = C_obs - C_eq17 = 1   (exact)

This means the observer adds exactly 1 unit to the C-factor compared to the Eq.17 prediction, regardless of lc. This has been verified both analytically and numerically.
