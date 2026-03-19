# 3-State PP Controller: Complete Derivation and C_pp Correction

## 1. Plant Dynamics (受控體動態)

### 1.1 Position dynamics (位置動態)

Simulink model: probe in fluid, discrete-time:

    z[k+1] = z[k] + a_x * (fd[k] + fT[k])

- z[k]: probe position at step k [um]
- a_x = Ts / gammaN: motion gain [um/pN]
  - Ts = 1/1600 s (sampling period)
  - gammaN = 0.0425 pN*s/um (Stokes drag)
  - a_x = 0.01471 um/pN
- fd[k]: control force [pN] (designed)
- fT[k]: thermal noise force [pN] (random, sigma_fT = sqrt(4*kB*T*gammaN_SI/Ts) * 1e12)

### 1.2 Tracking error dynamics (追蹤誤差動態)

Define tracking error: e[k] = z_d - z[k], where z_d is the target position (constant).

Derivation:

    e[k+1] = z_d - z[k+1]                            (definition)
           = z_d - [z[k] + a_x*(fd[k] + fT[k])]     (substitute plant)
           = (z_d - z[k]) - a_x*(fd[k] + fT[k])      (rearrange)
           = e[k] - a_x*fd[k] - a_x*fT[k]            (use e[k] = z_d - z[k])

Physical meaning: next-step error = current error - control displacement - noise displacement.

### 1.3 Measurement delay (量測延遲)

The optical sensor measures position with a 2-step processing delay:

    p_measured[k] = z[k-2]    (position from 2 steps ago)

Measured tracking error:

    dzm[k] = z_d - p_measured[k]
           = z_d - z[k-2]
           = e[k-2]              (we observe the error from 2 steps ago)

### 1.4 The delay problem (延遲的問題)

Without delay, a simple controller fd = (1/a_x)*(1-lc)*e[k] gives:

    e[k+1] = lc*e[k] - a_x*fT[k]    (perfect AR(1), pole at lc)

But we only have e[k-2], not e[k]. Naive control fd = (1/a_x)*(1-lc)*e[k-2]
creates a 3rd-order system with poles that deviate from lc.

Two approaches to compensate this delay:
- **Eq. 17**: direct delay compensation using measurement + control history
- **3-state PP**: observer-based delay compensation using state estimation


## 2. Delay Chain State Space (延遲鏈狀態空間)

Define the state vector to capture the delay chain:

    x1[k] = e[k-2]    (= dzm[k], the measurement)
    x2[k] = e[k-1]    (one step ahead of measurement)
    x3[k] = e[k]      (current error, not directly measurable)

State dynamics:

    x1[k+1] = x2[k]                                  (pure delay, shift)
    x2[k+1] = x3[k]                                  (pure delay, shift)
    x3[k+1] = x3[k] - a_x*fd[k] - a_x*fT[k]       (plant dynamics)

Output (measurement):

    y[k] = x1[k] = dzm[k]

Matrix form:

    A_plant = [0 1 0;  0 0 1;  0 0 1]
    B_u = [0; 0; -a_x]     (control input)
    B_w = [0; 0; -a_x]     (noise input)
    C_y = [1 0 0]           (output = x1)


## 3. Three-State Estimator (三狀態估測器)

### 3.1 Estimator structure (估測器結構)

The estimator tracks x_hat = [dz1_hat, dz2_hat, dz3_hat], estimating [e[k-2], e[k-1], e[k]].

**Control law** (using known a_x):

    fd[k] = (1/a_x) * (1 - lc) * dz3_hat[k]

When substituted into the x3 prediction:

    x3_hat_pred = dz3_hat - a_x * fd
               = dz3_hat - (1-lc) * dz3_hat
               = lc * dz3_hat

This is the key: the prediction for x3 is lc * dz3_hat (closed-loop prediction).

### 3.2 Innovation and update (創新值與更新)

Innovation (prediction error):

    innov[k] = y[k] - dz1_hat[k] = dzm[k] - dz1_hat[k]

Combined predict + update step:

    dz1_hat[k+1] = dz2_hat[k]      + L1 * innov[k]    (predict: shift x2 -> x1)
    dz2_hat[k+1] = dz3_hat[k]      + L2 * innov[k]    (predict: shift x3 -> x2)
    dz3_hat[k+1] = lc*dz3_hat[k]   + L3 * innov[k]    (predict: lc*x3, closed-loop)

### 3.3 Observer gain design: pole placement (極點配置)

The estimation error dynamics (see Section 4) must have eigenvalues at z = le (triple).
Characteristic polynomial: (z - le)^3 = z^3 - 3*le*z^2 + 3*le^2*z - le^3

Matching coefficients gives:

    L1 = 1 - 3*le
    L2 = 1 - 3*le + 3*le^2
    L3 = (1 - le)^3

For le = 0.3: L1 = 0.1, L2 = 0.37, L3 = 0.343

### 3.4 Summary of one loop iteration (code)

```matlab
innov = dzm - dz1_hat;                          % innovation
fd_k  = (1/a_x) * (1 - lc) * dz3_hat;          % control
dz1_hat = dz2_hat     + L1 * innov;             % update
dz2_hat = dz3_hat     + L2 * innov;
dz3_hat = lc*dz3_hat  + L3 * innov;
z_new   = z + a_x * (fd_k + fT(k));            % plant
```


## 4. Estimation Error Dynamics (估測誤差動態)

### 4.1 Derivation (推導)

Define estimation errors: et_i[k] = x_i[k] - x_i_hat[k]

**et1** (error in e[k-2] estimate):

    et1[k+1] = x1[k+1] - x1_hat[k+1]
             = x2[k] - (dz2_hat[k] + L1*(x1[k] - dz1_hat[k]))
             = (x2[k] - dz2_hat[k]) - L1*(x1[k] - dz1_hat[k])
             = et2[k] - L1*et1[k]

**et2** (error in e[k-1] estimate):

    et2[k+1] = x2[k+1] - x2_hat[k+1]
             = x3[k] - (dz3_hat[k] + L2*et1[k])
             = et3[k] - L2*et1[k]

**et3** (error in e[k] estimate):

    x3[k+1]     = lc*x3[k] + (1-lc)*et3[k] - a_x*fT[k]   (see Section 5.1)
    x3_hat[k+1] = lc*dz3_hat[k] + L3*et1[k]

    et3[k+1] = lc*x3[k] + (1-lc)*et3[k] - a_x*fT[k] - lc*dz3_hat[k] - L3*et1[k]
             = lc*et3[k] + (1-lc)*et3[k] - L3*et1[k] - a_x*fT[k]
             = et3[k] - L3*et1[k] - a_x*fT[k]

### 4.2 Matrix form

    [et1]       [-L1  1  0] [et1]   [ 0  ]
    [et2]     = [-L2  0  1] [et2] + [ 0  ] * fT[k]
    [et3]k+1    [-L3  0  1] [et3]k  [-a_x]

    A_obs = [-L1  1  0;  -L2  0  1;  -L3  0  1]

Eigenvalues of A_obs = {le, le, le} (by construction via L1, L2, L3).

Key point: the observer error is driven DIRECTLY by thermal noise fT through B_obs = [0; 0; -a_x].
The observer error dynamics are independent of the true state x (separation principle).


## 5. Full Closed-Loop System (完整閉迴路系統)

### 5.1 How et3 couples into x3

The true plant dynamics with estimator-based control:

    x3[k+1] = x3[k] - a_x*fd[k] - a_x*fT[k]
            = x3[k] - a_x * (1/a_x)*(1-lc)*dz3_hat[k] - a_x*fT[k]
            = x3[k] - (1-lc)*dz3_hat[k] - a_x*fT[k]

Since dz3_hat = x3 - et3:

            = x3[k] - (1-lc)*(x3[k] - et3[k]) - a_x*fT[k]
            = lc*x3[k] + (1-lc)*et3[k] - a_x*fT[k]
                          ^^^^^^^^^^^^^^
                          observer noise injection!

The term (1-lc)*et3[k] is the estimation noise coupling into the plant.
Smaller lc (more aggressive control) amplifies this coupling.

### 5.2 Full 6-state system

State: [x1, x2, x3, et1, et2, et3]

           [ 0   1   0    0    0    0   ]       [ 0  ]
           [ 0   0   1    0    0    0   ]       [ 0  ]
    A_cl = [ 0   0   lc   0    0  1-lc ] , B = [-a_x]
           [ 0   0   0   -L1   1    0   ]       [ 0  ]
           [ 0   0   0   -L2   0    1   ]       [ 0  ]
           [ 0   0   0   -L3   0    1   ]       [-a_x]

    C_out = [1  0  0  0  0  0]   (output = dzm = x1)

Closed-loop eigenvalues = {0, 0, lc, le, le, le}  (separation principle)

### 5.3 Two noise paths

Thermal noise fT reaches dzm through TWO paths:

    Path 1 (direct):  fT --(-a_x)--> x3 --> x2 --> x1 = dzm
    Path 2 (indirect): fT --(-a_x)--> et3 --(1-lc)--> x3 --> x2 --> x1 = dzm

Eq. 17 has only Path 1. The 3-state PP has both paths.
Path 2 is the "observer noise" that increases dzm variance.


## 6. Transfer Function (轉移函數)

### 6.1 From fT to et3

From the observer error dynamics:

    ET3(z) = -a_x * (z^2 + L1*z + L2) / (z - le)^3 * FT(z)

### 6.2 From fT to x3

    (z - lc)*X3(z) = (1-lc)*ET3(z) - a_x*FT(z)

    X3(z) = -a_x * N(z) / [(z-le)^3 * (z-lc)] * FT(z)

where N(z) = (1-lc)*(z^2 + L1*z + L2) + (z-le)^3

### 6.3 From fT to dzm

    DZM(z) = z^{-2} * X3(z)

    H_pp(z) = DZM/FT = -a_x * N(z) / [z^2 * (z-le)^3 * (z-lc)]

Poles: {0, 0, le, le, le, lc}  (6 poles)

### 6.4 Comparison with Eq. 17

    H_eq17(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^2 * (z-lc)]

Poles: {0, 0, lc}  (3 poles only)

The 3-state PP has 3 extra poles at z = le that add variance.


## 7. C Factor via Lyapunov Equation (C 因子)

### 7.1 Definition

    Var(dzm) = sigma_fT^2 * a_x^2 * C

    C = sum_{k=0}^{inf} h[k]^2    (impulse response energy)

### 7.2 Lyapunov equation

C_pp = P(1,1) where P solves the discrete Lyapunov equation:

    P = A_cl * P * A_cl' + B_norm * B_norm'

with B_norm = B / a_x = [0; 0; -1; 0; 0; -1]

Since P(1,1) = P(2,2) = P(3,3) (delay chain), we reduce to the 4-state
subsystem [x3, et1, et2, et3]:

    A4 = [lc    0    0   1-lc;
          0    -L1   1    0;
          0    -L2   0    1;
          0    -L3   0    1]

    B4 = [-1; 0; 0; -1]

    C_pp = P4(1,1)  where  P4 = A4 * P4 * A4' + B4 * B4'

### 7.3 Decomposition

From the Lyapunov equation, P4(1,1) satisfies:

    C_pp*(1-lc^2) = (1-lc)^2*Var(et3) + 2*lc*(1-lc)*Cov(x3,et3) + 1

Therefore:

    C_pp = 1/(1-lc^2) + (1-lc)^2*Var(et3)/(1-lc^2) + 2*lc*(1-lc)*Cov(x3,et3)/(1-lc^2)
           ^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           AR(1) term   observer contribution (always positive)

### 7.4 Comparison: C_eq17 vs C_pp

For Eq. 17 (no observer):

    C_eq17 = 2 + 1/(1-lc^2)

Impulse response: h = [0, 0, 1, 1, 1, lc, lc^2, lc^3, ...]
- Three initial values of 1 (from delay compensation feedforward)
- Then geometric decay with rate lc

For 3-state PP:

    C_pp = C_eq17 + Delta_C(lc, le)

The observer adds Delta_C > 0 for all le > 0.

### 7.5 Fundamental observer cost

In the deadbeat limit (le = 0):

    Delta_C(lc, le=0) = 1    (exactly)

Even the fastest possible observer adds exactly 1 to the C factor.
This is the irreducible cost of using an estimator.

Physical reason: each fT impulse creates a 3-step transient of estimation error
(observer order = 3), which couples into the plant through (1-lc)*et3.

### 7.6 Numerical results (le = 0.3)

    lc    C_eq17    C_pp      Delta_C   Delta/C_eq17
    0.4   3.190     5.040     1.849     58.0%
    0.5   3.333     5.246     1.913     57.4%
    0.6   3.562     5.540     1.978     55.5%
    0.7   3.961     6.007     2.046     51.7%
    0.8   4.778     6.897     2.119     44.4%
    0.9   7.263     9.462     2.199     30.3%

### 7.7 MATLAB computation

C_pp can be computed exactly via `dlyap()`:

```matlab
le = 0.3;
L1 = 1-3*le;  L2 = 1-3*le+3*le^2;  L3 = (1-le)^3;
A4 = [lc 0 0 1-lc; 0 -L1 1 0; 0 -L2 0 1; 0 -L3 0 1];
B4 = [-1; 0; 0; -1];
P4 = dlyap(A4, B4*B4');
C_pp = P4(1,1);
```


## 8. Corrected Eq. 13 (修正的 Eq. 13)

### 8.1 Original Eq. 13 (for Eq. 17 controller)

    a_xm = sigma^2_dzr / (4*kB*T * C_eq17(lc))

where C_eq17(lc) = 2 + 1/(1-lc^2).

This is exact ONLY for the Eq. 17 controller.

### 8.2 Corrected Eq. 13 (for 3-state PP)

    a_xm = sigma^2_dzr / (4*kB*T * C_pp(lc, le))

where C_pp(lc, le) is computed via the Lyapunov equation.

### 8.3 Verification results

Using C_pp correction (le = 0.3):

    lc    axm (C_eq17)  err%    axm (C_pp)   err%
    0.4   0.02292       55.8    0.01451      1.3
    0.5   0.02335       58.8    0.01484      0.9
    0.6   0.02322       57.9    0.01493      1.5
    0.7   0.02260       53.6    0.01490      1.3
    0.8   0.02206       50.0    0.01528      3.9
    0.9   0.01912       30.0    0.01468      0.2

Error reduced from ~55% to < 4%.


## 9. Summary (總結)

1. Eq. 17 and 3-state PP both achieve the same plant pole lc.

2. But the transfer function from fT to dzm has different orders:
   - Eq. 17: 3 poles {0, 0, lc}
   - 3-state PP: 6 poles {0, 0, lc, le, le, le}

3. The extra observer poles add variance through the indirect noise path:
   fT -> et3 -> (1-lc)*et3 -> x3 -> dzm

4. C_pp(lc, le) > C_eq17(lc) for all le > 0. Even deadbeat (le=0) adds exactly 1.

5. Using C_pp in Eq. 13 correctly recovers a_x to < 4% error.

6. The C_pp formula is computed via the discrete Lyapunov equation of the
   4-state subsystem [x3, et1, et2, et3].
