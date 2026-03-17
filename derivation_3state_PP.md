# 3-State Pole-Placement Controller — Complete Derivation

## Part 1: Plant and Problem Setup

**Plant dynamics** (discrete time, sampling time Ts):

```
p[k+1] = p[k] + a_x * (fd[k] + fT[k])
```

- `p[k]`  : probe position [um]
- `a_x`   : motion gain = Ts/gamma [um/pN]
- `fd[k]` : control force [pN]
- `fT[k]` : thermal noise force [pN], iid N(0, sigma_fT^2)
  - sigma_fT^2 = 4*kB*T*gamma / Ts  (Fluctuation-Dissipation theorem)

**Tracking error** (stationary trajectory z_d = const):

```
e[k] = z_d - p[k]

e[k+1] = z_d - p[k+1]
       = e[k] - a_x*fd[k] - a_x*fT[k]        ... (1)
```

**Measurement delay** d = 2:

```
At step k, only p[k-2] is available.
Measurable signal: y[k] = z_d - p[k-2] = e[k-2]     ... (2)
```

**Control objective**:

```
e[k+1] = lc * e[k] + noise                    ... (3)
0 < lc < 1, smaller lc = faster decay
```

---

## Part 2: Control Law Derivation

From (1) and (3):

```
lc*e[k] = e[k] - a_x*fd[k]
a_x*fd[k] = (1 - lc)*e[k]
```

**Ideal control law**:
```
fd[k] = (1/a_x) * (1 - lc) * e[k]
```

Problem: e[k] cannot be measured directly (delay d=2).
Solution: use observer estimate dz3_hat[k].

**Actual control law**:
```
fd[k] = (1/a_x) * (1 - lc) * dz3_hat[k]
```

---

## Part 3: 3-State Observer Design

### State definition

```
x[k] = [dz1[k], dz2[k], dz3[k]]^T

dz1[k] ~ e[k-2]   (measurable delayed error)
dz2[k] ~ e[k-1]   (1-step-ago error)
dz3[k] ~ e[k]     (current error, used for control)
```

### State transition model

Closed-loop assumption (control law fd = (1/a_x)*(1-lc)*dz3 already applied):

```
e[k+1] = lc*e[k] + noise
```

Therefore:
```
dz1[k+1] = e[k-1] = dz2[k]      (time shift)
dz2[k+1] = e[k]   = dz3[k]      (time shift)
dz3[k+1] = e[k+1] = lc*dz3[k]   (closed-loop dynamics)
```

**State matrix**:
```
     [0  1   0]
A =  [0  0   1]       ... (4)
     [0  0  lc]
```

### Measurement model

```
y[k] = e[k-2] = dz1[k]
C = [1  0  0]                ... (5)
```

---

## Part 4: Luenberger Observer Update Law

**Structure** (combined predict + correct):

```
x_hat[k] = A * x_hat[k-1] + L * (y[k] - C * x_hat[k-1])
         = A * x_hat[k-1] + L * innovation[k]

innovation[k] = y[k] - dz1_hat[k-1]
```

**Expanded**:
```
dz1_hat[k] = dz2_hat[k-1] + L1 * innov[k]
dz2_hat[k] = dz3_hat[k-1] + L2 * innov[k]
dz3_hat[k] = lc*dz3_hat[k-1] + L3 * innov[k]

innov[k] = y[k] - dz1_hat[k-1]
```

Physical meaning:
- `A*x_hat`: model-based prediction
  - dz1 <- dz2 (shift), dz2 <- dz3 (shift), dz3 <- lc*dz3 (decay)
- `L*innov`: measurement-based correction
  - innov > 0 means actual error is larger than predicted, correct upward

---

## Part 5: Observer Gain via Pole Placement

### Estimation error dynamics

```
eps[k] = x[k] - x_hat[k]

eps[k] = (A - L*C) * eps[k-1] + w[k-1]
```

Convergence rate determined by eigenvalues of (A - LC).

### Goal

Place all three eigenvalues of (A - LC) at le (observer pole):

```
det(zI - A + LC) = (z - le)^3 = z^3 - 3*le*z^2 + 3*le^2*z - le^3
```

### Compute det(zI - A + LC)

```
A - LC = [-L1  1   0 ]
         [-L2  0   1 ]
         [-L3  0   lc]

det(zI - (A-LC)) = det([z+L1  -1    0   ])
                       [L2     z   -1   ]
                       [L3     0   z-lc ]

Expand (along third column):
= (z+L1)*z*(z-lc) + L2*(z-lc) + L3
= z^3 + (L1-lc)*z^2 + (-L1*lc+L2)*z + (-L2*lc+L3)
```

### Match coefficients

```
z^2:  L1 - lc       = -3*le     =>  L1 = lc - 3*le
z^1:  -L1*lc + L2   = 3*le^2    =>  L2 = lc^2 - 3*lc*le + 3*le^2
z^0:  -L2*lc + L3   = -le^3     =>  L3 = (lc - le)^3
```

### Observer gains (correct formulas)

```
L1 = lc - 3*le
L2 = lc^2 - 3*lc*le + 3*le^2
L3 = (lc - le)^3

Condition: le < lc (observer must be faster than controller)
```

**Note**: `verify_eq13_estimator.m` originally used `L1 = 1 - 3*le`, etc.,
which assumes A(3,3)=1 (open-loop). This does NOT match the closed-loop
prediction model A(3,3)=lc, causing the observer poles to deviate from
their intended positions.

### Numerical examples (le = 0.3)

| lc  | L1     | L2    | L3     | Observer poles |
|-----|--------|-------|--------|----------------|
| 0.5 | -0.400 | 0.070 | 0.008  | 0.3, 0.3, 0.3 |
| 0.7 | -0.200 | 0.130 | 0.064  | 0.3, 0.3, 0.3 |
| 0.9 | +0.000 | 0.270 | 0.216  | 0.3, 0.3, 0.3 |

---

## Part 6: Closed-Loop Error Analysis

### Actual closed-loop with estimation error

Define estimation error: `eps[k] = e[k] - dz3_hat[k]`

```
fd[k] = (1/a_x)*(1-lc)*dz3_hat[k]
      = (1/a_x)*(1-lc)*(e[k] - eps[k])
```

Substitute into Plant (1):

```
e[k+1] = e[k] - (1-lc)*(e[k] - eps[k]) - a_x*fT[k]

e[k+1] = lc*e[k] + (1-lc)*eps[k] - a_x*fT[k]
          ^^^^^^^   ^^^^^^^^^^^^^^   ^^^^^^^^^^^
          desired    observer bias    thermal noise
```

### Steady-state variance

Assuming e, eps, fT approximately independent:

```
Var(e) = lc^2*Var(e) + (1-lc)^2*Var(eps) + a_x^2*sigma_fT^2

Var(e)*(1-lc^2) = (1-lc)^2*Var(eps) + a_x^2*sigma_fT^2

Var(e) = a_x^2*sigma_fT^2/(1-lc^2) + ((1-lc)/(1+lc))*Var(eps)
```

Including the full C(lc) = 2 + 1/(1-lc^2) factor from d=2 delay:

```
Var(e) = C(lc)*4kT*a_x + ((1-lc)/(1+lc))*Var(eps)
         ^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^
         Eq.12 theory     observer extra term (>=0)
```

### Impact on Eq.13

```
Eq.13: a_xm = Var(e) / (4kT*C(lc))

a_xm = a_x + bias

bias = ((1-lc)/(1+lc)) * Var(eps) / (4kT*C(lc))

bias >= 0  =>  a_xm >= a_x  (always overestimates)
eps = 0    =>  bias = 0      (Eq.17 case, exact)
```

---

## Part 7: Complete Algorithm (per-step execution)

```matlab
% Initialization
dz1_hat = 0;  dz2_hat = 0;  dz3_hat = 0;
L1 = lc - 3*le;
L2 = lc^2 - 3*lc*le + 3*le^2;
L3 = (lc - le)^3;

% Each step k:

% (1) Measurement
y_k = z_d - p(k-2);              % 2-step delayed measurement

% (2) Innovation
innov = y_k - dz1_hat;           % measurement - prediction

% (3) Control force
fd_k = (1/a_x)*(1-lc)*dz3_hat;  % use current estimate

% (4) Observer update
dz1_new = dz2_hat     + L1*innov;
dz2_new = dz3_hat     + L2*innov;
dz3_new = lc*dz3_hat  + L3*innov;
dz1_hat = dz1_new;
dz2_hat = dz2_new;
dz3_hat = dz3_new;

% (5) Plant update
p(k+1) = p(k) + a_x*(fd_k + fT(k));
```

---

## Part 8: Extension to Moving Trajectory

For moving trajectory z_d[k] != const:

```
Plant: e[k+1] = (z_d[k+1]-z_d[k]) + e[k] - a_x*fd[k] + noise

Ideal:   fd[k] = (1/a_x)*[(z_d[k+1]-z_d[k]) + (1-lc)*e[k]]
                            ^^^^^^^^^^^^^^^^^ requires future value!

Approx:  fd[k] = (1/a_x)*[(z_d[k]-z_d[k-1]) + (1-lc)*dz3_hat[k]]
                            ^^^^^^^^^^^^^^^^^ use past value

Approximation error = z_d[k+1] - 2*z_d[k] + z_d[k-1]
                    = trajectory second-order difference (acceleration)

Stationary: 0     (exact)
Constant velocity: 0     (exact)
Accelerating: != 0 (small tracking bias)
```

With disturbance compensation (7-state EKF version):
```
fd[k] = (1/a_x)*[(z_d[k]-z_d[k-1]) + (1-lc)*dz3_hat - zD1_hat]
```

---

## Part 9: Comparison of Three Controllers

| Property        | Eq.17          | 3-state PP     | 7-state EKF    |
|-----------------|----------------|----------------|----------------|
| Observer        | Not needed     | 3-state fixed  | 7-state adaptive |
| eps             | 0 (exact)      | Large          | Smaller        |
| Disturbance     | None           | None           | Yes (zD1_hat)  |
| a_x known?      | Yes            | Yes            | No (estimates) |
| Eq.13 accuracy  | <1%            | 30-56%         | 5-8%           |
| Moving traj     | Modify Eq.17   | Add feedforward| Built-in       |
| Complexity      | Low            | Medium         | High           |

### Key conclusions

1. **Control law is correct**:
   `fd = (1/a_x)*[(zd[k]-zd[k-1]) + (1-lc)*dz3_hat - D_hat]`

2. **Any observer introduces eps > 0**:
   `e[k+1] = lc*e[k] + (1-lc)*eps[k] + noise`
   => `Var(e) > C(lc)*4kT*a_x`
   => Eq.13 always overestimates a_xm

3. **EKF reduces eps** (adaptive Kalman gain):
   `Var(eps)_EKF / Var(eps)_3PP ~ 0.55` (about half)
   But eps still != 0, so Eq.13 still has bias

4. **Only Eq.17 (no observer) can exactly verify Eq.13**
   Because eps = 0, closed-loop is pure AR(1)
