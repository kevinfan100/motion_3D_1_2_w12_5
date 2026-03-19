# Tracking Error Variance -- Observer (3-State PP)

Plain-text companion to `derivation_eq12_variance_PP.tex`.

---

## 1. Observer Gains (Open-Loop Design, A(3,3)=1)

Place all three eigenvalues of the error dynamics at le:

    L1 = 1 - 3*le
    L2 = 1 - 3*le + 3*le^2
    L3 = (1 - le)^3

Gains are independent of lc (open-loop design).

Deadbeat (le=0): L1 = L2 = L3 = 1.

---

## 2. 4-State Augmented System

State: s = [dz, e_dz1, e_dz2, e_dz3]'

    s[k+1] = A4 * s[k] + B4 * f_T[k]

    A4 = [ lc,    0,   0,  1-lc ]      B4 = [ -a_z ]
         [  0,  -L1,   1,   0   ]           [   0  ]
         [  0,  -L2,   0,   1   ]           [   0  ]
         [  0,  -L3,   0,   1   ]           [ -a_z ]

Eigenvalues of A4: {lc, le, le, le}.

---

## 3. Transfer Function H(z)

H(z) = C_out * (zI - A4)^(-1) * B4,  C_out = [1,0,0,0].

z-form:

    H(z) = -a_z * N(z) / [(z - lc)(z - le)^3]

    N(z) = z^3 + (L1 - lc)*z^2 + (lc*L1 - L2)*z + (lc*L2 - L3)

z^(-1)-form:

    H(z^(-1)) = -a_z * [z^(-1) + (L1-lc)*z^(-2) + (lc*L1-L2)*z^(-3) + (lc*L2-L3)*z^(-4)]
                / [(1 - lc*z^(-1))(1 - le*z^(-1))^3]

---

## 4. Variance Formula

    sigma^2_dz = C_obs * 4*kB*T*a_z

C_obs from Lyapunov equation:

    P = A4 * P * A4' + Bn * Bn'       (Bn = B4/a_z = [-1; 0; 0; -1])
    C_obs = P(1,1)

C_obs from Parseval:

    C_obs = (1/pi) * integral_0^pi |N(e^(j*theta))|^2 / [|e^(j*theta) - lc|^2 * |e^(j*theta) - le|^6] d(theta)

---

## 5. Deadbeat Closed Form (le = 0)

With le=0: L1=L2=L3=1.

    C_obs(lc, 0) = (4 - 3*lc^2) / (1 - lc^2) = 3 + 1/(1 - lc^2)

Compared with Eq.17 (no observer):

    C_eq17(lc) = (3 - 2*lc^2) / (1 - lc^2) = 2 + 1/(1 - lc^2)

    DeltaC = C_obs(lc,0) - C_eq17(lc) = 1       (exact, for all lc)

Physical meaning: the observer adds exactly one step of pipeline delay,
contributing exactly 1 to the variance coefficient regardless of lc.

---

## 6. Execution Order Comparison

Prediction form (control executes first):   C_pred
Filtering form  (correct executes first):   C_filt < C_pred

Deadbeat filtering result: C_filt(lc,0) = C_eq17(lc).
In the deadbeat filtering case the observer cost vanishes entirely.

---

## 7. Modified Eq.13

When using the PP observer, replace C_eq17 with C_obs:

    a_z^(m) = sigma^2_dz / (4 * kB * T * C_obs)

---

## 8. Numerical Table (le = 0.3)

    lc      C_pred    C_filt    C_eq17
    0.50     5.25      4.52      3.33
    0.60     5.54      4.85      3.56
    0.70     6.01      5.36      3.96
    0.80     6.90      6.33      4.78
    0.90     9.46      9.00      7.26
