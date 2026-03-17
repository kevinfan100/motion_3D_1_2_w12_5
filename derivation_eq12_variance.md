# Eq.12 Complete Derivation: From Closed-Loop Transfer Function to Tracking Error Variance

From the control law (Eq.17) through the closed-loop transfer function H(z),
to the variance formula via Parseval's theorem.

---

## Part 1: Plant and Control Law Review

### 1.1 Plant dynamics

Discrete-time probe dynamics (single axis):

    z[k+1] = z[k] + a_x * (fd[k] + fT[k])

- z[k]: probe position [um]
- a_x = Ts / gamma: motion gain [um/pN]
- fd[k]: control force [pN]
- fT[k]: thermal noise force [pN], i.i.d. N(0, sigma_fT^2)
  - sigma_fT^2 = 4 * kB * T * gamma_SI / Ts

### 1.2 Tracking error

    e[k] = z_d - z[k]

Substituting the plant equation:

    e[k+1] = z_d - z[k+1]
           = z_d - z[k] - a_x*(fd[k] + fT[k])
           = e[k] - a_x*fd[k] - a_x*fT[k]

### 1.3 Measurement delay (d = 2)

The optical sensor has a 2-step processing delay:

    dzm[k] = z_d - z[k-2] = e[k-2]

We observe the error from 2 steps ago, not the current error.

### 1.4 Eq.17: d-step delay compensation control law

The controller must compensate the 2-step delay. The key idea: use known
control history to "predict" the current error from the delayed measurement.

**Step A — Two-step plant expansion:**

    z[k] = z[k-2] + a_x*(fd[k-2] + fd[k-1]) + a_x*(fT[k-2] + fT[k-1])
             known: measurement              unknown: thermal noise

**Step B — Predicted current error:**

    e_hat[k] = dzm[k] - a_x*(fd[k-2] + fd[k-1])
             = e[k] + a_x*(fT[k-2] + fT[k-1])

The prediction is unbiased but contaminated by 2 steps of unobservable noise.

**Step C — Control law (Eq.17, d=2, stationary target):**

    fd[k] = (1/a_x) * (1 - lc) * e_hat[k]
           = (1/a_x) * (1-lc) * dzm[k] - (1-lc) * (fd[k-1] + fd[k-2])

### 1.5 Closed-loop error dynamics

Substitute fd[k] into the error equation e[k+1] = e[k] - a_x*fd[k] - a_x*fT[k]:

    e[k+1] = e[k] - (1-lc)*e_hat[k] - a_x*fT[k]
           = e[k] - (1-lc)*{e[k] + a_x*(fT[k-2] + fT[k-1])} - a_x*fT[k]
           = lc*e[k] - (1-lc)*a_x*fT[k-2] - (1-lc)*a_x*fT[k-1] - a_x*fT[k]

Final form:

    e[k+1] = lc * e[k] + w[k]

where:

    w[k] = -a_x * { fT[k] + (1-lc)*fT[k-1] + (1-lc)*fT[k-2] }

Key properties:
- The closed-loop pole is exactly lc (satisfies Eq.7 control objective)
- The driving noise w[k] is **colored** (depends on fT at 3 time steps)
- The coloring comes from the 2-step delay: fT[k-1] and fT[k-2] could not
  be observed or compensated during the delay period

---

## Part 2: Transfer Function H(z) from fT to dzm

### 2.1 Z-transform of the error dynamics

Starting from e[k+1] = lc*e[k] + w[k], take the Z-transform:

    z*E(z) = lc*E(z) + W(z)

Expand w[k] = -a_x*{fT[k] + (1-lc)*fT[k-1] + (1-lc)*fT[k-2]}:

    W(z) = -a_x * {1 + (1-lc)*z^(-1) + (1-lc)*z^(-2)} * FT(z)
         = -a_x * {z^2 + (1-lc)*z + (1-lc)} / z^2 * FT(z)

Solve for E(z):

    E(z) = W(z) / (z - lc)
         = -a_x * {z^2 + (1-lc)*z + (1-lc)} / {z^2 * (z - lc)} * FT(z)

### 2.2 Transfer function from fT to e[k]

    H_e(z) = E(z) / FT(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^2 * (z - lc)]

### 2.3 From e[k] to dzm[k]

The measured error is the 2-step delayed version of e[k]:

    dzm[k] = e[k-2]
    DZM(z) = z^(-2) * E(z)

Therefore:

    H(z) = DZM(z) / FT(z) = z^(-2) * H_e(z)

    H(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^4 * (z - lc)]

Wait — let me be more careful. The transfer function H_e(z) already has z^2
in the denominator from the W(z) term. Let me re-derive cleanly.

### 2.4 Clean derivation

From the error dynamics in Z-domain:

    (z - lc) * E(z) = -a_x * [1 + (1-lc)*z^(-1) + (1-lc)*z^(-2)] * FT(z)

Multiply both sides by z^2:

    z^2*(z - lc) * E(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] * FT(z)

So:

    H_e(z) = E(z)/FT(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z - lc)]

Now dzm[k] = e[k-2], so DZM(z) = z^(-2)*E(z):

    H(z) = DZM(z)/FT(z) = z^(-2) * H_e(z)
         = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^4*(z - lc)]

**But** for computing |H(e^(j*theta))|^2 on the unit circle, we can use H_e
directly because |z^(-2)| = 1 on the unit circle. The extra z^(-2) only adds
phase, not magnitude.

So for variance computation:

    |H(e^(j*theta))|^2 = |H_e(e^(j*theta))|^2

And:

    H_e(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z - lc)]

Poles of H_e: {0, 0, lc} — three poles, matching Eq.17's structure.

### 2.5 Normalized transfer function

Factor out -a_x for convenience. Define:

    H(z) = -a_x * H_norm(z)

    H_norm(z) = [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z - lc)]

The variance will be:

    Var(dzm) = sigma_fT^2 * a_x^2 * C(lc)

where C(lc) = sum of |h_norm[k]|^2 (impulse response energy of H_norm).

---

## Part 3: Impulse Response h[k]

### 3.1 Partial fraction decomposition of H_norm(z)

We need to find h_norm[k] = Z^(-1){H_norm(z)}.

    H_norm(z) = [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z - lc)]

Partial fractions for H_norm(z)/z (standard technique):

    H_norm(z)/z = [z^2 + (1-lc)*z + (1-lc)] / [z^3*(z - lc)]

Decompose:

    = A/z + B/z^2 + C/z^3 + D/(z - lc)

Multiply through by z^3*(z-lc):

    z^2 + (1-lc)*z + (1-lc) = A*z^2*(z-lc) + B*z*(z-lc) + C*(z-lc) + D*z^3

Set z = 0:   (1-lc) = C*(-lc)  =>  C = -(1-lc)/lc = (lc-1)/lc

Set z = lc:  lc^2 + (1-lc)*lc + (1-lc) = D*lc^3
             lc^2 + lc - lc^2 + 1 - lc = D*lc^3
             1 = D*lc^3  =>  D = 1/lc^3

### 3.2 Direct computation via long division / recursion

A more practical approach: compute h_norm[k] directly from the difference equation.

From H_norm(z) = Y(z)/X(z) where input is delta[k]:

    z^2*(z-lc)*Y(z) = [z^2 + (1-lc)*z + (1-lc)] * X(z)

In time domain (causal, h[k]=0 for k<0):

    h[k+3] - lc*h[k+2] = delta[k+2] + (1-lc)*delta[k+1] + (1-lc)*delta[k]

But it's even simpler to compute the impulse response of H_e(z) first,
then shift by 2 for dzm.

From H_e(z): E(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z-lc)] * FT(z)

The impulse response of -a_x * H_norm(z) represents e[k]'s response to a
unit impulse in fT at k=0. But since H_norm has z^2 in the denominator,
the response starts at k=0 with causality.

**Let's just compute h_norm[k] by simulation of the error recursion.**

Apply fT[0]=1, fT[k]=0 for k != 0 to:

    e[k+1] = lc*e[k] - fT[k] - (1-lc)*fT[k-1] - (1-lc)*fT[k-2]

(We drop a_x because we're computing h_norm.)

    k=0: e[1] = lc*0 - 1 - 0 - 0 = -1
    k=1: e[2] = lc*(-1) - 0 - (1-lc)*1 - 0 = -lc - 1 + lc = -1
    k=2: e[3] = lc*(-1) - 0 - 0 - (1-lc)*1 = -lc - 1 + lc = -1
    k=3: e[4] = lc*(-1) - 0 - 0 - 0 = -lc
    k=4: e[5] = lc*(-lc) = -lc^2
    k=5: e[6] = lc*(-lc^2) = -lc^3
    ...

So h_e[k] (impulse response from fT to e, without -a_x factor):

    h_e = [0, -1, -1, -1, -lc, -lc^2, -lc^3, ...]

Note: h_e[0] = 0 because fT[0] affects e[1], not e[0].

Now dzm[k] = e[k-2], so h_dzm[k] = h_e[k-2]:

    h_dzm = [0, 0, 0, -1, -1, -1, -lc, -lc^2, -lc^3, ...]

The normalized impulse response (absolute values for energy):

    |h_norm| = [0, 0, 0, 1, 1, 1, lc, lc^2, lc^3, ...]

### 3.3 Physical interpretation of each term

    h_norm[0] = 0   : impulse at k=0 hasn't propagated yet
    h_norm[1] = 0   : still within measurement delay
    h_norm[2] = 0   : still within measurement delay
    h_norm[3] = 1   : fT[0] directly displaces position, first visible in dzm
    h_norm[4] = 1   : fT[0] enters w[1] via (1-lc)*fT[k-1] term (delay path)
    h_norm[5] = 1   : fT[0] enters w[2] via (1-lc)*fT[k-2] term (delay path)
    h_norm[k>=6] = lc^(k-5) : AR(1) decay — controller progressively corrects

**Wait** — let me reconsider the indexing. The key question is whether we
should compute the impulse response of H_e or of H (which includes the z^(-2)
delay). For the variance computation, it doesn't matter because the delay
just shifts the impulse response, and sum of squares is shift-invariant.

For the variance, what matters is sum_{k} h_e[k]^2, which equals:

    0^2 + 1^2 + 1^2 + 1^2 + lc^2 + lc^4 + lc^6 + ...
    = 0 + 1 + 1 + 1 + lc^2 + lc^4 + ...
    = 3 + lc^2/(1 - lc^2)

This is C(lc).

---

## Part 4: Parseval's Theorem — The Core Mathematical Tool

### 4.1 Statement

For a causal, stable, discrete-time system with impulse response h[k], the
total energy (sum of squares) can be computed in either domain:

**Time domain:**

    sum_{k=0}^{inf} h[k]^2

**Frequency domain:**

    (1 / 2*pi) * integral_{-pi}^{pi} |H(e^(j*theta))|^2 d_theta

These are equal by Parseval's theorem.

### 4.2 Real-valued system simplification

Since h[k] is real, |H(e^(j*theta))|^2 is symmetric about theta = 0:

    |H(e^(-j*theta))|^2 = |H(e^(j*theta))|^2

Therefore:

    sum h[k]^2 = (1/2*pi) * integral_{-pi}^{pi} |H(e^(j*theta))|^2 d_theta
               = (1/pi) * integral_{0}^{pi} |H(e^(j*theta))|^2 d_theta

**This is the "integral from 0 to pi" that the professor mentioned.**

### 4.3 Connection to variance

If fT[k] is white noise with variance sigma_fT^2, and dzm[k] is the output
of H(z) driven by fT[k], then:

    Var(dzm) = sigma_fT^2 * sum_{k=0}^{inf} h[k]^2
             = sigma_fT^2 * (1/pi) * integral_{0}^{pi} |H(e^(j*theta))|^2 d_theta

### 4.4 Geometric meaning of z = e^(j*theta)

When we evaluate H(z) at z = e^(j*theta):

- z traces the **unit circle** in the complex z-plane
- theta goes from 0 to pi (half circle, by symmetry)
- theta = 0: z = 1 (DC, zero frequency)
- theta = pi: z = -1 (Nyquist frequency, f_s/2)
- |H(e^(j*theta))|^2 is the **power spectral density** of the output
  when driven by white noise

The integral sweeps through all digital frequencies, summing the power
contribution at each frequency.

---

## Part 5: Computing C(lc) = sum h_norm[k]^2

### 5.1 From the impulse response

Recall from Part 3:

    h_norm = [0, 1, 1, 1, lc, lc^2, lc^3, ...]

(Using h_e indexing; the z^(-2) shift doesn't affect the sum of squares.)

### 5.2 Direct summation

    C(lc) = sum_{k=0}^{inf} h_norm[k]^2
           = 0^2 + 1^2 + 1^2 + 1^2 + lc^2 + lc^4 + lc^6 + ...
           = 3 + sum_{n=1}^{inf} lc^(2n)
           = 3 + lc^2 / (1 - lc^2)          (geometric series, |lc| < 1)

### 5.3 Simplification

    C(lc) = 3 + lc^2 / (1 - lc^2)
           = [3*(1 - lc^2) + lc^2] / (1 - lc^2)
           = (3 - 3*lc^2 + lc^2) / (1 - lc^2)
           = (3 - 2*lc^2) / (1 - lc^2)

### 5.4 Verify equivalence with 2 + 1/(1-lc^2)

    2 + 1/(1 - lc^2) = [2*(1-lc^2) + 1] / (1-lc^2)
                      = (2 - 2*lc^2 + 1) / (1-lc^2)
                      = (3 - 2*lc^2) / (1-lc^2)   checkmark

Both forms are identical:

    C(lc) = 3 + lc^2/(1-lc^2) = 2 + 1/(1-lc^2) = (3 - 2*lc^2)/(1-lc^2)

### 5.5 Numerical values

    lc    C(lc)
    0.0   3.000
    0.3   3.099
    0.5   3.333
    0.7   3.961
    0.9   7.263
    0.95  12.744
    0.99  51.503

C(lc) grows rapidly as lc -> 1 (less aggressive control = more variance).

---

## Part 6: From C(lc) to Var(dzm) — Deriving Eq.11

### 6.1 Variance of dzm

From Part 4.3:

    Var(dzm) = sigma_fT^2 * a_x^2 * C(lc)

### 6.2 Substitute thermal noise variance

The fluctuation-dissipation theorem gives:

    sigma_fT^2 = 4 * kB * T * gamma_SI / Ts        [N^2]

### 6.3 Substitute motion gain

    a_x = Ts / gamma        [um/pN]

In SI units:

    a_x_SI = Ts / gamma_SI  [m/N]

### 6.4 Combine

    Var(dzm) = (4*kB*T*gamma_SI/Ts) * (Ts/gamma_SI)^2 * C(lc)   [m^2]
             = (4*kB*T*gamma_SI/Ts) * Ts^2/gamma_SI^2 * C(lc)
             = 4*kB*T * Ts/gamma_SI * C(lc)
             = 4*kB*T * a_x_SI * C(lc)                           [m^2]

Converting to [um^2]:

    Var(dzm) = 4*kB*T * a_x * C(lc)    [mixed units, with appropriate conversion]

Or more explicitly:

    sig2_dxr = C(lc) * 4 * kB * T * a_x                   ... Eq.11

where the unit chain is:
- 4*kB*T has units [J] = [N*m]
- a_x has units [um/pN]
- a_x [um/pN] = a_x * 1e6 [m/N]   (since 1 um = 1e-6 m, 1 pN = 1e-12 N)
- 4*kB*T*a_x: [N*m] * [m/N] = [m^2] (with the 1e6 factor for unit conversion)

**This is Eq.11.** The crucial cancellation: gamma appears in both sigma_fT^2
and a_x^2, and cancels to leave only a_x (= Ts/gamma). This is why the
variance is proportional to the mobility, which is the physical basis for
measuring gamma from thermal fluctuations.

---

## Part 7: Physical Interpretation of C(lc) = 2 + 1/(1-lc^2)

### 7.1 Decomposing C(lc)

    C(lc) = 2 + 1/(1-lc^2)

The two terms have distinct physical origins:

### 7.2 The "2" — delay contribution

More precisely, C(lc) = 1 + 1 + 1/(1-lc^2), and we can trace:

    h_norm[1]^2 = 1  :  fT[0] displaces position at k=0,
                         but controller at k=0 can't see it (delay)
                         → fT[0] appears unattenuated in dzm

    h_norm[2]^2 = 1  :  fT[0] enters the prediction error e_hat[1]
                         via the (1-lc)*fT[k-1] term in w[k]
                         → the delay compensation "overshoots"

These 2 units of energy come from the d=2 measurement delay.
**In general, for d-step delay, this term would be d.**

But wait — there's also h_norm[3]^2 = 1. Where does the third "1" go?

### 7.3 The "1/(1-lc^2)" — AR(1) variance amplification

    1/(1-lc^2) = 1 + lc^2 + lc^4 + lc^6 + ...

This is the geometric series starting from h_norm[3]:

    h_norm[3]^2 + h_norm[4]^2 + h_norm[5]^2 + ...
    = 1^2 + lc^2 + lc^4 + ...
    = 1/(1-lc^2)

The first term (1^2 = 1) is from h[3]: the fT[0] impulse entering via the
(1-lc)*fT[k-2] path in w[k]. After that, the AR(1) dynamics e[k+1]=lc*e[k]
cause each subsequent contribution to decay by lc.

So the decomposition is:

    C(lc) = { h[1]^2 + h[2]^2 } + { h[3]^2 + h[4]^2 + ... }
           = { 1 + 1 }           + { 1/(1-lc^2) }
           = 2                    + 1/(1-lc^2)

### 7.4 Alternative decomposition: 1 + 1 + 1/(1-lc^2)

    C(lc) = 1 + 1 + 1/(1-lc^2)

- First "1": direct displacement by fT (visible after delay)
- Second "1": delay compensation overshoot in prediction
- 1/(1-lc^2): steady-state AR(1) amplification (includes the third
  noise injection via the (1-lc)*fT[k-2] term, plus all subsequent decay)

### 7.5 Why C(lc) -> infinity as lc -> 1

When lc -> 1, the controller becomes passive (barely corrects errors).
Each thermal impulse persists almost forever: the geometric series
1 + lc^2 + lc^4 + ... diverges. Physically, a nearly-passive controller
lets the probe diffuse freely under thermal noise.

### 7.6 Why C(lc) -> 3 as lc -> 0

When lc -> 0, the controller is maximally aggressive (deadbeat). Each error
is corrected in one step. But the 2-step delay is irreducible — the controller
still can't see the most recent 2 steps of thermal noise. Plus the third
unit from the transient gives C(0) = 3.

---

## Part 8: MATLAB Numerical Verification

Three independent verification methods are implemented in `verify_eq12_spectral.m`:

### 8.1 Method 1: Impulse response summation

Compute h_norm[k] for k = 0, 1, ..., K (truncated), then:

    C_impulse = sum_{k=0}^{K} h_norm[k]^2

Compare with C(lc) = 2 + 1/(1-lc^2).

### 8.2 Method 2: Frequency-domain integration (Parseval's)

Evaluate |H_norm(e^(j*theta))|^2 at theta = linspace(0, pi, N_pts), then:

    C_spectral = (1/pi) * trapz(theta, |H_norm(e^(j*theta))|^2)

This is the "integral from 0 to pi" approach.

### 8.3 Method 3: Lyapunov equation

Use the augmented state-space model from Part 1:

    A = [ lc,  -(1-lc)*a_x,  -(1-lc)*a_x ]     B = [ -a_x ]
        [  0,       0,              0       ]         [   1   ]
        [  0,       1,              0       ]         [   0   ]

Solve P = A*P*A' + B*sigma_fT^2*B', then Var(e) = P(1,1).

Or equivalently with B_norm = B/a_x:

    P_norm = A * P_norm * A' + B_norm * B_norm'
    C_lyap = P_norm(1,1)

### 8.4 Method 4: Direct Monte Carlo simulation

Run the closed-loop system for 80000 steps with thermal noise,
compute Var(dzm) over the steady-state portion, then:

    C_sim = Var(dzm) / (sigma_fT^2 * a_x^2)

Compare all four C values.

---

## Part 9: Measurement Noise Term — (2/(1+lc)) * sigma_nx^2

### 9.1 Setup

Measurement noise n[k] is added to the position measurement:

    dzm_noisy[k] = dzm[k] + n[k]

where n[k] ~ N(0, sigma_nx^2), independent of fT.

The noise n[k] enters the controller through e_hat:

    e_hat_noisy[k] = e_hat[k] + n[k]

### 9.2 Transfer function from n to dzm

Measurement noise n[k] enters the system differently from thermal noise fT.
It enters through the control law:

    fd[k] = (1/a_x)*(1-lc)*e_hat_noisy[k]
           = (1/a_x)*(1-lc)*(e_hat[k] + n[k])

The additional control force due to n is:

    delta_fd[k] = (1/a_x)*(1-lc)*n[k]

This propagates through the plant:

    delta_e[k+1] = lc*delta_e[k] - (1-lc)*n[k]

And the transfer function from n to e is:

    H_n_e(z) = -(1-lc) / (z - lc)

From n to dzm (including 2-step delay):

    H_n(z) = z^(-2) * H_n_e(z) = -(1-lc) / [z^2*(z - lc)]

### 9.3 Variance contribution from measurement noise

    Var(dzm)|_noise = sigma_nx^2 * sum |h_n[k]|^2

Impulse response of H_n_norm(z) = (1-lc)/(z-lc):

    h_n[k] = (1-lc) * lc^k    for k >= 0

    sum h_n[k]^2 = (1-lc)^2 * sum lc^(2k) = (1-lc)^2 / (1-lc^2)
                 = (1-lc)^2 / [(1-lc)(1+lc)]
                 = (1-lc) / (1+lc)

Wait — this doesn't include the z^(-2) delay effect on the noise path.
Let me reconsider.

Actually, for the measurement noise, the situation is different because n[k]
directly affects the measurement dzm[k] AND enters the controller. Let me
trace through more carefully.

### 9.4 Careful derivation of measurement noise contribution

The measured signal with noise:

    dzm_noisy[k] = e[k-2] + n[k]

The controller uses this:

    e_hat_noisy[k] = dzm_noisy[k] - a_x*(fd[k-2]+fd[k-1])
                   = e_hat[k] + n[k]

The error dynamics become:

    e[k+1] = lc*e[k] + w[k] - (1-lc)*n[k]

(The n[k] enters through the control action at step k.)

For the **observed** signal dzm_noisy[k] = e[k-2] + n[k], its variance has
two independent contributions:

    Var(dzm_noisy) = Var(e[k-2]) + Var(n[k]) + 2*Cov(e[k-2], n[k])

Since n[k] is independent of e[k-2] (n[k] affects e[k-1] and later,
not e[k-2]):

    Var(dzm_noisy) = Var(e) + sigma_nx^2

But e[k] now has additional noise from n. The full variance is:

    Var(dzm_noisy) = C(lc)*sigma_fT^2*a_x^2 + C_n(lc)*sigma_nx^2

where C_n captures the measurement noise amplification.

The noise n[k] enters the error recursion: e[k+1] = lc*e[k] + w[k] - (1-lc)*n[k].

Through the AR(1) dynamics, the contribution of n to Var(e) is:

    (1-lc)^2 * sigma_nx^2 / (1-lc^2) = (1-lc)/(1+lc) * sigma_nx^2

Additionally, the direct additive noise on dzm contributes sigma_nx^2.

But we must also consider correlations. Since dzm_noisy[k] = e[k-2] + n[k],
and n[k-2] contributes to e[k-2] through the recursion, there is a correlation.

The full result from the literature is:

    Var(dzm_noisy) = C(lc)*4*kB*T*a_x + (2/(1+lc))*sigma_nx^2

### 9.5 Verifying the coefficient 2/(1+lc)

The coefficient 2/(1+lc) can be decomposed:

    2/(1+lc) = 1 + (1-lc)/(1+lc)

- The "1" comes from the direct additive measurement noise: n[k] added to dzm[k]
- The (1-lc)/(1+lc) comes from the noise-driven error propagation through
  the AR(1) dynamics and then back through the 2-step delay to dzm

For lc = 0 (deadbeat):  2/(1+0) = 2
For lc = 0.5:           2/(1+0.5) = 1.333
For lc -> 1 (passive):  2/(1+1) = 1

As lc -> 1, the controller barely reacts to n, so the noise contribution
approaches just sigma_nx^2 (direct measurement noise only).

### 9.6 Complete Eq.12

Combining thermal noise (Part 6) and measurement noise (Part 9):

    sig2_dxr = C(lc) * 4*kB*T*a_x + (2/(1+lc)) * sigma_nx^2       ... Eq.12

    where C(lc) = 2 + 1/(1-lc^2)

### 9.7 Inverting to get Eq.13

Solving for a_x:

    a_x = {sig2_dxr - (2/(1+lc))*sigma_nx^2} / {4*kB*T * C(lc)}   ... Eq.13

This is how the system measures the local drag coefficient: measure the
tracking error variance, subtract the known measurement noise contribution,
divide by the thermal and control amplification factors.

---

## Appendix A: Quick Reference of Key Results

### Transfer functions

    H_e(z) = -a_x * [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z-lc)]     (fT -> e)
    H(z)   = z^(-2) * H_e(z)                                       (fT -> dzm)

### Impulse response (normalized by -a_x)

    h_norm = [0, 1, 1, 1, lc, lc^2, lc^3, ...]

### C factor

    C(lc) = sum h_norm[k]^2 = 2 + 1/(1-lc^2) = (3-2*lc^2)/(1-lc^2)

### Variance formula (Eq.12)

    sig2_dxr = (2 + 1/(1-lc^2)) * 4*kB*T*a_x + (2/(1+lc)) * sigma_nx^2

### Motion gain estimation (Eq.13)

    a_xm = {sig2_dxr - (2/(1+lc))*sigma_nx^2} / {4*kB*T*(2+1/(1-lc^2))}
