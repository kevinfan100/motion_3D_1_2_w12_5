# Case 1: Eq.13 Verification — Fixed Trajectory, Known Gamma

## Experiment Setup

- **Model**: `motion_3D_1_2_w12_5.slx`
- **Branch**: `case1-known-gamma`
- **Date**: 2026-03-17

### Trajectory
- Fixed point: x_d=0, y_d=0, z_d=25 um (stationary)

### Wall Geometry
- pz (wall position) = -15 um
- R (probe radius) = 2.25 um
- h = z_d - pz = 40 um (probe center to wall)
- h - R = 37.75 um (probe surface to wall)
- hbar = h/R = 17.78
- c_perp(hbar) = 1.0675

### Control Law
- Uses `mgain_z = Ts/gamma_z` (known plant motion gain, not EKF estimate)
- fd_Z = (1/mgain_z) * (zd_k - zd_k1 + (1-lamdaC)*dz3_hat_k - zD1_hat_k1)
- EKF estimator still runs (provides dz3_hat_k, zD1_hat_k1)

### Parameters
- Ts = 1/1600 s
- gammaN = 0.0425 pN*s/um
- Avar = 0.005, Avar22 = 0.05, Avar3 = 0.05
- Am_scaling = 1
- StopTime = 50 s
- Steady-state window: 25~50 s

### EKF Input
- e_z1_k = dz_k2 - dz1_hat_k1
- dz_k2 = zd[k-2] - pz[k] (2-step delayed tracking error)

## Eq.13 Formula

```
C(lc) = 2 + 1/(1 - lc^2)

Eq.12:  sigma^2_dz = C(lc) * 4*kb*T * a_x     [um^2]
Eq.13:  a_xm = sigma^2_dz / (4*kb*T * C(lc))   [um/pN]

Unit conversion:
  a_xm [um/pN] = (Var(dz_k2) [um^2] * 1e-12) / (4*kb*T*C(lc) [N*m]) * 1e-6

With noise correction:
  sigma^2_corrected = Var(dz_k2) - (2/(1+lc)) * sigma^2_nz
  sigma^2_nz = (2.3e-3)^2 um^2
```

## Results (Steady-state 25~50 s)

### Reference Values
- Ts/gammaN (nominal)   = 0.014706 um/pN
- Ts/gamma_z (plant)    = 0.013777 um/pN  (with c_perp=1.0675)

### Multi-lambda_c Sweep

| lc  | C(lc)  | Var(dz_k2) [um^2] | a_xm (Eq.13) | vs Ts/gammaN | vs Ts/gamma_z |
|-----|--------|--------------------|---------------|-------------|--------------|
| 0.6 | 3.562  | 9.171e-04          | 0.014920      | +1.5%       | +8.3%        |
| 0.7 | 3.961  | 1.017e-03          | 0.014903      | +1.3%       | +8.2%        |
| 0.8 | 4.778  | 1.224e-03          | 0.014883      | +1.2%       | +8.0%        |
| 0.9 | 7.263  | 1.808e-03          | 0.014492      | -1.5%       | +5.2%        |

### Key Observations

1. **Eq.13 recovers Ts/gammaN within +/-1.5%** across all tested lc values.
   The blue data points cluster around the nominal Ts/gammaN line.

2. **Eq.13 does NOT recover Ts/gamma_z** (the c_perp-corrected value).
   At hbar=17.8, c_perp=1.0675 causes a 6.7% gap between Ts/gammaN and Ts/gamma_z.
   The 5-8% positive error vs Ts/gamma_z is almost entirely this c_perp effect.

3. **a_xm is approximately constant across lc** (~0.0149 for lc=0.6-0.8),
   confirming C(lc) correctly compensates the lc-dependent variance.
   lc=0.9 is slightly lower (0.01449), possibly due to finite-sample variance at high lc.

4. **Stability**: lc <= 0.5 diverges with the current EKF+mgain_z control law.
   Only lc = 0.6, 0.7, 0.8, 0.9 were stable for 50s.

## Figures

- `fig_eq13_multi_lc_v2.png` — a_xm vs lc with Ts/gammaN reference line
- `fig_eq13_verify_fixed_traj.png` — Time series: Eq.13/IIR/EKF vs theory (lc=0.9)
- `fig_eq13_verify_dzk2.png` — Tracking error dz_k2 time series (lc=0.9)
