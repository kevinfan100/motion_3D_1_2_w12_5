# Project Architecture

## Directory Structure

```
motion_3D_1_2_w12_5/
  run_case1.m            Simulink driver script (sets params, runs model, post-processes)
  path_control.m         Path/trajectory control logic
  motion_3D_1_2_w12_5.slx  Simulink model (closed-loop microprobe control)
  startup.m              MATLAB startup configuration
  verify/
    verify_eq12_spectral.m       C(lc) 4-way check (analytic, sum(h^2), PSD, Monte Carlo)
    verify_eq12_spectral_PP.m    C_obs with 3-state observer (prediction vs filtering)
    verify_eq13.m                Basic Eq.12/13 check (sigma^2 vs lc)
    verify_eq13_estimator.m      Estimator variance vs Eq.17 theory
    verify_eq13_unified.m        3-method comparison: Eq.17 / PP observer / EKF
  docs/                  Derivation documents
  figures/               Output figures (.png)
  ref/                   Reference papers
  agent_docs/            Claude reference documents (this directory)
```

## Data Flow

1. **run_case1.m** sets up physical and control parameters, constructs `ParametersBus`.
2. Simulink model `motion_3D_1_2_w12_5.slx` runs the closed-loop simulation.
3. Model exports workspace variables:
   - `dz_k2` — position increments (2-step)
   - `azm_k` — measured acceleration
   - `az_hat_k` — estimated acceleration
   - `mgain_z` — mobility gain estimate
4. **run_case1.m** post-processes workspace data: computes running variance, plots figures.
5. Standalone verification scripts in `verify/` run independently of Simulink (Monte Carlo + analytic).

## ParametersBus

34-element `Simulink.Bus` defined in `run_case1.m`. Fields:

| Field | Description |
|---|---|
| Ts | Sampling time (1/1600 s) |
| lamdaC | Controller bandwidth lambda_c |
| theta, phi | Trajectory angles |
| pz | Path parameter |
| R | Probe radius |
| kb | Boltzmann constant |
| T | Temperature |
| ax_normal | Normal-direction mobility |
| az_normal | z-direction mobility |
| gammaN | Nominal Stokes drag |
| Avar | IIR filter coefficient (stationary) |
| Avar2 | IIR filter coefficient 2 |
| Avar22 | IIR filter coefficient 22 |
| Avar3 | IIR filter coefficient 3 |
| Am_scaling | Amplitude scaling factor |
| beta | EKF beta parameter |
| lamdaF | EKF lambda_F |
| Pfz_11..77 | EKF initial covariance diagonal (7 entries) |
| rz, qz scalings | EKF Q/R scaling parameters |

## Script Roles

| Script | Purpose |
|---|---|
| `run_case1.m` | Simulink driver: configure parameters, run model, post-process results |
| `verify/verify_eq12_spectral.m` | C(lc) 4-way verification: analytic formula, sum(h^2), PSD integration, Monte Carlo |
| `verify/verify_eq12_spectral_PP.m` | C_obs with 3-state observer: prediction-form vs filtering-form vs Eq.17 |
| `verify/verify_eq13.m` | Basic Eq.12/13 check: sigma^2 vs lc with theory, direct, and IIR curves |
| `verify/verify_eq13_estimator.m` | Estimator variance vs Eq.17 analytic prediction |
| `verify/verify_eq13_unified.m` | Unified 3-method comparison: Eq.17 analytic, PP observer, EKF |
