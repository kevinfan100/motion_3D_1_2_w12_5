# Parameters Reference

## Physical Constants

| Parameter | Value | Units | Description |
|---|---|---|---|
| Ts | 1/1600 = 6.25e-4 | s | Sampling time |
| kb | 1.3806503e-23 | J/K | Boltzmann constant |
| T_temp | 310.15 | K | Temperature (37 C body temp) |
| R | 2.25e-6 | m | Probe radius |
| eta | 0.001 | Pa*s | Fluid viscosity (water at 37 C) |
| gammaN | 0.0425 | pN*s/um | Nominal Stokes drag = 6*pi*eta*R |
| a_x | Ts / gammaN | um/pN | Mobility (position increment per unit force per step) |
| sigma_fT | sqrt(4*kb*T*gamma_SI/Ts) * 1e12 | pN | Thermal force standard deviation |

## Control Parameters

| Parameter | Typical Range | Description |
|---|---|---|
| lc (lambda_c) | 0.1 - 0.95 | Controller bandwidth / closed-loop pole |
| le (lambda_e) | 0.3 (typical) | Observer pole |
| L1 | 1 - 3*le | Observer gain 1 |
| L2 | 1 - 3*le + 3*le^2 | Observer gain 2 |
| L3 | (1 - le)^3 | Observer gain 3 |

## EKF Parameters

| Parameter | Value | Description |
|---|---|---|
| beta | 0.5 | EKF forgetting factor |
| lamdaF | 1.0 | EKF lambda_F scaling |
| Pfz_11..77 | (diagonal entries) | EKF initial error covariance |
| qz11..qz77 | (scaling values) | Process noise Q diagonal scaling |
| rz11, rz22 | (scaling values) | Measurement noise R diagonal scaling |

## IIR Filter Parameters

| Parameter | Value | Context |
|---|---|---|
| Avar | 0.05 or 0.45 | IIR smoothing coefficient (stationary vs moving) |
| Avar22 | 0.05 | Secondary IIR coefficient |
| Avar3 | 0.05 | Tertiary IIR coefficient |
| Am_scaling | 10 | Amplitude scaling factor |

## Unit System

The project uses two unit systems in parallel:

- **SI units**: m, N, s, K (used in physical constants)
- **Working units**: um, pN (used in simulation and control)

### Key Conversions

| Conversion | Factor |
|---|---|
| 1 um | 1e-6 m |
| 1 pN | 1e-12 N |
| gamma: pN*s/um to N*s/m | multiply by 1e-6 |
| a_x [um/pN] to [m/N] | multiply by 1e6 |

### Variance Scaling

The thermal variance in working units:

    4 * kb * T * a_x * 1e18   [um^2]

The 1e18 factor comes from converting m^2 to um^2 (1e12) and N to pN (1e12), combined with the a_x unit conversion. This factor appears frequently in verification scripts when comparing analytic predictions to simulation results.
