---
globs: "*.m"
---

# MATLAB Rules

- Execute .m files using MATLAB MCP, never via shell
- Do not modify rng seeds — they ensure reproducible Monte Carlo results
- Monte Carlo simulations: N >= 80000 steps, steady_start at 50% of N
- Save all figures to the figures/ directory
- Physical constants must stay consistent across scripts:
  Ts = 1/1600, kb = 1.3806503e-23, T_temp = 310.15 K, gammaN = 0.0425 pN*s/um
- Use a_x = Ts/gammaN for mobility (um/pN)
- Variance units: um^2 (multiply by 1e18 from SI)
