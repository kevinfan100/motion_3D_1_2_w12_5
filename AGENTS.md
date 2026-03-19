Use MATLAB MCP as the default tool for this workspace when the task involves MATLAB, Simulink, numerical verification, plotting, or simulation.

For this repository:
- Prefer MATLAB MCP for `.m` files, `.slx` models, verification scripts, figure generation, and toolbox inspection.
- When invoking MATLAB tools for repository work, use `C:\Users\PME406_01\Desktop\code\motion_3D_1_2_w12_5` as the project path / working folder.
- Fall back to shell or Python only for file management, Git operations, text editing, or when MATLAB MCP cannot perform the task.
- Do not ask the user to explicitly remind you to use MATLAB MCP for normal local MATLAB work in this repository.
- Respect Codex and Windows security boundaries; if a task cannot run under the current permissions, report that clearly instead of trying to bypass it.

## Project Structure

- **Root**: Simulink model (`.slx`) + driver scripts (`run_case1.m`, `path_control.m`)
- **verify/**: Standalone verification scripts (no Simulink dependency)
- **docs/**: Derivation documents (`.md`, `.tex`)
- **figures/**: All generated figure files (`.png`)
- **ref/**: Reference papers and presentations
- **agent_docs/**: Claude reference documents (architecture, equations, parameters, figure map)

## Key Conventions

- Figures are saved to `figures/` directory
- Verification scripts in `verify/` use `fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures')` for output paths
- Physical constants: `Ts=1/1600`, `kb=1.3806503e-23`, `T_temp=310.15`, `gammaN=0.0425`
- Do not modify rng seeds or the `.slx` model file
