% run_case1.m — Case 1: Known Gamma Verification
%
% Runs the Simulink model with known nominal gamma (mgain_z = Ts/gammaN)
% and verifies Eq.13 using three methods:
%   Method A: Direct Var(dz_k2) over steady-state (gold standard)
%   Method B: IIR filter azm_k / Am_scaling (Simulink internal)
%   Theory:   mgain_z = Ts/gammaN (nominal)
%
% Before running: ensure the Simulink model has Case 1 modifications
% (Trajectory fixed at z=25, Control law uses mgain_z, Avar=0.005).

clear; clc; close all;

%% ===== Configurable Parameters =====
lamdaC = 0.6;       % controller bandwidth (change to 0.9 for second run)

%% ===== Physical Constants =====
Ts      = 1/1600;
kb      = 1.3806503e-23;       % Boltzmann constant [J/K]
T_temp  = 273.15 + 37;         % temperature [K]
gammaN  = 0.0425;              % nominal Stokes drag [pN*s/um]
Am_scaling = 10;               % scaling factor used in Simulink

a_x_nominal = Ts / gammaN;     % theoretical motion gain [um/pN]
sigma2_nz   = (2.3e-3)^2;      % measurement noise variance [um^2]

%% ===== Setup Model =====
model = 'motion_3D_1_2_w12_5';
load_system(model);
open_system(model);
set_param(model, 'StopTime', '50');

% Programmatically set lamdaC in the Parameters block
rt = sfroot;
chart = rt.find('-isa', 'Stateflow.EMChart', 'Path', [model '/Parameters']);
script = chart.Script;
lines = strsplit(script, newline);
for i = 1:numel(lines)
    if contains(lines{i}, 'Parameters.lamdaC') && contains(lines{i}, '=')
        lines{i} = sprintf('    Parameters.lamdaC = %.2f;', lamdaC);
        break;
    end
end
script = strjoin(lines, newline);
chart.Script = script;
fprintf('Set lamdaC = %.2f in Parameters block.\n', lamdaC);

%% ===== Define ParametersBus (required by model) =====
clear elms;
elms(1) = Simulink.BusElement; elms(1).Name='Ts'; elms(1).DataType='double'; elms(1).Dimensions=1;
elms(2) = Simulink.BusElement; elms(2).Name='lamdaC'; elms(2).DataType='double'; elms(2).Dimensions=1;
elms(3) = Simulink.BusElement; elms(3).Name='theta'; elms(3).DataType='double'; elms(3).Dimensions=1;
elms(4) = Simulink.BusElement; elms(4).Name='phi'; elms(4).DataType='double'; elms(4).Dimensions=1;
elms(5) = Simulink.BusElement; elms(5).Name='pz'; elms(5).DataType='double'; elms(5).Dimensions=1;
elms(6) = Simulink.BusElement; elms(6).Name='R'; elms(6).DataType='double'; elms(6).Dimensions=1;
elms(7) = Simulink.BusElement; elms(7).Name='kb'; elms(7).DataType='double'; elms(7).Dimensions=1;
elms(8) = Simulink.BusElement; elms(8).Name='T'; elms(8).DataType='double'; elms(8).Dimensions=1;
elms(9) = Simulink.BusElement; elms(9).Name='ax_normal'; elms(9).DataType='double'; elms(9).Dimensions=1;
elms(10) = Simulink.BusElement; elms(10).Name='az_normal'; elms(10).DataType='double'; elms(10).Dimensions=1;
elms(11) = Simulink.BusElement; elms(11).Name='gammaN'; elms(11).DataType='double'; elms(11).Dimensions=1;
elms(12) = Simulink.BusElement; elms(12).Name='Avar'; elms(12).DataType='double'; elms(12).Dimensions=1;
elms(13) = Simulink.BusElement; elms(13).Name='Avar2'; elms(13).DataType='double'; elms(13).Dimensions=1;
elms(14) = Simulink.BusElement; elms(14).Name='Avar22'; elms(14).DataType='double'; elms(14).Dimensions=1;
elms(15) = Simulink.BusElement; elms(15).Name='Avar3'; elms(15).DataType='double'; elms(15).Dimensions=1;
elms(16) = Simulink.BusElement; elms(16).Name='Am_scaling'; elms(16).DataType='double'; elms(16).Dimensions=1;
elms(17) = Simulink.BusElement; elms(17).Name='beta'; elms(17).DataType='double'; elms(17).Dimensions=1;
elms(18) = Simulink.BusElement; elms(18).Name='lamdaF'; elms(18).DataType='double'; elms(18).Dimensions=1;
elms(19) = Simulink.BusElement; elms(19).Name='Pfz_11'; elms(19).DataType='double'; elms(19).Dimensions=1;
elms(20) = Simulink.BusElement; elms(20).Name='Pfz_22'; elms(20).DataType='double'; elms(20).Dimensions=1;
elms(21) = Simulink.BusElement; elms(21).Name='Pfz_33'; elms(21).DataType='double'; elms(21).Dimensions=1;
elms(22) = Simulink.BusElement; elms(22).Name='Pfz_44'; elms(22).DataType='double'; elms(22).Dimensions=1;
elms(23) = Simulink.BusElement; elms(23).Name='Pfz_55'; elms(23).DataType='double'; elms(23).Dimensions=1;
elms(24) = Simulink.BusElement; elms(24).Name='Pfz_66'; elms(24).DataType='double'; elms(24).Dimensions=1;
elms(25) = Simulink.BusElement; elms(25).Name='Pfz_77'; elms(25).DataType='double'; elms(25).Dimensions=1;
elms(26) = Simulink.BusElement; elms(26).Name='rz11_scaling'; elms(26).DataType='double'; elms(26).Dimensions=1;
elms(27) = Simulink.BusElement; elms(27).Name='rz22_scaling'; elms(27).DataType='double'; elms(27).Dimensions=1;
elms(28) = Simulink.BusElement; elms(28).Name='qz11_scaling'; elms(28).DataType='double'; elms(28).Dimensions=1;
elms(29) = Simulink.BusElement; elms(29).Name='qz22_scaling'; elms(29).DataType='double'; elms(29).Dimensions=1;
elms(30) = Simulink.BusElement; elms(30).Name='qz33_scaling'; elms(30).DataType='double'; elms(30).Dimensions=1;
elms(31) = Simulink.BusElement; elms(31).Name='qz44_scaling'; elms(31).DataType='double'; elms(31).Dimensions=1;
elms(32) = Simulink.BusElement; elms(32).Name='qz55_scaling'; elms(32).DataType='double'; elms(32).Dimensions=1;
elms(33) = Simulink.BusElement; elms(33).Name='qz66_scaling'; elms(33).DataType='double'; elms(33).Dimensions=1;
elms(34) = Simulink.BusElement; elms(34).Name='qz77_scaling'; elms(34).DataType='double'; elms(34).Dimensions=1;
ParametersBus = Simulink.Bus;
ParametersBus.Elements = elms;
assignin('base', 'ParametersBus', ParametersBus);

%% ===== Run Simulation =====
fprintf('Running simulation with lamdaC = %.2f (StopTime = 50s)...\n', lamdaC);
out = sim(model);
fprintf('Simulation complete.\n');

%% ===== Extract Data =====
N = length(azm_k);
t = (0:N-1)' * Ts;
ss = round(N/2):N;    % steady-state: last 25 seconds

%% ===== Eq.13 Verification =====
C_lc = 2 + 1 / (1 - lamdaC^2);
den_eq13 = 4 * kb * T_temp * C_lc;    % [N*m]

% --- Method A: Direct Var(dz_k2) ---
sigma2_direct = var(dz_k2(ss));        % [um^2]
% Eq.13 with measurement noise correction
sigma2_corrected = sigma2_direct - (2/(1 + lamdaC)) * sigma2_nz;
if sigma2_corrected < 0, sigma2_corrected = sigma2_direct; end
axm_direct = (sigma2_corrected * 1e-12) / den_eq13 * 1e-6;   % [um/pN]

% --- Method A (no noise correction) ---
axm_direct_raw = (sigma2_direct * 1e-12) / den_eq13 * 1e-6;

% --- Method B: IIR filter (azm_k / Am_scaling) ---
axm_iir = mean(azm_k(ss)) / Am_scaling;

% --- Method C: EKF az_hat_k ---
axm_ekf = mean(az_hat_k(ss));

% --- Theory ---
axm_theory = a_x_nominal;              % Ts/gammaN

% --- Also compute actual mgain_z from plant (mean of logged mgain_z) ---
axm_plant = mean(mgain_z(ss));         % Ts/gamma_z (includes c_perp)

%% ===== Eq.12 check: sigma^2 theory vs measured =====
sigma2_theory = C_lc * 4 * kb * T_temp * axm_plant * 1e18;  % [um^2], use plant a_x

%% ===== Print Results =====
fprintf('\n============================================================\n');
fprintf('  Case 1 Verification (lamdaC = %.2f)\n', lamdaC);
fprintf('  C(lc) = %.4f\n', C_lc);
fprintf('============================================================\n\n');

fprintf('--- Eq.12: Tracking error variance ---\n');
fprintf('  sigma^2 theory (plant a_x):  %.4e um^2\n', sigma2_theory);
fprintf('  sigma^2 measured Var(dz_k2): %.4e um^2\n', sigma2_direct);
fprintf('  sigma^2 corrected (-noise):  %.4e um^2\n', sigma2_corrected);
fprintf('  Ratio measured/theory:       %.3f\n\n', sigma2_direct/sigma2_theory);

fprintf('--- Eq.13: Recovered motion gain a_xm ---\n');
fprintf('  %-35s  %.6f um/pN\n', 'Theory (Ts/gammaN):', axm_theory);
fprintf('  %-35s  %.6f um/pN\n', 'Plant actual (Ts/gamma_z):', axm_plant);
fprintf('  %-35s  %.6f um/pN  (err: %.1f%%)\n', ...
    'Method A: Direct Var (raw):', axm_direct_raw, ...
    abs(axm_direct_raw - axm_plant)/axm_plant*100);
fprintf('  %-35s  %.6f um/pN  (err: %.1f%%)\n', ...
    'Method A: Direct Var (corrected):', axm_direct, ...
    abs(axm_direct - axm_plant)/axm_plant*100);
fprintf('  %-35s  %.6f um/pN  (err: %.1f%%)\n', ...
    'Method B: IIR (azm_k/Am_scaling):', axm_iir, ...
    abs(axm_iir - axm_plant)/axm_plant*100);
fprintf('  %-35s  %.6f um/pN  (err: %.1f%%)\n', ...
    'Method C: EKF az_hat_k:', axm_ekf, ...
    abs(axm_ekf - axm_plant)/axm_plant*100);

fprintf('\n--- Stability ---\n');
fprintf('  Simulation stable: %s\n', string(~any(isnan(dz_k2) | isinf(dz_k2))));
fprintf('  std(dz_k2) steady-state: %.4e um\n', std(dz_k2(ss)));

%% ===== Figure: Running variance vs time =====
% Compute running Var(dz_k2) [um^2]
run_var_ts = zeros(N, 1);
dz_sum = 0; dz2_sum = 0;
for k = 1:N
    if k > round(N/4)   % skip initial transient
        n_so_far = k - round(N/4);
        dz_sum  = dz_sum + dz_k2(k);
        dz2_sum = dz2_sum + dz_k2(k)^2;
        run_var_ts(k) = dz2_sum/n_so_far - (dz_sum/n_so_far)^2;
    end
end

fig1 = figure('Position', [50 50 1000 550], 'Color', 'w');
plot(t, run_var_ts, 'b', 'LineWidth', 1.5, ...
    'DisplayName', 'Running Var(dz_{k2})');
hold on;
yline(sigma2_theory, 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Theory \\sigma^2 = %.4e', sigma2_theory));
hold off;

xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\sigma^2_{dxr}  (\mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Case 1: Running Variance (\\lambda_c = %.1f)', lamdaC), ...
    'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;

lc_str = strrep(sprintf('%.1f', lamdaC), '.', '');
fname1 = sprintf('fig_case1_variance_lc%s.png', lc_str);
saveas(fig1, fullfile('figures', fname1));
fprintf('\nFigure saved: %s\n', fname1);
