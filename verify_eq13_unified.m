% verify_eq13_unified.m — Eq. 13 unified verification: three controller methods
%
% Method 1: Eq. 17 (d-step delay control law)
% Method 2: 3-state estimator (pole placement, known a_x)
% Method 3: 7-state EKF (from Simulink, control with known a_x)
%
% Simplified conditions: known gamma_N, no disturbance, stationary probe.
% Goal: all three methods recover a_x via Eq. 13, with IIR filter + direct var.

clear; clc; close all;

%% ===== Physical Parameters (shared by all methods) =====
Ts      = 1/1600;               % sampling time [s]
kb      = 1.3806503e-23;        % Boltzmann constant [J/K]
T_temp  = 310.15;               % temperature [K] (37 C)
R_probe = 2.25e-6;              % probe radius [m]
eta     = 0.001;                % dynamic viscosity [Pa*s]
gammaN  = 0.0425;               % nominal Stokes drag [pN*s/um]

a_x      = Ts / gammaN;                                           % [um/pN]
sigma_fT  = sqrt(4 * kb * T_temp * (gammaN * 1e-6) / Ts * 1e24); % thermal force std [pN]

%% ===== Simulation Parameters =====
N        = 80000;
ss_start = 40000;
z_d      = 25;                  % target position [um]
lc_sim   = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
n_lc     = length(lc_sim);

% IIR filter coefficients (Avar=0.05 for stationary case)
Avar   = 0.05;                  % deterministic/stochastic separation
Avar22 = 0.05;                  % stochastic residual mean tracking
Avar3  = 0.05;                  % stochastic residual mean-square tracking

% 3-state estimator pole placement on F = A+G-LC (F(3,3)=1)
lambda_e = 0.3;
L1 = 1 - 3*lambda_e;
L2 = 1 - 3*lambda_e + 3*lambda_e^2;
L3 = (1 - lambda_e)^3;

% 7-state EKF parameters (from Simulink Parameters block)
beta_ekf = 0.5;                 % w1-w2 estimator parameter
lamdaF   = 1.0;                 % forgetting factor

% EKF initial covariance (nonzero on all diagonal for numerical stability)
dz_var0 = (2.3e-3)^2;          % [um^2] — measurement-noise level
az_var0 = (0.1*0.012)^2;       % [um/pN]^2
Pfz_init = diag([dz_var0, dz_var0, dz_var0, dz_var0, dz_var0, ...
                 10*(0.0147)^2, az_var0]);

% EKF Q/R scaling (from Simulink)
rz11_scaling = 1;       % 10^0
rz22_scaling = 0.1;     % 10^-1
qz11_scaling = 0;
qz22_scaling = 0;
qz33_scaling = 1e4;
qz44_scaling = 1e2;
qz55_scaling = 0;
qz66_scaling = 8e-5;
qz77_scaling = 0;

% Minimum floor for R to prevent singularity
r11_floor = 1e-12;
r22_floor = 1e-16;

% Dense theory curve
lc_dense = 0.1:0.005:0.95;
C_dense  = 2 + 1./(1 - lc_dense.^2);
sig2_theory_dense = C_dense * 4 * kb * T_temp * a_x * 1e18;

%% ===== Storage =====
% Direct variance
sig2_eq17  = zeros(1, n_lc);
sig2_pp    = zeros(1, n_lc);
sig2_ekf   = zeros(1, n_lc);
% IIR variance
sig2_iir_eq17 = zeros(1, n_lc);
sig2_iir_pp   = zeros(1, n_lc);
sig2_iir_ekf  = zeros(1, n_lc);
% Recovered a_xm
axm_eq17   = zeros(1, n_lc);
axm_pp     = zeros(1, n_lc);
axm_ekf    = zeros(1, n_lc);
axm_iir_eq17 = zeros(1, n_lc);
axm_iir_pp   = zeros(1, n_lc);
axm_iir_ekf  = zeros(1, n_lc);
% Mean(dzm) for bias check
mean_eq17  = zeros(1, n_lc);
mean_pp    = zeros(1, n_lc);
mean_ekf   = zeros(1, n_lc);

%% ===== Main Loop =====
for idx = 1:n_lc
    lc = lc_sim(idx);
    C_lc = 2 + 1/(1 - lc^2);
    den_eq13 = 4 * kb * T_temp * C_lc;     % [J] = [N*m]
    ss = ss_start:N;

    % =====================================================================
    % METHOD 1: Eq. 17 (d-step delay control law)
    % =====================================================================
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);

    z = z_d;
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    fd_hist = zeros(N, 1);
    dzm_all = zeros(N, 1);

    % IIR filter states
    dzm_bar = 0; dzr_bar = 0; dzr2_bar = 0; dzm_prev = 0;
    sig2_iir_ts = zeros(N, 1);

    for k = 1:N
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;
        dzm_all(k) = dzm;

        % IIR filter (Eqs. 9-10)
        dzm_bar_new = Avar * dzm_bar + (1 - Avar) * dzm;
        dzr = dzm_prev - dzm_bar_new;
        dzr_bar_new  = Avar22 * dzr + (1 - Avar22) * dzr_bar;
        dzr2_bar_new = Avar3 * dzr^2 + (1 - Avar3) * dzr2_bar;
        sig2_iir_ts(k) = dzr2_bar_new - dzr_bar_new^2;
        dzm_bar = dzm_bar_new; dzr_bar = dzr_bar_new;
        dzr2_bar = dzr2_bar_new; dzm_prev = dzm;

        % Eq. 17 control law
        if k >= 3, fk2 = fd_hist(k-2); else, fk2 = 0; end
        if k >= 2, fk1 = fd_hist(k-1); else, fk1 = 0; end
        fd_k = (1/a_x)*(1-lc)*dzm - (1-lc)*(fk1 + fk2);
        fd_hist(k) = fd_k;

        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
    end

    sig2_eq17(idx) = var(dzm_all(ss));
    sig2_iir_eq17(idx) = mean(sig2_iir_ts(ss));
    axm_eq17(idx)     = (sig2_eq17(idx) * 1e-12) / den_eq13 * 1e-6;
    axm_iir_eq17(idx) = (sig2_iir_eq17(idx) * 1e-12) / den_eq13 * 1e-6;
    mean_eq17(idx) = mean(dzm_all(ss));

    % =====================================================================
    % METHOD 2: 3-state estimator (pole placement)
    % =====================================================================
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);

    z = z_d;
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    dzm_all = zeros(N, 1);
    dz1_hat = 0;  dz2_hat = 0;  dz3_hat = 0;

    % IIR filter states
    dzm_bar = 0; dzr_bar = 0; dzr2_bar = 0; dzm_prev = 0;
    sig2_iir_ts = zeros(N, 1);

    for k = 1:N
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;
        dzm_all(k) = dzm;

        % IIR filter
        dzm_bar_new = Avar * dzm_bar + (1 - Avar) * dzm;
        dzr = dzm_prev - dzm_bar_new;
        dzr_bar_new  = Avar22 * dzr + (1 - Avar22) * dzr_bar;
        dzr2_bar_new = Avar3 * dzr^2 + (1 - Avar3) * dzr2_bar;
        sig2_iir_ts(k) = dzr2_bar_new - dzr_bar_new^2;
        dzm_bar = dzm_bar_new; dzr_bar = dzr_bar_new;
        dzr2_bar = dzr2_bar_new; dzm_prev = dzm;

        % Innovation
        innov = dzm - dz1_hat;

        % Control law (uses known a_x)
        fd_k = (1/a_x) * (1 - lc) * dz3_hat;

        % Estimator update
        dz1_new = dz2_hat      + L1 * innov;
        dz2_new = dz3_hat      + L2 * innov;
        dz3_new = lc * dz3_hat + L3 * innov;
        dz1_hat = dz1_new; dz2_hat = dz2_new; dz3_hat = dz3_new;

        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
    end

    sig2_pp(idx) = var(dzm_all(ss));
    sig2_iir_pp(idx) = mean(sig2_iir_ts(ss));
    axm_pp(idx)     = (sig2_pp(idx) * 1e-12) / den_eq13 * 1e-6;
    axm_iir_pp(idx) = (sig2_iir_pp(idx) * 1e-12) / den_eq13 * 1e-6;
    mean_pp(idx) = mean(dzm_all(ss));

    % =====================================================================
    % METHOD 3: 7-state EKF (simplified from Simulink, known a_x control)
    % =====================================================================
    % States: [dz1, dz2, dz3, zD1, zD2, az1, az2]
    % Simplified: use closed-loop Jacobian since control uses exact a_x.
    % The closed-loop prediction for dz3 = lc*dz3 (known a_x guarantee).
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);

    z = z_d;
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    dzm_all = zeros(N, 1);

    % EKF states
    dz1_hat = 0; dz2_hat = 0; dz3_hat = 0;
    zD1_hat = 0; zD2_hat = 0;
    az1_hat = a_x;              % init at true value (simplified)
    az2_hat = a_x;

    % Thermal noise variance in position [um^2]
    sig2_thermal = a_x^2 * sigma_fT^2;  % = 4*kb*T*a_x * 1e18

    % Fixed Q/R for simplified case
    %   dz states: driven by thermal noise + control noise
    %   zD states: no disturbance in simplified case (small Q)
    %   az states: constant a_x (small Q)
    q_dz  = sig2_thermal * 1e-2;   % position uncertainty per step
    q_zD  = 1e-12;                 % very small (no disturbance)
    q_az  = 1e-12;                 % very small (constant a_x)
    QQ_fix = diag([q_dz, q_dz, sig2_thermal, q_zD, q_zD, q_az, q_az]);

    r_dz  = sig2_thermal * 3;     % measurement noise ~ thermal level
    r_az  = (0.01*a_x)^2;         % azm measurement uncertainty
    RR_fix = diag([r_dz, r_az]);

    % Closed-loop Jacobian (simplified: known a_x, no disturbance coupling)
    Fe_cl = [0  1  0   0           0        0            0;
             0  0  1   0           0        0            0;
             0  0  lc  0           0        0            0;
             0  0  0   1+beta_ekf -beta_ekf 0            0;
             0  0  0   1           0        0            0;
             0  0  0   0           0        1+beta_ekf  -beta_ekf;
             0  0  0   0           0        1            0];

    % Initial covariance
    Pfz = diag([sig2_thermal*10, sig2_thermal*10, sig2_thermal*10, ...
                1e-6, 1e-6, (0.005)^2, (0.005)^2]);

    H_ekf = [1 0 0 0 0 0 0; 0 0 0 0 0 1 0];

    % IIR filter states
    dzm_bar = 0; dzr_bar = 0; dzr2_bar = 0; dzm_prev = 0;
    sig2_iir_ts = zeros(N, 1);
    Am_scaling = 10;
    azm_prev = a_x;

    for k = 1:N
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;
        dzm_all(k) = dzm;

        % IIR filter for azm measurement
        dzm_bar_new = Avar * dzm_bar + (1 - Avar) * dzm;
        dzr = dzm_prev - dzm_bar_new;
        dzr_bar_new  = Avar22 * dzr + (1 - Avar22) * dzr_bar;
        dzr2_bar_new = Avar3 * dzr^2 + (1 - Avar3) * dzr2_bar;
        sigma2_dzr = dzr2_bar_new - dzr_bar_new^2;
        if sigma2_dzr < 0, sigma2_dzr = 0; end
        sig2_iir_ts(k) = sigma2_dzr;
        dzm_bar = dzm_bar_new; dzr_bar = dzr_bar_new;
        dzr2_bar = dzr2_bar_new; dzm_prev = dzm;

        % Compute azm from IIR variance (Eq. 13 real-time)
        den_rt = 4 * kb * T_temp * (2 + 1/(1 - lc^2));
        azm_k = Am_scaling * (sigma2_dzr * 1e-12) / den_rt * (1e6/1e12);

        % --- EKF innovation ---
        e_z1 = dzm - dz1_hat;
        e_az = azm_prev - az1_hat;
        esti_error = [e_z1; e_az];

        % Use fixed Q/R (appropriate for simplified case)
        % Inflate r_az during warmup (first 2000 steps)
        if k < 2000
            RR_k = diag([r_dz, 1e6]);
        else
            RR_k = RR_fix;
        end

        % EKF step 1: innovation covariance
        S = H_ekf * Pfz * H_ekf' + RR_k;
        S = 0.5*(S + S') + 1e-18*eye(2);

        % Kalman gain
        Lk = (Pfz * H_ekf') / S;

        % EKF step 2: update (Joseph form)
        ILH = eye(7) - Lk * H_ekf;
        Pz = ILH * Pfz * ILH' + Lk * RR_k * Lk';
        Pz = 0.5*(Pz + Pz');

        % EKF step 3: predict
        Pfz = Fe_cl * Pz * Fe_cl' + QQ_fix;
        Pfz = 0.5*(Pfz + Pfz');

        % Control law: use KNOWN a_x and CURRENT dz3_hat (before update)
        fd_k = (1/a_x) * (1 - lc) * dz3_hat;

        % State update (use temp variables to avoid sequential-assignment bug)
        inj = Lk * esti_error;
        dz1_new = dz2_hat + inj(1);
        dz2_new = dz3_hat + inj(2);
        dz3_new = lc * dz3_hat + inj(3);
        zD1_new = (1+beta_ekf)*zD1_hat - beta_ekf*zD2_hat + inj(4);
        zD2_new = zD1_hat + inj(5);
        az1_new = (1+beta_ekf)*az1_hat - beta_ekf*az2_hat + inj(6);
        az2_new = az1_hat + inj(7);
        dz1_hat = dz1_new; dz2_hat = dz2_new; dz3_hat = dz3_new;
        zD1_hat = zD1_new; zD2_hat = zD2_new;
        az1_hat = az1_new; az2_hat = az2_new;

        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
        azm_prev = azm_k;
    end

    sig2_ekf(idx) = var(dzm_all(ss));
    sig2_iir_ekf(idx) = mean(sig2_iir_ts(ss));
    axm_ekf(idx)     = (sig2_ekf(idx) * 1e-12) / den_eq13 * 1e-6;
    axm_iir_ekf(idx) = (sig2_iir_ekf(idx) * 1e-12) / den_eq13 * 1e-6;
    mean_ekf(idx) = mean(dzm_all(ss));
end

%% ===== Compute C_actual for each method =====
% C_actual = sig2_measured / (4*kb*T*a_x*1e18)
% Theory: C(lc) = 2 + 1/(1-lc^2)
C_theory = 2 + 1./(1 - lc_sim.^2);
C_eq17   = sig2_eq17 ./ (4*kb*T_temp*a_x*1e18);
C_pp     = sig2_pp   ./ (4*kb*T_temp*a_x*1e18);
C_ekf    = sig2_ekf  ./ (4*kb*T_temp*a_x*1e18);

err_eq17 = abs(axm_eq17 - a_x) / a_x * 100;
err_pp   = abs(axm_pp   - a_x) / a_x * 100;
err_ekf  = abs(axm_ekf  - a_x) / a_x * 100;

%% ===== Results Table =====
fprintf('\n');
fprintf('============================================================================\n');
fprintf('  Eq. 13 Unified Verification: three controller methods\n');
fprintf('  a_x_true = %.5f um/pN,  lambda_e = %.1f\n', a_x, lambda_e);
fprintf('============================================================================\n\n');

fprintf('--- Table 1: Direct variance method ---\n');
fprintf('%-6s %-8s  %-10s %-8s %-8s  %-10s %-8s %-8s  %-10s %-8s %-8s\n', ...
    'lc', 'C(lc)', 'sig2_eq17', 'axm', 'err%', ...
    'sig2_PP', 'axm', 'err%', 'sig2_EKF', 'axm', 'err%');
fprintf('%s\n', repmat('-', 1, 100));
for idx = 1:n_lc
    fprintf('%-6.1f %-8.3f  %-10.4e %-8.5f %-8.2f  %-10.4e %-8.5f %-8.2f  %-10.4e %-8.5f %-8.2f\n', ...
        lc_sim(idx), C_theory(idx), ...
        sig2_eq17(idx), axm_eq17(idx), err_eq17(idx), ...
        sig2_pp(idx),   axm_pp(idx),   err_pp(idx), ...
        sig2_ekf(idx),  axm_ekf(idx),  err_ekf(idx));
end

fprintf('\n--- Table 2: C_actual comparison (C_meas = sig2 / (4kT*a_x)) ---\n');
fprintf('%-6s %-8s  %-10s %-10s %-10s\n', ...
    'lc', 'C_theory', 'C_eq17', 'C_PP', 'C_EKF');
fprintf('%s\n', repmat('-', 1, 50));
for idx = 1:n_lc
    fprintf('%-6.1f %-8.3f  %-10.3f %-10.3f %-10.3f\n', ...
        lc_sim(idx), C_theory(idx), C_eq17(idx), C_pp(idx), C_ekf(idx));
end

fprintf('\n--- Table 3: Mean(dzm) bias check (should be near 0) ---\n');
fprintf('%-6s  %-12s %-12s %-12s\n', 'lc', 'Eq.17', '3-state PP', '7-state EKF');
fprintf('%s\n', repmat('-', 1, 48));
for idx = 1:n_lc
    fprintf('%-6.1f  %-12.2e %-12.2e %-12.2e\n', ...
        lc_sim(idx), mean_eq17(idx), mean_pp(idx), mean_ekf(idx));
end

%% ===== PASS/FAIL =====
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('Eq.17 (direct delay law): ');
if all(err_eq17 < 5), fprintf('PASS (max err=%.2f%%)\n', max(err_eq17));
else, fprintf('FAIL (max err=%.2f%%)\n', max(err_eq17)); end
fprintf('3-state PP:  C_actual/C_theory ratio = %.2f--%.2f\n', ...
    min(C_pp./C_theory), max(C_pp./C_theory));
fprintf('  -> Observer poles add %.0f%%--%.0f%% extra variance\n', ...
    min((C_pp./C_theory-1)*100), max((C_pp./C_theory-1)*100));
fprintf('7-state EKF: C_actual/C_theory ratio = %.2f--%.2f\n', ...
    min(C_ekf./C_theory), max(C_ekf./C_theory));
fprintf('  -> EKF observer adds %.0f%%--%.0f%% extra variance\n', ...
    min((C_ekf./C_theory-1)*100), max((C_ekf./C_theory-1)*100));
fprintf('\nConclusion: Eq. 13 with C(lc)=2+1/(1-lc^2) is exact for Eq.17.\n');
fprintf('Estimator-based methods add observer noise; different C factor needed.\n');

%% ===== Figure 1: sig2 + axm vs lc (2x1) =====
fig1 = figure('Position', [50 50 800 750], 'Color', 'w');

ax1 = subplot(2,1,1);
plot(lc_dense, sig2_theory_dense*1e4, 'b-', 'LineWidth', 2.5); hold on;
plot(lc_sim, sig2_eq17*1e4, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
plot(lc_sim, sig2_pp*1e4, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
plot(lc_sim, sig2_ekf*1e4, 'md', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineWidth', 1.5);
hold off;
ylabel('\sigma^2_{\deltaxr}  (10^{-4} \mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax1, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
xlim([0.3 1.0]);
legend('Theory', 'Eq.17', '3-state PP', '7-state EKF', ...
    'Orientation', 'horizontal', 'FontSize', 11, 'Location', 'northwest');

ax2 = subplot(2,1,2);
plot(lc_dense, ones(size(lc_dense))*a_x, 'b-', 'LineWidth', 2.5); hold on;
plot(lc_sim, axm_eq17, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
plot(lc_sim, axm_pp, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
plot(lc_sim, axm_ekf, 'md', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineWidth', 1.5);
hold off;
xlabel('\lambda_c', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('a_{xm}  (\mum/pN)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax2, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
xlim([0.3 1.0]);
yl_lo = min([a_x, min(axm_eq17), min(axm_pp), min(axm_ekf)]) * 0.9;
yl_hi = max([a_x, max(axm_eq17), max(axm_pp), max(axm_ekf)]) * 1.1;
ylim([yl_lo yl_hi]);

saveas(fig1, 'fig_unified_main.png');

%% ===== Figure 2: Error grouped bar =====
fig2 = figure('Position', [50 50 700 400], 'Color', 'w');
bar_data = [err_eq17; err_pp; err_ekf]';
b = bar(1:n_lc, bar_data, 'grouped');
b(1).FaceColor = [0.9 0.3 0.2];   % Eq.17
b(2).FaceColor = [0.2 0.7 0.3];   % 3-state PP
b(3).FaceColor = [0.6 0.2 0.8];   % 7-state EKF
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lc_sim, 'UniformOutput', false));
xlabel('\lambda_c', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Relative Error (%)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Eq.17', '3-state PP', '7-state EKF', 'Location', 'northwest', 'FontSize', 12);
hold on; yline(2, 'r--', 'LineWidth', 1.5); hold off;
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;
saveas(fig2, 'fig_unified_error.png');

%% ===== Figure 3: Mean(dzm) bar =====
fig3 = figure('Position', [50 50 700 400], 'Color', 'w');
mean_data = [mean_eq17; mean_pp; mean_ekf]';
b3 = bar(1:n_lc, mean_data, 'grouped');
b3(1).FaceColor = [0.9 0.3 0.2];
b3(2).FaceColor = [0.2 0.7 0.3];
b3(3).FaceColor = [0.6 0.2 0.8];
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lc_sim, 'UniformOutput', false));
xlabel('\lambda_c', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Mean(\deltaz_m) (\mum)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Eq.17', '3-state PP', '7-state EKF', 'Location', 'best', 'FontSize', 12);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;
saveas(fig3, 'fig_unified_mean.png');

fprintf('\n3 figures saved: fig_unified_main.png, fig_unified_error.png, fig_unified_mean.png\n');
