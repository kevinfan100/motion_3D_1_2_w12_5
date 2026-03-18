% verify_eq13_estimator.m — Eq. 13 verification: two delay-compensation methods
%
% Method 1: 3-state estimator (dz1->dz2->dz3 + pole placement)
% Method 2: d-step delay control law (Eq. 17, expanded form)
%
% Simplified conditions: known gamma_N, no disturbance, stationary probe.
% Goal: verify both methods recover a_x via Eq. 13, and compare accuracy.

clear; clc; close all;

%% ===== Physical Parameters =====
Ts      = 1/1600;            % sampling time [s]
kb      = 1.3806503e-23;     % Boltzmann constant [J/K]
T_temp  = 310.15;            % temperature [K]
R_probe = 2.25e-6;           % probe radius [m]
eta     = 0.001;             % viscosity [Pa*s]
gammaN  = 0.0425;            % friction coefficient [pN*s/um]

a_x      = Ts / gammaN;                                    % mobility [um/pN]
sigma_fT  = sqrt(4 * kb * T_temp * (gammaN * 1e-6) / Ts * 1e24); % thermal force std [pN]

%% ===== Estimator: Pole Placement =====
% Pole placement on TRUE error dynamics F = A + G - LC, where
% A+G has (3,3) = lc + (1-lc) = 1 due to controller-observer coupling.
% Gains are independent of lc.
lambda_e = 0.3;
L1 = 1 - 3*lambda_e;
L2 = 1 - 3*lambda_e + 3*lambda_e^2;
L3 = (1 - lambda_e)^3;
fprintf('Estimator pole placement: lambda_e = %.1f\n\n', lambda_e);

%% ===== Simulation Parameters =====
lc_sim   = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9];  % lambda_c > lambda_e = 0.3
n_lc     = length(lc_sim);
N        = 80000;
ss_start = 40000;
z_d      = 25;               % desired position [um]

% Dense theory curve
lc_dense = 0.1:0.005:0.95;
C_dense  = 2 + 1./(1 - lc_dense.^2);
sig2_theory_dense = C_dense * 4 * kb * T_temp * a_x * 1e18;  % [um^2]

%% ===== Storage =====
sig2_est   = zeros(1, n_lc);   % estimator method
sig2_eq17  = zeros(1, n_lc);   % d-step delay control law (Eq. 17)
axm_est    = zeros(1, n_lc);   % recovered a_x from estimator [um/pN]
axm_eq17   = zeros(1, n_lc);   % recovered a_x from Eq. 17 [um/pN]

%% ===== Main Loop =====
for idx = 1:n_lc
    lc = lc_sim(idx);
    C_lc = 2 + 1/(1 - lc^2);
    den_eq13 = 4 * kb * T_temp * C_lc;

    % --- Method 1: 3-state estimator ---
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);          % thermal force [pN]

    z = z_d;                               % position [um]
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    dzm_all_est = zeros(N, 1);

    % Estimator states: dz = z_d - z (displacement error [um])
    dz1_hat = 0;  dz2_hat = 0;  dz3_hat = 0;

    for k = 1:N
        % Delayed measurement (d = 2 steps)
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;                % [um]
        dzm_all_est(k) = dzm;

        % Innovation
        innov = dzm - dz1_hat;

        % Control law: fd = (1/a_x)*[(z_d[k]-z_d[k-1]) + (1-lc)*dz3_hat]
        % Feedforward term is 0 for stationary z_d.
        fd_k = (1/a_x) * (1 - lc) * dz3_hat;   % [pN]

        % Estimator update
        dz1_new = dz2_hat      + L1 * innov;
        dz2_new = dz3_hat      + L2 * innov;
        dz3_new = lc * dz3_hat + L3 * innov;
        dz1_hat = dz1_new;
        dz2_hat = dz2_new;
        dz3_hat = dz3_new;

        % Plant dynamics: dz = a_x * (fd + fT) [um]
        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
    end

    ss = ss_start:N;
    sig2_est(idx) = var(dzm_all_est(ss));                     % [um^2]
    axm_est(idx) = (sig2_est(idx) * 1e-12) / den_eq13 * 1e-6; % [um/pN]

    % --- Method 2: d-step delay control law, Eq. 17 (same noise) ---
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);          % thermal force [pN]

    z = z_d;
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    fd_hist = zeros(N, 1);
    dzm_all_eq17 = zeros(N, 1);

    for k = 1:N
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;                % [um]
        dzm_all_eq17(k) = dzm;

        % Eq. 17: d-step delay control law (d=2, stationary)
        %   fd[k] = (1/a_x)*(1-lc)*dzm[k] - (1-lc)*sum(fd[k-i], i=1..d)
        if k >= 3, fk2 = fd_hist(k-2); else, fk2 = 0; end
        if k >= 2, fk1 = fd_hist(k-1); else, fk1 = 0; end
        fd_k = (1/a_x)*(1-lc)*dzm - (1-lc)*(fk1 + fk2);  % [pN]
        fd_hist(k) = fd_k;

        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
    end

    sig2_eq17(idx) = var(dzm_all_eq17(ss));                     % [um^2]
    axm_eq17(idx) = (sig2_eq17(idx) * 1e-12) / den_eq13 * 1e-6; % [um/pN]
end

%% ===== Results Table =====
% axm_est / axm_eq17: a_x recovered via Eq. 13 [um/pN]
% err_est / err_eq17: relative error vs true a_x [%]
rel_err_est  = abs(axm_est  - a_x) / a_x * 100;
rel_err_eq17 = abs(axm_eq17 - a_x) / a_x * 100;

fprintf('==========================================================\n');
fprintf('  Eq. 13 verification: 3-state estimator vs Eq. 17\n');
fprintf('==========================================================\n\n');
fprintf('%-6s  %-12s %-12s %-10s %-10s\n', ...
    'lc', 'axm_est', 'axm_eq17', 'err_est', 'err_eq17');
fprintf('%s\n', repmat('-', 1, 54));
for idx = 1:n_lc
    fprintf('%-6.1f  %-12.5f %-12.5f %-10.2f %-10.2f\n', ...
        lc_sim(idx), axm_est(idx), axm_eq17(idx), ...
        rel_err_est(idx), rel_err_eq17(idx));
end
fprintf('\na_x_true = %.5f um/pN\n', a_x);

%% ===== Lambda_e sensitivity test (lc = 0.7) =====
% Test how estimator pole lambda_e affects a_x estimation accuracy
lc_test = 0.7;
le_list = [0.1, 0.2, 0.3, 0.5, 0.7];
axm_le = zeros(size(le_list));
C_test = 2 + 1/(1 - lc_test^2);
den_test = 4 * kb * T_temp * C_test;

for j = 1:length(le_list)
    le = le_list(j);
    L1t = 1 - 3*le;
    L2t = 1 - 3*le + 3*le^2;
    L3t = (1 - le)^3;

    rng(99);
    fT = sigma_fT * randn(N, 1);
    z = z_d; z_hist = zeros(N,1); z_hist(1) = z_d;
    dzm_all = zeros(N,1);
    dz1h = 0; dz2h = 0; dz3h = 0;

    for k = 1:N
        if k>=3, z_del=z_hist(k-2); else, z_del=z_d; end
        dzm = z_d - z_del; dzm_all(k) = dzm;
        innov = dzm - dz1h;
        fd_k = (1/a_x)*(1-lc_test)*dz3h;
        dz1h_n = dz2h + L1t*innov;
        dz2h_n = dz3h + L2t*innov;
        dz3h_n = lc_test*dz3h + L3t*innov;
        dz1h = dz1h_n; dz2h = dz2h_n; dz3h = dz3h_n;
        z_new = z + a_x*(fd_k + fT(k));
        if k<N, z_hist(k+1)=z_new; end; z=z_new;
    end
    axm_le(j) = (var(dzm_all(ss))*1e-12)/den_test*1e-6;
end

fprintf('\n--- Lambda_e sensitivity (lc = 0.7) ---\n');
fprintf('%-10s %-12s %-10s\n', 'lambda_e', 'a_xm', 'error');
for j = 1:length(le_list)
    fprintf('%-10.1f %-12.5f %-10.2f%%\n', le_list(j), axm_le(j), ...
        abs(axm_le(j)-a_x)/a_x*100);
end

%% ===== Figure 1: Main result =====
fig1 = figure('Position', [50 50 700 750], 'Color', 'w');

ax1 = subplot(2,1,1);
h1a = plot(lc_dense, sig2_theory_dense*1e4, 'b-', 'LineWidth', 2.5);
hold on;
h1b = plot(lc_sim, sig2_est*1e4, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
hold off;
xlabel('\lambda_c', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('\sigma^2_{\deltaxr}  (10^{-4} \mum^2)', 'FontSize', 16, 'FontWeight', 'bold');
set(ax1, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2, ...
    'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
xlim([0.1 1.0]);
set(ax1, 'Position', [0.13 0.53 0.82 0.38]);

ax2 = subplot(2,1,2);
plot(lc_dense, ones(size(lc_dense))*a_x, 'b-', 'LineWidth', 2.5);
hold on;
plot(lc_sim, axm_est, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
hold off;
xlabel('\lambda_c', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('a_{xm}  (\mum/pN)', 'FontSize', 16, 'FontWeight', 'bold');
set(ax2, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2, ...
    'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
xlim([0.1 1.0]);
% Auto ylim to include both theory line and simulation data
yl_lo = min([a_x, min(axm_est)]) * 0.9;
yl_hi = max([a_x, max(axm_est)]) * 1.1;
ylim([yl_lo yl_hi]);
set(ax2, 'Position', [0.13 0.08 0.82 0.38]);

lg = legend(ax1, [h1a, h1b], {'Theory', 'Simulation'}, ...
    'Orientation', 'horizontal', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Box', 'on', 'LineWidth', 1.5);
set(lg, 'Position', [0.32 0.93 0.36 0.04]);

saveas(fig1, 'fig_estimator_main.png');

%% ===== Figure 2: Error bar chart =====
fig2 = figure('Position', [50 50 600 350], 'Color', 'w');
bar(lc_sim, rel_err_est, 0.6, 'FaceColor', [0.2 0.5 0.8]);
hold on; yline(2, 'r--', 'LineWidth', 1.5);
text(0.85, 2.3, '2%', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
hold off;
xlabel('\lambda_c', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Relative Error (%)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, ...
    'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
ylim([0 max(rel_err_est)*1.2]);
saveas(fig2, 'fig_estimator_error.png');

%% ===== Figure 3: Lambda_e sensitivity =====
fig3 = figure('Position', [50 50 600 350], 'Color', 'w');
plot(le_list, axm_le, 'ro-', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2);
hold on;
yline(a_x, 'b-', 'LineWidth', 2.5);
hold off;
xlabel('\lambda_e', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('a_{xm}  (\mum/pN)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, ...
    'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
ylim([0.012 max(axm_le)*1.1]);
lg3 = legend('Simulation', 'a_x true', 'Location', 'north', ...
    'Orientation', 'horizontal', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Box', 'on', 'LineWidth', 1.5);
saveas(fig3, 'fig_estimator_lambda_e.png');

fprintf('\n3 figures saved.\n');
