% verify_eq13_estimator.m — Eq. 13 verification using 3-state estimator
%
% Uses the paper's estimator structure (dx1→dx2→dx3) with pole placement
% to compensate the 2-step measurement delay.
% Simplified conditions: known gamma_N, no disturbance, stationary probe.
%
% Comparison with Smith Predictor (verify_eq13.m) is included.

clear; clc; close all;

%% ===== Physical Parameters =====
Ts      = 1/1600;
kb      = 1.3806503e-23;
T_temp  = 310.15;
R_probe = 2.25e-6;
eta     = 0.001;
gammaN  = 0.0425;

a_x      = Ts / gammaN;
gammaN_SI = 6 * pi * eta * R_probe;
sigma_fT  = sqrt(4 * kb * T_temp * gammaN_SI / Ts * 1e24);

%% ===== Estimator: Pole Placement =====
lambda_e = 0.3;   % estimator pole

L1 = 1 - 3*lambda_e;
L2 = 1 - 3*lambda_e + 3*lambda_e^2;
L3 = (1 - lambda_e)^3;

fprintf('Estimator pole placement: lambda_e = %.1f\n', lambda_e);
fprintf('L = [%.3f, %.3f, %.3f]\n\n', L1, L2, L3);

%% ===== Simulation Parameters =====
lc_sim   = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
n_lc     = length(lc_sim);
N        = 80000;
ss_start = 40000;
z_d      = 25;

% Dense theory curve
lc_dense = 0.1:0.005:0.95;
C_dense  = 2 + 1./(1 - lc_dense.^2);
sig2_theory_dense = C_dense * 4 * kb * T_temp * a_x * 1e18;

%% ===== Storage =====
sig2_est  = zeros(1, n_lc);   % estimator method
sig2_sp   = zeros(1, n_lc);   % Smith Predictor (comparison)
axm_est   = zeros(1, n_lc);
axm_sp    = zeros(1, n_lc);

%% ===== Main Loop =====
for idx = 1:n_lc
    lc = lc_sim(idx);
    C_lc = 2 + 1/(1 - lc^2);
    den_eq13 = 4 * kb * T_temp * C_lc;

    % --- Method 1: 3-state estimator ---
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);

    z = z_d;
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    dzm_all_est = zeros(N, 1);

    % Estimator states
    dx1_hat = 0;  dx2_hat = 0;  dx3_hat = 0;

    for k = 1:N
        % Delayed measurement
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;
        dzm_all_est(k) = dzm;

        % Innovation
        innov = dzm - dx1_hat;

        % Control law (uses dx3_hat)
        fd_k = (1/a_x) * (1 - lc) * dx3_hat;

        % Estimator update
        dx1_new = dx2_hat      + L1 * innov;
        dx2_new = dx3_hat      + L2 * innov;
        dx3_new = lc * dx3_hat + L3 * innov;
        dx1_hat = dx1_new;
        dx2_hat = dx2_new;
        dx3_hat = dx3_new;

        % Plant
        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
    end

    ss = ss_start:N;
    sig2_est(idx) = var(dzm_all_est(ss));
    axm_est(idx) = (sig2_est(idx) * 1e-12) / den_eq13 * 1e-6;

    % --- Method 2: Smith Predictor (comparison, same noise) ---
    rng(42 + idx);
    fT = sigma_fT * randn(N, 1);

    z = z_d;
    z_hist = zeros(N, 1);  z_hist(1) = z_d;
    fd_hist = zeros(N, 1);
    dzm_all_sp = zeros(N, 1);

    for k = 1:N
        if k >= 3, z_del = z_hist(k-2); else, z_del = z_d; end
        dzm = z_d - z_del;
        dzm_all_sp(k) = dzm;

        if k >= 3, fk2 = fd_hist(k-2); else, fk2 = 0; end
        if k >= 2, fk1 = fd_hist(k-1); else, fk1 = 0; end
        e_hat = dzm - a_x * (fk2 + fk1);
        fd_k = (1/a_x) * (1 - lc) * e_hat;
        fd_hist(k) = fd_k;

        z_new = z + a_x * (fd_k + fT(k));
        if k < N, z_hist(k+1) = z_new; end
        z = z_new;
    end

    sig2_sp(idx) = var(dzm_all_sp(ss));
    axm_sp(idx) = (sig2_sp(idx) * 1e-12) / den_eq13 * 1e-6;
end

%% ===== Results Table =====
rel_err_est = abs(axm_est - a_x) / a_x * 100;
rel_err_sp  = abs(axm_sp  - a_x) / a_x * 100;

fprintf('=====================================================\n');
fprintf('  Eq. 13 Verification: Estimator vs Smith Predictor\n');
fprintf('=====================================================\n\n');
fprintf('%-6s  %-12s %-12s %-10s %-10s\n', ...
    'lc', 'axm_est', 'axm_SP', 'err_est', 'err_SP');
fprintf('%s\n', repmat('-', 1, 52));
for idx = 1:n_lc
    fprintf('%-6.1f  %-12.5f %-12.5f %-10.2f %-10.2f\n', ...
        lc_sim(idx), axm_est(idx), axm_sp(idx), ...
        rel_err_est(idx), rel_err_sp(idx));
end
fprintf('\na_x_true = %.5f um/pN\n', a_x);

%% ===== Lambda_e sensitivity test (lc = 0.7) =====
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
    dx1h = 0; dx2h = 0; dx3h = 0;

    for k = 1:N
        if k>=3, z_del=z_hist(k-2); else, z_del=z_d; end
        dzm = z_d - z_del; dzm_all(k) = dzm;
        innov = dzm - dx1h;
        fd_k = (1/a_x)*(1-lc_test)*dx3h;
        dx1h_n = dx2h + L1t*innov;
        dx2h_n = dx3h + L2t*innov;
        dx3h_n = lc_test*dx3h + L3t*innov;
        dx1h = dx1h_n; dx2h = dx2h_n; dx3h = dx3h_n;
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

%% ===== Figure 1: Main result (same format as before) =====
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
xlim([0.1 1.0]); ylim([0.012 0.018]);
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
ylim([0 3]);
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
ylim([0.012 0.018]);
lg3 = legend('Simulation', 'a_x true', 'Location', 'north', ...
    'Orientation', 'horizontal', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Box', 'on', 'LineWidth', 1.5);
saveas(fig3, 'fig_estimator_lambda_e.png');

fprintf('\n3 figures saved.\n');
