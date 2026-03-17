% verify_eq13.m — Standalone verification of Eq. 12 and Eq. 13
%
% Purpose:
%   Under simplified conditions (known gamma_N, sigma^2_nx = 0, probe
%   stationary at 25 um), verify that:
%   1. sigma^2_dxr varies with lambda_c as predicted by Eq. 12
%   2. Eq. 13 correctly compensates this variation, yielding a_xm = a_x
%      = Ts/gamma_N for all lambda_c values.
%
% Simplified equations (sigma^2_nx = 0):
%   Eq. 12: sigma^2_dxr = (2 + 1/(1-lambda_c^2)) * 4*kb*T*a_x
%   Eq. 13: a_xm = sigma^2_dxr / (4*kb*T*(2 + 1/(1-lambda_c^2)))
%
% Controller: Predictor-based controller that compensates for 2-step
%   measurement delay. Without prediction, the naive controller
%   fd = (1/a_x)*(1-lc)*dzm creates 3rd-order closed-loop dynamics
%   with poles that deviate from lambda_c. The predictor restores the
%   intended AR(1) error dynamics: e[k+1] = lambda_c*e[k] + noise.
%
% Two methods for sigma^2_dxr:
%   A) Direct: Var(dzm) over steady-state window (gold standard)
%   B) IIR filter: Eqs. 9-10 stochastic component variance estimate
%      (practical method used in full system; less accurate in stationary
%       case because the filter absorbs stochastic energy into the
%       "deterministic" estimate when no actual deterministic trend exists)

clear; clc; close all;

%% ===== Physical Parameters =====
Ts      = 1/1600;               % sampling time [s]
kb      = 1.3806503e-23;        % Boltzmann constant [J/K]
T_temp  = 310.15;               % temperature [K] (37 C)
R       = 2.25e-6;              % probe radius [m]
eta     = 0.001;                % dynamic viscosity [Pa*s]
gammaN  = 0.0425;               % nominal Stokes drag [pN*s/um]

% Derived quantities
a_x_true = Ts / gammaN;         % theoretical motion gain [um/pN]

% Thermal noise force (use gammaN * 1e-6 for consistency)
sigma2_fT_SI = 4 * kb * T_temp * (gammaN * 1e-6) / Ts;  % [N^2]
sigma_fT_pN  = sqrt(sigma2_fT_SI * 1e24);                % [pN]

% IIR filter coefficients (from Stateflow chart, for Method B)
Avar   = 0.45;    % deterministic component filter
Avar22 = 0.05;    % stochastic component mean filter
Avar3  = 0.05;    % stochastic component mean-square filter

%% ===== Simulation Parameters =====
lambda_c_list = [0.3, 0.5, 0.7, 0.9];
n_lc = length(lambda_c_list);

N_steps      = 80000;                    % 50 seconds at 1600 Hz
steady_start = 40000;                    % use last 25 s for steady-state
z_d          = 25;                       % target position [um]

%% ===== Storage =====
sigma2_direct  = zeros(1, n_lc);   % Method A: Var(dzm)
sigma2_iir     = zeros(1, n_lc);   % Method B: IIR filter estimate
sigma2_theory  = zeros(1, n_lc);   % Eq. 12 theory
axm_direct     = zeros(1, n_lc);   % a_xm from Method A
axm_iir        = zeros(1, n_lc);   % a_xm from Method B
axm_ts_direct  = cell(1, n_lc);    % time series (running variance)
axm_ts_iir     = cell(1, n_lc);    % time series (IIR)

%% ===== Main Loop over lambda_c =====
for idx = 1:n_lc
    lc = lambda_c_list(idx);
    C_lc = 2 + 1 / (1 - lc^2);
    den_eq13 = 4 * kb * T_temp * C_lc;    % [J] = [N*m]

    % Eq. 12 theoretical sigma^2_dxr [um^2]
    % Unit chain: 4*kb*T [J=N*m] * a_x [um/pN = 1e6 m / 1e12 N = 1e-6 m/N]
    %   = 4*kb*T*a_x*1e-6 [m^2] = 4*kb*T*a_x*1e-6*1e12 [um^2]
    %   = 4*kb*T*a_x*1e6... NO. Let me be precise:
    % a_x [um/pN] * (1e-6 m/um) / (1e-12 N/pN) = a_x * 1e6 [m/N]
    % 4*kb*T [N*m] * a_x*1e6 [m/N] = 4*kb*T*a_x*1e6 [m^2]
    % [m^2] * (1e6 um/m)^2 = * 1e12 => 4*kb*T*a_x*1e18 [um^2]
    sigma2_theory(idx) = C_lc * 4 * kb * T_temp * a_x_true * 1e18;

    % Generate thermal noise
    rng(42 + idx);
    fT = sigma_fT_pN * randn(N_steps, 1);  % [pN]

    % State initialization
    z      = z_d;
    z_hist = zeros(N_steps, 1);  z_hist(1) = z_d;
    fd_hist = zeros(N_steps, 1);
    dzm_all = zeros(N_steps, 1);

    % IIR filter states
    dzm_bar_prev  = 0;
    dzr_bar_prev  = 0;
    dzr2_bar_prev = 0;
    dzm_prev      = 0;
    sig2_iir_ts   = zeros(N_steps, 1);

    % Running variance accumulators (for Method A time series)
    dzm_sum  = 0;
    dzm2_sum = 0;
    run_var_ts = zeros(N_steps, 1);

    for k = 1:N_steps
        % --- Measurement with two-step delay ---
        if k >= 3
            z_delayed = z_hist(k-2);
        elseif k == 2
            z_delayed = z_hist(1);
        else
            z_delayed = z_d;
        end
        dzm = z_d - z_delayed;
        dzm_all(k) = dzm;

        % --- Method B: IIR Filters (Eqs. 9-10) ---
        dzm_bar = Avar * dzm_bar_prev + (1 - Avar) * dzm;
        dzr = dzm_prev - dzm_bar;
        dzr_bar  = Avar22 * dzr + (1 - Avar22) * dzr_bar_prev;
        dzr2_bar = Avar3 * dzr^2 + (1 - Avar3) * dzr2_bar_prev;
        sigma2_dzr = dzr2_bar - dzr_bar^2;
        sig2_iir_ts(k) = sigma2_dzr;
        dzm_bar_prev  = dzm_bar;
        dzr_bar_prev  = dzr_bar;
        dzr2_bar_prev = dzr2_bar;
        dzm_prev      = dzm;

        % --- Method A: Running variance of dzm ---
        if k > 100  % wait for transient to pass
            n_so_far = k - 100;
            dzm_sum  = dzm_sum + dzm;
            dzm2_sum = dzm2_sum + dzm^2;
            run_var_ts(k) = dzm2_sum/n_so_far - (dzm_sum/n_so_far)^2;
        end

        % --- Predictor-based controller (compensates 2-step delay) ---
        % Predict current state from delayed measurement + known controls
        %   z_hat[k] = z[k-2] + a_x*(fd[k-2] + fd[k-1])
        %   e_hat[k] = z_d - z_hat[k]
        if k >= 3, fd_km2 = fd_hist(k-2); else, fd_km2 = 0; end
        if k >= 2, fd_km1 = fd_hist(k-1); else, fd_km1 = 0; end
        e_hat = dzm - a_x_true * (fd_km2 + fd_km1);
        fd_k  = (1 / a_x_true) * (1 - lc) * e_hat;
        fd_hist(k) = fd_k;

        % --- Plant dynamics ---
        z_new = z + a_x_true * (fd_k + fT(k));
        if k < N_steps
            z_hist(k+1) = z_new;
        end
        z = z_new;
    end

    % --- Steady-state analysis ---
    ss = steady_start:N_steps;

    % Method A: direct variance
    sigma2_direct(idx) = var(dzm_all(ss));
    axm_direct(idx)    = (sigma2_direct(idx) * 1e-12) / den_eq13 * 1e-6;

    % Method B: IIR filter mean
    sigma2_iir(idx) = mean(sig2_iir_ts(ss));
    axm_iir(idx)    = (sigma2_iir(idx) * 1e-12) / den_eq13 * 1e-6;

    % Time series for Eq. 13 (convert running variance to a_xm)
    axm_ts_direct{idx} = (run_var_ts * 1e-12) / den_eq13 * 1e-6;
    axm_ts_iir{idx}    = (sig2_iir_ts * 1e-12) / den_eq13 * 1e-6;

    fprintf('--- lambda_c = %.1f  C(lc) = %.3f ---\n', lc, C_lc);
    fprintf('  sigma2 theory=%.4e  direct=%.4e  IIR=%.4e  [um^2]\n', ...
        sigma2_theory(idx), sigma2_direct(idx), sigma2_iir(idx));
    fprintf('  a_xm   theory=%.5f  direct=%.5f (%.1f%%)  IIR=%.5f (%.1f%%)\n', ...
        a_x_true, axm_direct(idx), ...
        abs(axm_direct(idx)-a_x_true)/a_x_true*100, ...
        axm_iir(idx), ...
        abs(axm_iir(idx)-a_x_true)/a_x_true*100);
end

%% ===== Figure 1: sigma^2_dxr vs lambda_c (Verify Eq. 12) =====
t_vec = (0:N_steps-1)' * Ts;

figure('Name', 'Eq.12 Verification: sigma^2_dxr vs lambda_c');
bar_data = [sigma2_theory; sigma2_direct]';
b = bar(1:n_lc, bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.9 0.3 0.2];
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lambda_c_list, 'UniformOutput', false));
xlabel('\lambda_c');
ylabel('\sigma^2_{dxr} (\mum^2)');
title('Eq. 12 Verification: \sigma^2_{dxr} vs \lambda_c');
legend('Theory (Eq. 12)', 'Measured Var(dzm)', 'Location', 'northwest');
grid on;

%% ===== Figure 2: a_xm time series (Verify Eq. 13 convergence) =====
figure('Name', 'Eq.13 Verification: a_xm convergence');
colors = {'b', 'r', [0 0.6 0], [0.8 0.4 0]};
hold on;
for idx = 1:n_lc
    plot(t_vec, axm_ts_direct{idx}, 'Color', colors{idx}, ...
        'DisplayName', sprintf('\\lambda_c=%.1f', lambda_c_list(idx)));
end
yline(a_x_true, 'k--', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('a_x true = %.5f', a_x_true));
hold off;
xlabel('Time (s)');
ylabel('a_{xm} (\mum/pN)');
title('Eq. 13 Verification: a_{xm} convergence (direct variance method)');
legend('Location', 'best');
grid on;
ylim([0, a_x_true * 3]);

%% ===== Figure 3: Steady-state a_xm bar chart + error =====
figure('Name', 'Eq.13 Summary: steady-state a_xm');

subplot(2,1,1);
bar_axm = [repmat(a_x_true, 1, n_lc); axm_direct]';
b2 = bar(1:n_lc, bar_axm);
b2(1).FaceColor = [0.2 0.4 0.8];
b2(2).FaceColor = [0.9 0.3 0.2];
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lambda_c_list, 'UniformOutput', false));
xlabel('\lambda_c');
ylabel('a_{xm} (\mum/pN)');
title('Steady-state a_{xm} vs theory (direct variance)');
legend('a_x true', 'a_{xm} measured', 'Location', 'best');
grid on;

subplot(2,1,2);
rel_err = abs(axm_direct - a_x_true) / a_x_true * 100;
bar(1:n_lc, rel_err, 'FaceColor', [0.9 0.6 0.1]);
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lambda_c_list, 'UniformOutput', false));
xlabel('\lambda_c');
ylabel('Relative Error (%)');
title('Relative error of a_{xm}');
yline(5, 'r--', 'LineWidth', 1.5);
grid on;

%% ===== Summary Table =====
fprintf('\n========== VERIFICATION SUMMARY ==========\n');
fprintf('Theoretical a_x = Ts/gammaN = %.5f um/pN\n\n', a_x_true);
fprintf('--- Method A: Direct Variance (Var of tracking error) ---\n');
fprintf('%-10s %-8s %-14s %-14s %-14s %-10s\n', ...
    'lambda_c', 'C(lc)', 'sig2_theory', 'sig2_meas', 'axm_meas', 'Error(%)');
fprintf('%s\n', repmat('-', 1, 70));
for idx = 1:n_lc
    lc = lambda_c_list(idx);
    C_lc = 2 + 1 / (1 - lc^2);
    fprintf('%-10.1f %-8.3f %-14.4e %-14.4e %-14.5f %-10.2f\n', ...
        lc, C_lc, sigma2_theory(idx), sigma2_direct(idx), ...
        axm_direct(idx), rel_err(idx));
end

fprintf('\n--- Method B: IIR Filter (Eqs. 9-10 stochastic residual) ---\n');
rel_err_iir = abs(axm_iir - a_x_true) / a_x_true * 100;
fprintf('%-10s %-14s %-14s %-10s\n', 'lambda_c', 'sig2_iir', 'axm_iir', 'Error(%)');
fprintf('%s\n', repmat('-', 1, 48));
for idx = 1:n_lc
    fprintf('%-10.1f %-14.4e %-14.5f %-10.2f\n', ...
        lambda_c_list(idx), sigma2_iir(idx), axm_iir(idx), rel_err_iir(idx));
end

fprintf('\nNote: IIR filter underestimates sigma^2_dxr in the stationary case\n');
fprintf('because it absorbs stochastic energy into the deterministic estimate.\n');
fprintf('In the full system with a moving probe, the filter correctly separates\n');
fprintf('deterministic and stochastic components.\n');

fprintf('%s\n', repmat('=', 1, 70));
if all(rel_err < 5)
    fprintf('RESULT: PASS - All direct-method relative errors < 5%%\n');
else
    fprintf('RESULT: PARTIAL - Some errors >= 5%%\n');
    fprintf('  Failed: '); fprintf('%.1f ', lambda_c_list(rel_err >= 5)); fprintf('\n');
end
