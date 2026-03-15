% verify_eq13.m — Standalone verification of Eq. 12 and Eq. 13
%
% Purpose:
%   Under simplified conditions (known gamma_N, sigma^2_nx = 0, probe
%   stationary at 25 um), verify that:
%   1. sigma^2_dxr varies with lambda_c as predicted by Eq. 12
%   2. Eq. 13 correctly compensates this variation, yielding a_xm ≈ a_x
%      = Ts/gamma_N for all lambda_c values.
%
% Simplified equations (sigma^2_nx = 0):
%   Eq. 12: sigma^2_dxr = (2 + 1/(1-lambda_c^2)) * 4*kb*T*a_x
%   Eq. 13: a_xm = sigma^2_dxr / (4*kb*T*(2 + 1/(1-lambda_c^2)))
%
% The control loop simulates 1D (z-direction) plant dynamics with:
%   - Two-step measurement delay
%   - IIR filters for deterministic/stochastic component separation
%   - Known gamma_N controller (no estimation needed)

clear; clc; close all;
rng(42);  % fixed seed for reproducibility

%% ===== Physical Parameters =====
Ts      = 1/1600;               % sampling time [s]
kb      = 1.3806503e-23;        % Boltzmann constant [J/K]
T_temp  = 310.15;               % temperature [K] (37 C)
R       = 2.25e-6;              % probe radius [m]
eta     = 0.001;                % dynamic viscosity [Pa*s]
gammaN  = 0.0425;               % nominal Stokes drag [pN*s/um]

% Derived quantities
a_x_true = Ts / gammaN;         % theoretical motion gain [um/pN] ≈ 0.01471

% Thermal noise force (SI units first, then convert)
gammaN_SI = 6 * pi * eta * R;   % [N*s/m] = 4.2412e-8
sigma2_fT_SI = 4 * kb * T_temp * gammaN_SI / Ts;  % [N^2]
sigma2_fT_pN = sigma2_fT_SI * 1e24;               % [pN^2]
sigma_fT_pN  = sqrt(sigma2_fT_pN);                % [pN]

% IIR filter coefficients (from Stateflow chart)
Avar   = 0.45;    % deterministic component filter
Avar22 = 0.05;    % stochastic component mean filter
Avar3  = 0.05;    % stochastic component mean-square filter

%% ===== Simulation Parameters =====
lambda_c_list = [0.3, 0.5, 0.7, 0.9];
n_lc = length(lambda_c_list);

sim_time  = 5;                   % [s]
N_steps   = round(sim_time / Ts); % 8000 steps
steady_start = round(3 / Ts);    % last 2 seconds for steady-state

% Target position [um]
z_d = 25;

%% ===== Storage for results =====
sigma2_dxr_meas   = zeros(1, n_lc);  % measured sigma^2_dxr (steady-state)
axm_meas          = zeros(1, n_lc);  % measured a_xm from Eq. 13 (steady-state)
sigma2_dxr_theory = zeros(1, n_lc);  % theoretical sigma^2_dxr from Eq. 12
axm_time_all      = cell(1, n_lc);   % full time series of a_xm

%% ===== Main Loop over lambda_c =====
for idx = 1:n_lc
    lc = lambda_c_list(idx);
    fprintf('--- lambda_c = %.1f ---\n', lc);

    % Theoretical values
    C_lc = 2 + 1 / (1 - lc^2);  % coefficient in Eq. 12/13
    % Eq. 12 theoretical: sigma^2_dxr in um^2
    %   4*kb*T*a_x has units: [J/K]*[K]*[um/pN] = [J]*[um/pN]
    %   Need to convert to um^2:
    %   [J]*[um/pN] = [N*m]*[um/pN] = [m]*[um] = 1e-6 * [um^2]
    %   Actually let's be careful:
    %   4*kb*T in SI = [J] = [N*m]
    %   a_x in um/pN = a_x * 1e-6 / 1e-12 [m/N] = a_x * 1e6 [m/N]
    %   So 4*kb*T*a_x [N*m]*[m/N] = [m^2]
    %   Convert to um^2: * 1e12
    sigma2_dxr_theory(idx) = C_lc * 4 * kb * T_temp * (a_x_true * 1e6) * 1e12;
    % Simplify: 4*kb*T * a_x_true * 1e6 * 1e12 = 4*kb*T * a_x_true * 1e18
    % But let's keep it clear:
    % a_x_true [um/pN] -> SI: a_x_true * 1e-6 / 1e-12 = a_x_true * 1e6 [m/N]
    % 4*kb*T [J] * a_x_true*1e6 [m/N] = 4*kb*T*a_x_true*1e6 [m^2]
    % to um^2: * 1e12 -> 4*kb*T*a_x_true*1e18 [um^2]
    sigma2_dxr_theory(idx) = C_lc * 4 * kb * T_temp * a_x_true * 1e18;

    % Pre-generate thermal noise for this run
    rng(42 + idx);  % different but reproducible per lambda_c
    fT_all = sigma_fT_pN * randn(N_steps, 1);  % [pN]

    % State variables
    z    = z_d;         % probe position [um], start at target
    z_hist = zeros(N_steps, 1);  % position history
    z_hist(1) = z;

    % Measurement delay buffer (need z[k-2])
    z_buf = [z_d; z_d; z_d];  % circular buffer for delay

    % IIR filter states
    dzm_bar_prev  = 0;   % filtered deterministic component
    dzr_bar_prev  = 0;   % filtered stochastic mean
    dzr2_bar_prev = 0;   % filtered stochastic mean-square
    dzm_prev      = 0;   % previous measurement error

    % Output storage
    sigma2_dzr_ts = zeros(N_steps, 1);  % sigma^2_dxr time series
    axm_ts        = zeros(N_steps, 1);  % a_xm time series

    % Eq. 13 denominator (constant for given lambda_c)
    den_eq13 = 4 * kb * T_temp * C_lc;  % [J] = [N*m]

    for k = 1:N_steps
        % --- Measurement with two-step delay ---
        % At step k, we measure z from step k-2
        if k >= 3
            z_delayed = z_hist(k-2);
        elseif k == 2
            z_delayed = z_hist(1);
        else
            z_delayed = z_d;
        end
        dzm = z_d - z_delayed;  % measurement error [um] (sigma^2_nx = 0)

        % --- IIR Filters (Eqs. 9-10) ---
        % Eq 9: deterministic component (low-pass filter on dzm)
        dzm_bar = Avar * dzm_bar_prev + (1 - Avar) * dzm;

        % Stochastic component: use dzm[k-1]
        dzr = dzm_prev - dzm_bar;

        % Eq 10: variance estimation
        dzr_bar  = Avar22 * dzr + (1 - Avar22) * dzr_bar_prev;   % mean
        dzr2_bar = Avar3 * dzr^2 + (1 - Avar3) * dzr2_bar_prev;  % mean-square
        sigma2_dzr = dzr2_bar - dzr_bar^2;                        % variance [um^2]

        % --- Eq. 13: Measured motion gain ---
        % sigma2_dzr [um^2] -> [m^2]: * 1e-12
        % den_eq13 [N*m]
        % ratio: [m^2]/[N*m] = [m/N]
        % to um/pN: [m/N] * (1e6 um/m) / (1e12 pN/N) = * 1e-6
        if k > 10  % avoid division issues at start
            axm_k = (sigma2_dzr * 1e-12) / den_eq13 * 1e-6;  % [um/pN]
        else
            axm_k = 0;
        end

        % Store
        sigma2_dzr_ts(k) = sigma2_dzr;
        axm_ts(k)        = axm_k;

        % Update IIR filter states
        dzm_bar_prev  = dzm_bar;
        dzr_bar_prev  = dzr_bar;
        dzr2_bar_prev = dzr2_bar;
        dzm_prev      = dzm;

        % --- Control law (Eq. 6, simplified with known a_x) ---
        fd = (1 / a_x_true) * (1 - lc) * dzm;  % [pN]

        % --- Plant dynamics (Eq. 5) ---
        % z[k+1] = z[k] + a_x * (fd[k] + fT[k])
        z_new = z + a_x_true * (fd + fT_all(k));

        % Store position
        if k < N_steps
            z_hist(k+1) = z_new;
        end
        z = z_new;
    end

    % Steady-state analysis (last 2 seconds)
    ss_idx = steady_start:N_steps;
    sigma2_dxr_meas(idx) = mean(sigma2_dzr_ts(ss_idx));
    axm_meas(idx)        = mean(axm_ts(ss_idx));
    axm_time_all{idx}    = axm_ts;

    fprintf('  C(lambda_c) = %.3f\n', C_lc);
    fprintf('  sigma2_dxr: theory = %.4e um^2, measured = %.4e um^2\n', ...
        sigma2_dxr_theory(idx), sigma2_dxr_meas(idx));
    fprintf('  a_xm: theory = %.5f um/pN, measured = %.5f um/pN\n', ...
        a_x_true, axm_meas(idx));
    fprintf('  Relative error: %.2f%%\n', ...
        abs(axm_meas(idx) - a_x_true) / a_x_true * 100);
end

%% ===== Figure 1: sigma^2_dxr vs lambda_c (Verify Eq. 12) =====
t_vec = (0:N_steps-1)' * Ts;

figure('Name', 'Eq.12 Verification: sigma^2_dxr vs lambda_c');
bar_data = [sigma2_dxr_theory; sigma2_dxr_meas]';
b = bar(1:n_lc, bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.9 0.3 0.2];
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lambda_c_list, 'UniformOutput', false));
xlabel('lambda_c');
ylabel('sigma^2_{dxr} (um^2)');
title('Eq. 12 Verification: sigma^2_{dxr} vs lambda_c');
legend('Theory (Eq. 12)', 'Measured (simulation)', 'Location', 'northwest');
grid on;

%% ===== Figure 2: a_xm time series (Verify Eq. 13 convergence) =====
figure('Name', 'Eq.13 Verification: a_xm convergence');
colors = {'b', 'r', [0 0.6 0], [0.8 0.4 0]};
hold on;
for idx = 1:n_lc
    plot(t_vec, axm_time_all{idx}, 'Color', colors{idx}, ...
        'DisplayName', sprintf('lambda_c=%.1f', lambda_c_list(idx)));
end
yline(a_x_true, 'k--', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('a_x true = %.5f', a_x_true));
hold off;
xlabel('Time (s)');
ylabel('a_{xm} (um/pN)');
title('Eq. 13 Verification: a_{xm} convergence for different lambda_c');
legend('Location', 'best');
grid on;
ylim([0, a_x_true * 3]);

%% ===== Figure 3: Steady-state a_xm bar chart + error =====
figure('Name', 'Eq.13 Summary: steady-state a_xm');

subplot(2,1,1);
bar_axm = [repmat(a_x_true, 1, n_lc); axm_meas]';
b2 = bar(1:n_lc, bar_axm);
b2(1).FaceColor = [0.2 0.4 0.8];
b2(2).FaceColor = [0.9 0.3 0.2];
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lambda_c_list, 'UniformOutput', false));
xlabel('lambda_c');
ylabel('a_{xm} (um/pN)');
title('Steady-state a_{xm} vs theory');
legend('a_x true', 'a_{xm} measured', 'Location', 'best');
grid on;

subplot(2,1,2);
rel_err = abs(axm_meas - a_x_true) / a_x_true * 100;
bar(1:n_lc, rel_err, 'FaceColor', [0.9 0.6 0.1]);
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), lambda_c_list, 'UniformOutput', false));
xlabel('lambda_c');
ylabel('Relative Error (%)');
title('Relative error of a_{xm}');
yline(5, 'r--', 'LineWidth', 1.5);
grid on;

%% ===== Summary Table =====
fprintf('\n========== VERIFICATION SUMMARY ==========\n');
fprintf('Theoretical a_x = Ts/gammaN = %.5f um/pN\n\n', a_x_true);
fprintf('%-10s %-12s %-18s %-18s %-15s %-12s\n', ...
    'lambda_c', 'C(lc)', 'sig2_theory', 'sig2_meas', 'axm_meas', 'Error(%)');
fprintf('%s\n', repmat('-', 1, 85));
for idx = 1:n_lc
    lc = lambda_c_list(idx);
    C_lc = 2 + 1 / (1 - lc^2);
    fprintf('%-10.1f %-12.3f %-18.4e %-18.4e %-15.5f %-12.2f\n', ...
        lc, C_lc, sigma2_dxr_theory(idx), sigma2_dxr_meas(idx), ...
        axm_meas(idx), rel_err(idx));
end
fprintf('%s\n', repmat('=', 1, 85));

% Pass/Fail check
if all(rel_err < 5)
    fprintf('RESULT: PASS - All relative errors < 5%%\n');
else
    fprintf('RESULT: PARTIAL - Some relative errors >= 5%%\n');
    fprintf('  Failed lambda_c values: ');
    fprintf('%.1f ', lambda_c_list(rel_err >= 5));
    fprintf('\n');
end
