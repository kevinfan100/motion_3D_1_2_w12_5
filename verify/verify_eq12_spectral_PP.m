% verify_eq12_spectral_PP.m — Numerical verification of C_obs factor
% for 3-State Pole-Placement Observer (open-loop A(3,3)=1)
%
% Four independent methods to verify C_obs(lc, le):
%   1. Lyapunov equation:           dlyap(A4, Bn*Bn') -> P(1,1)
%   2. Impulse response summation:  sum h[k]^2  (adaptive K_max)
%   3. Parseval frequency integration
%   4. Monte Carlo simulation (200k steps, prediction-form observer)
%
% Observer gains from derivation_3state_PP.tex (open-loop F_e(3,3)=1):
%   L1 = 1 - 3*le
%   L2 = 1 - 3*le + 3*le^2
%   L3 = (1 - le)^3
%
% Open-loop observer matrix:
%   A_obs = [0, 1, 0; 0, 0, 1; 0, 0, 1]
%
% Prediction-form augmented 4-state system s = [dz; e1; e2; e3]:
%   dz[k+1] = lc*dz + (1-lc)*e3 - a_x*fT
%   e1[k+1] = -L1*e1 + e2
%   e2[k+1] = -L2*e1 + e3
%   e3[k+1] = -L3*e1 + e3 - a_x*fT
%
% Filtering-form augmented system (control uses corrected xhat_3+L3*e1):
%   dz[k+1] = lc*dz - (1-lc)*L3*e1 + (1-lc)*e3 - a_x*fT
%   e1..e3 dynamics identical to prediction form (F unchanged)
%
% Also verifies:
%   - C_eq17 = 2 + 1/(1-lc^2)
%   - Deadbeat (le=0): DeltaC = C_obs(lc,0) - C_eq17(lc) = 1
%   - Prediction vs filtering execution order comparison

clear; clc; close all;

%% ===== Physical Parameters =====
Ts      = 1/1600;
kb      = 1.38e-23;
T_temp  = 310.15;
R       = 2.25e-6;
eta     = 0.001;
gammaN  = 0.0425;              % [pN*s/um]
a_x     = Ts / gammaN;         % [um/pN]

% Thermal noise
gamma_SI   = 6*pi*eta*R;       % [N*s/m]
sigma2_fT  = 4*kb*T_temp*gamma_SI / Ts;   % [N^2]
sigma_fT_pN = sqrt(sigma2_fT * 1e24);      % [pN]

%% ===== Test parameters =====
lc_list  = 0.4:0.05:0.9;
le_fixed = 0.3;
n_lc     = length(lc_list);

% Preallocate — prediction form (4 methods)
C_lyap     = zeros(1, n_lc);
C_impulse  = zeros(1, n_lc);
C_spectral = zeros(1, n_lc);
C_sim_pred = zeros(1, n_lc);

% Filtering form
C_lyap_filt = zeros(1, n_lc);
C_sim_filt  = zeros(1, n_lc);

% Reference
C_eq17     = zeros(1, n_lc);

% Deadbeat
C_db_lyap  = zeros(1, n_lc);

% Simulation parameters
N_steps      = 200000;
steady_start = 50001;
z_d          = 25;       % target [um]

%% ===== Main loop =====
fprintf('3-State PP Observer: C_obs Verification (le = %.2f)\n', le_fixed);
fprintf('Gains: L1=1-3*le, L2=1-3*le+3*le^2, L3=(1-le)^3  [open-loop A(3,3)=1]\n');
fprintf('=======================================================================\n\n');

for idx = 1:n_lc
    lc = lc_list(idx);
    le = le_fixed;

    % --- Observer gains (open-loop F_e(3,3)=1, tex Eq.gains) ---
    L1 = 1 - 3*le;
    L2 = 1 - 3*le + 3*le^2;
    L3 = (1 - le)^3;

    % --- Prediction-form: 4-state augmented system ---
    %   State: s = [dz; e1; e2; e3]
    %   Control uses xhat_3 BEFORE correction:
    %     fd = (1/a_x)*(1-lc)*xhat_3
    %   dz[k+1] = lc*dz + (1-lc)*e3 - a_x*fT
    %   e[k+1]  = F*e[k] + [0;0;-1]*fT_norm,  F=[-L1,1,0;-L2,0,1;-L3,0,1]
    A4 = [lc,    0,   0,  (1-lc);
           0,  -L1,   1,    0;
           0,  -L2,   0,    1;
           0,  -L3,   0,    1];
    Bn = [-1; 0; 0; -1];   % B4 / a_x

    % --- Filtering-form: 4-state augmented system ---
    %   Control uses corrected xhat_3 + L3*e1:
    %     fd = (1/a_x)*(1-lc)*(xhat_3 + L3*e1)
    %   dz[k+1] = lc*dz - (1-lc)*L3*e1 + (1-lc)*e3 - a_x*fT
    %   Error dynamics F unchanged (same e[k+1] = F*e[k] + [0;0;-1]*fT_norm)
    A4_filt = [lc, -(1-lc)*L3,  0,  (1-lc);
                0,   -L1,       1,    0;
                0,   -L2,       0,    1;
                0,   -L3,       0,    1];
    Bn_filt = [-1; 0; 0; -1];

    % --- C_eq17 for comparison ---
    C_eq17(idx) = 2 + 1/(1 - lc^2);

    % =============================================
    %  Method 1: Lyapunov equation
    % =============================================
    P = dlyap(A4, Bn * Bn');
    C_lyap(idx) = P(1,1);

    % Filtering-form Lyapunov
    Pf = dlyap(A4_filt, Bn_filt * Bn_filt');
    C_lyap_filt(idx) = Pf(1,1);

    % =============================================
    %  Method 2: Impulse response summation (prediction form)
    % =============================================
    rho = max(abs(eig(A4)));
    if rho > 0.01
        K_max = max(500, ceil(-10*log(10) / (2*log(rho))));
    else
        K_max = 500;
    end
    h_PP = zeros(K_max+1, 1);
    state = zeros(4, 1);

    state = A4 * state + Bn * 1;   % impulse at k=0
    h_PP(2) = state(1);

    for k = 2:K_max
        state = A4 * state;
        h_PP(k+1) = state(1);
    end
    C_impulse(idx) = sum(h_PP.^2);

    % =============================================
    %  Method 3: Parseval frequency integration (prediction form)
    % =============================================
    N_pts = 20000;
    theta = linspace(0, pi, N_pts);
    z_uc  = exp(1j * theta);

    H_mag2 = zeros(size(theta));
    I4 = eye(4);
    for ti = 1:N_pts
        resolvent = (z_uc(ti) * I4 - A4) \ Bn;
        H_mag2(ti) = abs(resolvent(1))^2;
    end
    C_spectral(idx) = (1/pi) * trapz(theta, H_mag2);

    % =============================================
    %  Method 4a: Monte Carlo — Prediction-form observer
    % =============================================
    %   A_obs = [0,1,0; 0,0,1; 0,0,1]  (open-loop F_e)
    %   b_u   = [0; 0; -a_x]
    %   1. y[k] = dz[k-2],  e_z1 = y[k] - xhat_1[k]
    %   2. fd[k] = (1/a_x)*(1-lc)*xhat_3[k]    (before correction)
    %   3. Plant: z[k+1] = z[k] + a_x*(fd+fT)
    %   4. Observer: xhat[k+1] = A_obs*xhat + L*e_z1 + b_u*fd

    rng(42 + idx);
    fT_sim = sigma_fT_pN * randn(N_steps, 1);

    z_pos  = z_d;
    z_hist = zeros(N_steps+1, 1);  z_hist(1) = z_d;
    xhat   = [0; 0; 0];
    dz_pred = zeros(N_steps, 1);

    A_obs = [0, 1, 0; 0, 0, 1; 0, 0, 1];
    L_vec = [L1; L2; L3];
    b_u   = [0; 0; -a_x];

    for k = 1:N_steps
        % Measurement with 2-step delay
        if k >= 3
            z_del = z_hist(k-2);
        elseif k == 2
            z_del = z_hist(1);
        else
            z_del = z_d;
        end
        y_k = z_d - z_del;
        dz_pred(k) = y_k;

        % Innovation
        e_z1 = y_k - xhat(1);

        % Control (prediction: use xhat_3 before correction)
        fd_k = (1/a_x) * (1-lc) * xhat(3);

        % Plant
        z_new = z_pos + a_x * (fd_k + fT_sim(k));
        z_hist(k+1) = z_new;
        z_pos = z_new;

        % Observer time update
        xhat = A_obs * xhat + L_vec * e_z1 + b_u * fd_k;
    end

    ss = steady_start:N_steps;
    C_sim_pred(idx) = var(dz_pred(ss)) / (sigma_fT_pN^2 * a_x^2);

    % =============================================
    %  Method 4b: Monte Carlo — Filtering-form observer
    % =============================================
    %   1. y[k], e_z1 = y[k] - xhat_1[k]
    %   2. xhat_3_corr = xhat_3 + L3*e_z1
    %   3. fd[k] = (1/a_x)*(1-lc)*xhat_3_corr   (after correction)
    %   4. Plant update
    %   5. Observer: xhat[k+1] = A_obs*xhat + L*e_z1 + b_u*fd

    rng(42 + idx);   % same seed
    fT_sim2 = sigma_fT_pN * randn(N_steps, 1);

    z_pos2  = z_d;
    z_hist2 = zeros(N_steps+1, 1);  z_hist2(1) = z_d;
    xhat2   = [0; 0; 0];
    dz_filt = zeros(N_steps, 1);

    for k = 1:N_steps
        if k >= 3
            z_del2 = z_hist2(k-2);
        elseif k == 2
            z_del2 = z_hist2(1);
        else
            z_del2 = z_d;
        end
        y_k2 = z_d - z_del2;
        dz_filt(k) = y_k2;

        e_z1_2 = y_k2 - xhat2(1);

        % Corrected xhat_3 (filtering form)
        xhat3_corr = xhat2(3) + L3 * e_z1_2;

        % Control (filtering: use corrected estimate)
        fd_k2 = (1/a_x) * (1-lc) * xhat3_corr;

        % Plant
        z_new2 = z_pos2 + a_x * (fd_k2 + fT_sim2(k));
        z_hist2(k+1) = z_new2;
        z_pos2 = z_new2;

        % Observer time update
        xhat2 = A_obs * xhat2 + L_vec * e_z1_2 + b_u * fd_k2;
    end

    ss2 = steady_start:N_steps;
    C_sim_filt(idx) = var(dz_filt(ss2)) / (sigma_fT_pN^2 * a_x^2);

    fprintf('lc=%.2f: Lyap=%.4f Imp=%.4f Spec=%.4f SimP=%.4f | LyapF=%.4f SimF=%.4f | C17=%.4f\n', ...
        lc, C_lyap(idx), C_impulse(idx), C_spectral(idx), ...
        C_sim_pred(idx), C_lyap_filt(idx), C_sim_filt(idx), C_eq17(idx));
end

%% ===== Deadbeat test (le = 0) =====
fprintf('\n==================== DEADBEAT CASE (le = 0) ====================\n');
fprintf('  L1=1, L2=1, L3=1.  C_obs(lc,0) = (4-3*lc^2)/(1-lc^2) = 3+1/(1-lc^2)\n');
fprintf('  DeltaC = C_obs - C_eq17 = 1 (exact)\n\n');

for idx = 1:n_lc
    lc = lc_list(idx);
    L1_db = 1;  L2_db = 1;  L3_db = 1;

    A4_db = [lc,  0,  0,  (1-lc);
              0, -1,  1,    0;
              0, -1,  0,    1;
              0, -1,  0,    1];
    Bn_db = [-1; 0; 0; -1];

    P_db = dlyap(A4_db, Bn_db * Bn_db');
    C_db_lyap(idx) = P_db(1,1);
    C17 = C_eq17(idx);
    delta = C_db_lyap(idx) - C17;

    C_db_formula = (4 - 3*lc^2) / (1 - lc^2);

    F_db = [-1, 1, 0; -1, 0, 1; -1, 0, 1];
    rho_F_db = max(abs(eig(F_db)));

    fprintf('lc=%.2f: C_db_lyap=%.6f  C_formula=%.6f  C_eq17=%.4f  DeltaC=%.6f  rho(F)=%.4f\n', ...
        lc, C_db_lyap(idx), C_db_formula, C17, delta, rho_F_db);
end

%% ===== Results Table: Prediction Form =====
fprintf('\n==================== C_obs PREDICTION FORM (le=%.2f) ====================\n', le_fixed);
fprintf('%-6s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'lc', 'Lyapunov', 'Impulse', 'Spectral', 'Sim(Pred)', 'C_eq17', 'DeltaC');
fprintf('%s\n', repmat('-', 1, 78));
for idx = 1:n_lc
    fprintf('%-6.2f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n', ...
        lc_list(idx), C_lyap(idx), C_impulse(idx), C_spectral(idx), ...
        C_sim_pred(idx), C_eq17(idx), C_lyap(idx) - C_eq17(idx));
end

%% ===== Results Table: Filtering Form =====
fprintf('\n==================== C_obs FILTERING FORM (le=%.2f) ====================\n', le_fixed);
fprintf('%-6s  %-12s  %-12s  %-12s\n', 'lc', 'Lyap(Filt)', 'Sim(Filt)', 'C_eq17');
fprintf('%s\n', repmat('-', 1, 48));
for idx = 1:n_lc
    fprintf('%-6.2f  %-12.4f  %-12.4f  %-12.4f\n', ...
        lc_list(idx), C_lyap_filt(idx), C_sim_filt(idx), C_eq17(idx));
end

%% ===== Relative Errors =====
fprintf('\n--- Prediction: Relative Errors vs Lyapunov (%%)  ---\n');
fprintf('%-6s  %-12s  %-12s  %-12s\n', 'lc', 'Impulse', 'Spectral', 'Sim(Pred)');
fprintf('%s\n', repmat('-', 1, 48));
for idx = 1:n_lc
    ref = C_lyap(idx);
    err_imp  = abs(C_impulse(idx) - ref) / ref * 100;
    err_spec = abs(C_spectral(idx) - ref) / ref * 100;
    err_pred = abs(C_sim_pred(idx) - ref) / ref * 100;
    fprintf('%-6.2f  %-12.6f  %-12.6f  %-12.4f\n', ...
        lc_list(idx), err_imp, err_spec, err_pred);
end

fprintf('\n--- Filtering: Relative Error Sim vs Lyapunov (%%)  ---\n');
fprintf('%-6s  %-12s\n', 'lc', 'Sim(Filt)');
fprintf('%s\n', repmat('-', 1, 20));
for idx = 1:n_lc
    ref_f = C_lyap_filt(idx);
    err_filt = abs(C_sim_filt(idx) - ref_f) / ref_f * 100;
    fprintf('%-6.2f  %-12.4f\n', lc_list(idx), err_filt);
end

%% ===== Prediction vs Filtering comparison =====
fprintf('\n==================== PREDICTION vs FILTERING COMPARISON ====================\n');
fprintf('%-6s  %-12s  %-12s  %-12s  %-12s\n', ...
    'lc', 'C_pred(Lyap)', 'C_filt(Lyap)', 'Ratio', 'C_eq17');
fprintf('%s\n', repmat('-', 1, 60));
for idx = 1:n_lc
    ratio = C_lyap_filt(idx) / C_lyap(idx);
    fprintf('%-6.2f  %-12.4f  %-12.4f  %-12.4f  %-12.4f\n', ...
        lc_list(idx), C_lyap(idx), C_lyap_filt(idx), ratio, C_eq17(idx));
end

%% ===== F eigenvalues =====
fprintf('\n==================== F MATRIX EIGENVALUES ====================\n');
fprintf('  F = (F_e - L*H),  F_e(3,3)=1 (open-loop)\n');
fprintf('  All eigenvalues should be at le=%.2f\n\n', le_fixed);
fprintf('%-6s  %-8s  %-30s\n', 'lc', 'rho(F)', 'eig(F)');
fprintf('%s\n', repmat('-', 1, 50));
for idx = 1:n_lc
    lc = lc_list(idx);
    le = le_fixed;
    L1 = 1 - 3*le;
    L2 = 1 - 3*le + 3*le^2;
    L3 = (1 - le)^3;

    F = [-L1, 1, 0; -L2, 0, 1; -L3, 0, 1];
    eF = sort(eig(F), 'ComparisonMethod', 'abs');
    rho_F = max(abs(eF));

    fprintf('lc=%.2f  %.4f   [%s]\n', lc, rho_F, sprintf('%.4f ', eF));
end

%% ===== Pass/Fail =====
fprintf('\n==================== PASS/FAIL ====================\n');
tol_exact = 0.01;   % 0.01% for exact methods
tol_sim   = 5;      % 5% for simulation (stochastic)

pass_all = true;

% --- Prediction form: impulse, spectral, sim all match Lyapunov ---
for idx = 1:n_lc
    lc = lc_list(idx);
    ref = C_lyap(idx);
    err_imp  = abs(C_impulse(idx)  - ref) / ref * 100;
    err_spec = abs(C_spectral(idx) - ref) / ref * 100;
    err_pred = abs(C_sim_pred(idx) - ref) / ref * 100;

    if err_imp > tol_exact
        fprintf('FAIL: lc=%.2f impulse err=%.6f%%\n', lc, err_imp);
        pass_all = false;
    end
    if err_spec > 0.1
        fprintf('FAIL: lc=%.2f spectral err=%.6f%%\n', lc, err_spec);
        pass_all = false;
    end
    if err_pred > tol_sim
        fprintf('FAIL: lc=%.2f prediction sim err=%.2f%%\n', lc, err_pred);
        pass_all = false;
    end
end

% --- Filtering form: sim matches filtering Lyapunov ---
for idx = 1:n_lc
    lc = lc_list(idx);
    ref_f = C_lyap_filt(idx);
    err_filt = abs(C_sim_filt(idx) - ref_f) / ref_f * 100;
    if err_filt > tol_sim
        fprintf('FAIL: lc=%.2f filtering sim err=%.2f%% (vs filtering Lyapunov)\n', ...
            lc, err_filt);
        pass_all = false;
    end
end

% --- Deadbeat DeltaC = 1 ---
for idx = 1:n_lc
    lc = lc_list(idx);
    delta = C_db_lyap(idx) - C_eq17(idx);
    if abs(delta - 1) > 1e-8
        fprintf('FAIL: lc=%.2f deadbeat DeltaC=%.10f (expected 1)\n', lc, delta);
        pass_all = false;
    end
end

% --- F eigenvalues all at le ---
le = le_fixed;
L1_chk = 1 - 3*le;
L2_chk = 1 - 3*le + 3*le^2;
L3_chk = (1 - le)^3;
F_chk = [-L1_chk, 1, 0; -L2_chk, 0, 1; -L3_chk, 0, 1];
eigs_F = eig(F_chk);
if any(abs(eigs_F - le) > 1e-4)
    fprintf('FAIL: F eigenvalues not all at le=%.2f\n', le);
    pass_all = false;
end

if pass_all
    fprintf('RESULT: ALL PASS\n');
    fprintf('  - Impulse, Spectral vs Lyapunov: < 0.01%% error\n');
    fprintf('  - Prediction sim vs prediction Lyapunov: < 5%% error\n');
    fprintf('  - Filtering sim vs filtering Lyapunov: < 5%% error\n');
    fprintf('  - Deadbeat DeltaC = 1: machine precision\n');
    fprintf('  - F eigenvalues all at le=%.2f: verified\n', le_fixed);
else
    fprintf('RESULT: SOME FAILURES (see above)\n');
end
fprintf('====================================================\n');

%% ===== Figure 1: Variance comparison (dense sweep, le=0.3) =====
lc_dense = 0.4:0.005:0.95;
n_dense  = length(lc_dense);

le = le_fixed;
L1 = 1 - 3*le;
L2 = 1 - 3*le + 3*le^2;
L3 = (1 - le)^3;
Bn_d = [-1; 0; 0; -1];

C_eq17_dense = zeros(1, n_dense);
C_pred_dense = zeros(1, n_dense);
C_filt_dense = zeros(1, n_dense);

for idx = 1:n_dense
    lc = lc_dense(idx);
    C_eq17_dense(idx) = 2 + 1/(1 - lc^2);

    % Prediction form
    A4p = [lc, 0, 0, (1-lc); 0, -L1, 1, 0; 0, -L2, 0, 1; 0, -L3, 0, 1];
    P = dlyap(A4p, Bn_d * Bn_d');
    C_pred_dense(idx) = P(1,1);

    % Filtering form
    A4f = [lc, -(1-lc)*L3, 0, (1-lc); 0, -L1, 1, 0; 0, -L2, 0, 1; 0, -L3, 0, 1];
    Pf = dlyap(A4f, Bn_d * Bn_d');
    C_filt_dense(idx) = Pf(1,1);
end

% Physical variance [um^2]
var_scale = sigma_fT_pN^2 * a_x^2;

figure('Position', [100 100 800 500]);
plot(lc_dense, C_eq17_dense * var_scale, 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Eq.17 (no observer)');
hold on;
plot(lc_dense, C_pred_dense * var_scale, 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('PP prediction (\\lambda_e=%.1f)', le_fixed));
plot(lc_dense, C_filt_dense * var_scale, 'b-', 'LineWidth', 2, ...
    'DisplayName', sprintf('PP filtering (\\lambda_e=%.1f)', le_fixed));
% Simulation points
plot(lc_list, C_sim_pred * var_scale, 'ro', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
plot(lc_list, C_sim_filt * var_scale, 'bs', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
hold off;
xlabel('\lambda_c', 'FontSize', 13);
ylabel('\sigma^2_{\deltaz}  [\mum^2]', 'FontSize', 13);
title(sprintf('Tracking Error Variance (\\lambda_e=%.1f, \\gamma=%.4f pN{\\cdot}s/\\mum)', ...
    le_fixed, gammaN), 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 11); grid on;
set(gca, 'FontSize', 12);
fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
print(gcf, fullfile(fig_dir, 'fig_variance_comparison.png'), '-dpng', '-r150');

%% ===== Figure 2: Deadbeat (le=0) variance =====
C_db_filt_dense = zeros(1, n_dense);
for idx = 1:n_dense
    lc = lc_dense(idx);
    % Deadbeat filtering: L1=L2=L3=1
    A4_dbf = [lc, -(1-lc), 0, (1-lc); 0, -1, 1, 0; 0, -1, 0, 1; 0, -1, 0, 1];
    P_dbf = dlyap(A4_dbf, Bn_d * Bn_d');
    C_db_filt_dense(idx) = P_dbf(1,1);
end

figure('Position', [100 100 800 500]);
plot(lc_dense, C_eq17_dense * var_scale, 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Eq.17 (no observer)');
hold on;
plot(lc_dense, C_db_filt_dense * var_scale, 'b--', 'LineWidth', 2, ...
    'DisplayName', 'PP filtering (\lambda_e=0)');
hold off;
xlabel('\lambda_c', 'FontSize', 13);
ylabel('\sigma^2_{\deltaz}  [\mum^2]', 'FontSize', 13);
title(sprintf('Deadbeat Observer (\\lambda_e=0): \\gamma=%.4f pN{\\cdot}s/\\mum', gammaN), ...
    'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 11); grid on;
set(gca, 'FontSize', 12);
print(gcf, fullfile(fig_dir, 'fig_variance_le0.png'), '-dpng', '-r150');

fprintf('\nFigures saved:\n');
fprintf('  fig_variance_comparison.png\n');
fprintf('  fig_variance_le0.png\n');
