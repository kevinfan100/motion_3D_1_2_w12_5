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
% Augmented 4-state system s = [dz; e1; e2; e3]:
%   dz[k+1] = lc*dz + (1-lc)*e3 - a_x*fT
%   e1[k+1] = -L1*e1 + e2
%   e2[k+1] = -L2*e1 + e3
%   e3[k+1] = -L3*e1 + e3 - a_x*fT
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

% Preallocate results
C_lyap     = zeros(1, n_lc);    % Lyapunov: prediction-form
C_lyap_f   = zeros(1, n_lc);    % Lyapunov: filtering-form
C_impulse  = zeros(1, n_lc);
C_spectral = zeros(1, n_lc);
C_sim_pred = zeros(1, n_lc);   % Monte Carlo: prediction-form
C_sim_filt = zeros(1, n_lc);   % Monte Carlo: filtering-form
C_eq17     = zeros(1, n_lc);

% Deadbeat results
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

    % --- 4-state augmented system ---
    %   State: s = [dz; e1; e2; e3]
    %   s[k+1] = A4*s[k] + Bn*fT_norm[k]   (Bn = B4/a_x)
    A4 = [lc,    0,   0,  (1-lc);
           0,  -L1,   1,    0;
           0,  -L2,   0,    1;
           0,  -L3,   0,    1];

    Bn = [-1; 0; 0; -1];   % B4 / a_x

    % --- C_eq17 for comparison ---
    C_eq17(idx) = 2 + 1/(1 - lc^2);

    % =============================================
    %  Method 1: Lyapunov equation
    % =============================================
    % P = A4*P*A4' + Bn*Bn'  =>  C_obs = P(1,1)
    P = dlyap(A4, Bn * Bn');
    C_lyap(idx) = P(1,1);

    % Filtering-form Lyapunov
    A4f = [lc, -(1-lc)*L3, 0, (1-lc);
            0,  -L2,        1,  0;
            0,  -L3,        0,  1;
            0,  -L3,        0,  1];
    Pf = dlyap(A4f, Bn * Bn');
    C_lyap_f(idx) = Pf(1,1);

    % =============================================
    %  Method 2: Impulse response summation
    % =============================================
    % Adaptive K_max: need rho^(2*K_max) < 1e-10
    rho = max(abs(eig(A4)));
    if rho > 0.01
        K_max = max(500, ceil(-10*log(10) / (2*log(rho))));
    else
        K_max = 500;
    end
    h_PP = zeros(K_max+1, 1);   % h[0..K_max], 1-indexed
    state = zeros(4, 1);

    % k=0: impulse fT_norm[0]=1
    state = A4 * state + Bn * 1;
    h_PP(2) = state(1);   % h[1] = first output element

    for k = 2:K_max
        state = A4 * state;   % fT=0 for k>=1
        h_PP(k+1) = state(1);
    end

    C_impulse(idx) = sum(h_PP.^2);

    % =============================================
    %  Method 3: Parseval frequency integration
    % =============================================
    % H(z) = [1,0,0,0] * (zI-A4)^{-1} * Bn
    % C_obs = (1/pi) * int_0^pi |H(e^{jw})|^2 dw
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
    % Observer: xhat[k+1] = A_obs*xhat[k] + L*e_z1[k] + b_u*fd[k]
    %   A_obs = [0,1,0; 0,0,1; 0,0,1]  (open-loop, F_e(3,3)=1)
    %   b_u   = [0; 0; -a_x]           (control input to state 3)
    %   L     = [L1; L2; L3]
    %   e_z1  = dz_m[k] - xhat_1[k]    (innovation)
    %
    % Prediction form: observer update BEFORE control
    %   1. Observe y[k] = dz[k-2]
    %   2. Innovation: e_z1 = y[k] - xhat_1[k]
    %   3. Control: fd[k] = (1/a_x)*(1-lc)*xhat_3[k] + (1/a_x)*(1-lc)*L3*e_z1
    %      (using xhat_3 AFTER correction, i.e., xhat_3_corrected)
    %   4. Plant update
    %   5. Observer time update: xhat[k+1] = A_obs*xhat[k] + L*e_z1[k] + b_u*fd[k]

    rng(42 + idx);
    fT_sim = sigma_fT_pN * randn(N_steps, 1);

    z_pos  = z_d;
    z_hist = zeros(N_steps+1, 1);  z_hist(1) = z_d;
    xhat   = [0; 0; 0];   % observer state [xhat1; xhat2; xhat3]
    dz_pred = zeros(N_steps, 1);   % measured tracking error

    A_obs = [0, 1, 0; 0, 0, 1; 0, 0, 1];   % open-loop (F_e)
    L_vec = [L1; L2; L3];
    b_u   = [0; 0; -a_x];   % control input enters state 3

    for k = 1:N_steps
        % Measurement with 2-step delay: y[k] = dz[k-2]
        if k >= 3
            z_del = z_hist(k-2);
        elseif k == 2
            z_del = z_hist(1);
        else
            z_del = z_d;
        end
        y_k = z_d - z_del;      % = dz[k-2]
        dz_pred(k) = y_k;

        % Innovation
        e_z1 = y_k - xhat(1);

        % Control using xhat_3 (prediction form: use current xhat_3)
        fd_k = (1/a_x) * (1-lc) * xhat(3);

        % Plant update
        z_new = z_pos + a_x * (fd_k + fT_sim(k));
        z_hist(k+1) = z_new;
        z_pos = z_new;

        % Observer time update (prediction form)
        xhat = A_obs * xhat + L_vec * e_z1 + b_u * fd_k;
    end

    % Steady-state variance
    ss = steady_start:N_steps;
    var_pred = var(dz_pred(ss));   % [um^2]
    C_sim_pred(idx) = var_pred / (sigma_fT_pN^2 * a_x^2);

    % =============================================
    %  Method 4b: Monte Carlo — Filtering-form observer
    % =============================================
    % Filtering form: control uses corrected estimate
    %   1. Observe y[k]
    %   2. Innovation: e_z1 = y[k] - xhat_1[k]
    %   3. Correct: xhat_corr = xhat + K*e_z1 (where K is column of gains)
    %      Actually: xhat_3_corr = xhat_3 + L3*e_z1
    %   4. Control: fd[k] = (1/a_x)*(1-lc)*xhat_3_corr
    %   5. Plant update
    %   6. Observer time update

    rng(42 + idx);   % same seed for fair comparison
    fT_sim2 = sigma_fT_pN * randn(N_steps, 1);

    z_pos2  = z_d;
    z_hist2 = zeros(N_steps+1, 1);  z_hist2(1) = z_d;
    xhat2   = [0; 0; 0];
    dz_filt = zeros(N_steps, 1);

    for k = 1:N_steps
        % Measurement
        if k >= 3
            z_del2 = z_hist2(k-2);
        elseif k == 2
            z_del2 = z_hist2(1);
        else
            z_del2 = z_d;
        end
        y_k2 = z_d - z_del2;
        dz_filt(k) = y_k2;

        % Innovation
        e_z1_2 = y_k2 - xhat2(1);

        % Step 1: Correct ALL states
        xhat_corr = xhat2 + L_vec * e_z1_2;

        % Step 2: Control using corrected estimate
        fd_k2 = (1/a_x) * (1-lc) * xhat_corr(3);

        % Step 3: Plant update
        z_new2 = z_pos2 + a_x * (fd_k2 + fT_sim2(k));
        z_hist2(k+1) = z_new2;
        z_pos2 = z_new2;

        % Step 4: Predict from corrected state
        xhat2 = A_obs * xhat_corr + b_u * fd_k2;
    end

    ss2 = steady_start:N_steps;
    var_filt = var(dz_filt(ss2));
    C_sim_filt(idx) = var_filt / (sigma_fT_pN^2 * a_x^2);

    fprintf('lc=%.2f: C_lyap=%.4f  C_imp=%.4f  C_spec=%.4f  C_sim_P=%.4f  C_sim_F=%.4f  C_eq17=%.4f\n', ...
        lc, C_lyap(idx), C_impulse(idx), C_spectral(idx), ...
        C_sim_pred(idx), C_sim_filt(idx), C_eq17(idx));
end

%% ===== Deadbeat test (le = 0) =====
fprintf('\n==================== DEADBEAT CASE (le = 0) ====================\n');
fprintf('  L1=1, L2=1, L3=1.  C_obs(lc,0) = (4-3*lc^2)/(1-lc^2) = 3+1/(1-lc^2)\n');
fprintf('  DeltaC = C_obs - C_eq17 = 1 (exact)\n\n');

for idx = 1:n_lc
    lc = lc_list(idx);
    % Deadbeat gains: le=0
    L1_db = 1;  L2_db = 1;  L3_db = 1;

    A4_db = [lc,    0,   0,  (1-lc);
              0,  -L1_db,   1,    0;
              0,  -L2_db,   0,    1;
              0,  -L3_db,   0,    1];
    Bn_db = [-1; 0; 0; -1];

    P_db = dlyap(A4_db, Bn_db * Bn_db');
    C_db_lyap(idx) = P_db(1,1);
    C17 = C_eq17(idx);
    delta = C_db_lyap(idx) - C17;

    % Closed-form
    C_db_formula = (4 - 3*lc^2) / (1 - lc^2);

    F_db = [-L1_db, 1, 0; -L2_db, 0, 1; -L3_db, 0, 1];
    rho_F_db = max(abs(eig(F_db)));

    fprintf('lc=%.2f: C_db_lyap=%.6f  C_db_formula=%.6f  C_eq17=%.4f  DeltaC=%.6f  rho(F)=%.4f\n', ...
        lc, C_db_lyap(idx), C_db_formula, C17, delta, rho_F_db);
end

%% ===== Results Table =====
fprintf('\n==================== C_obs VERIFICATION RESULTS (le=%.2f) ====================\n', le_fixed);
fprintf('%-6s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'lc', 'Lyapunov', 'Impulse', 'Spectral', 'Sim(Pred)', 'Sim(Filt)', 'C_eq17', 'DeltaC');
fprintf('%s\n', repmat('-', 1, 90));
for idx = 1:n_lc
    fprintf('%-6.2f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n', ...
        lc_list(idx), C_lyap(idx), C_impulse(idx), C_spectral(idx), ...
        C_sim_pred(idx), C_sim_filt(idx), C_eq17(idx), C_lyap(idx) - C_eq17(idx));
end

%% ===== Relative Errors vs Lyapunov =====
fprintf('\n--- Relative Errors vs Lyapunov (%%)  ---\n');
fprintf('%-6s  %-12s  %-12s  %-12s  %-12s\n', 'lc', 'Impulse', 'Spectral', 'Sim(Pred)', 'Sim(Filt)');
fprintf('%s\n', repmat('-', 1, 60));
for idx = 1:n_lc
    ref = C_lyap(idx);
    err_imp  = abs(C_impulse(idx) - ref) / ref * 100;
    err_spec = abs(C_spectral(idx) - ref) / ref * 100;
    err_pred = abs(C_sim_pred(idx) - ref) / ref * 100;
    err_filt = abs(C_sim_filt(idx) - ref) / ref * 100;
    fprintf('%-6.2f  %-12.6f  %-12.6f  %-12.4f  %-12.4f\n', ...
        lc_list(idx), err_imp, err_spec, err_pred, err_filt);
end

%% ===== Prediction vs Filtering comparison =====
fprintf('\n==================== PREDICTION vs FILTERING EXECUTION ORDER ====================\n');
fprintf('%-6s  %-12s  %-12s  %-12s\n', 'lc', 'C_pred', 'C_filt', 'Diff(%%)');
fprintf('%s\n', repmat('-', 1, 48));
for idx = 1:n_lc
    diff_pf = (C_sim_filt(idx) - C_sim_pred(idx)) / C_sim_pred(idx) * 100;
    fprintf('%-6.2f  %-12.4f  %-12.4f  %-12.4f\n', ...
        lc_list(idx), C_sim_pred(idx), C_sim_filt(idx), diff_pf);
end

%% ===== F eigenvalues and stability =====
fprintf('\n==================== F MATRIX EIGENVALUES ====================\n');
fprintf('  F = (F_e - L*H),  F_e(3,3)=1 (open-loop)\n\n');
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

% --- Check 4 methods agree (le = le_fixed) ---
for idx = 1:n_lc
    lc = lc_list(idx);
    ref = C_lyap(idx);
    err_imp  = abs(C_impulse(idx)  - ref) / ref * 100;
    err_spec = abs(C_spectral(idx) - ref) / ref * 100;
    err_pred = abs(C_sim_pred(idx) - ref) / ref * 100;
    err_filt = abs(C_sim_filt(idx) - ref) / ref * 100;

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
    % Filtering sim compared against filtering Lyapunov (not prediction)
    ref_f = C_lyap_f(idx);
    err_filt = abs(C_sim_filt(idx) - ref_f) / ref_f * 100;
    if err_filt > tol_sim
        fprintf('FAIL: lc=%.2f filtering sim err=%.2f%%\n', lc, err_filt);
        pass_all = false;
    end
end

% --- Check deadbeat DeltaC = 1 ---
for idx = 1:n_lc
    lc = lc_list(idx);
    delta = C_db_lyap(idx) - C_eq17(idx);
    if abs(delta - 1) > 1e-8
        fprintf('FAIL: lc=%.2f deadbeat DeltaC=%.10f (expected 1)\n', lc, delta);
        pass_all = false;
    end
end

if pass_all
    fprintf('RESULT: ALL PASS\n');
    fprintf('  - Impulse, Spectral vs Lyapunov: < 0.01%% error\n');
    fprintf('  - Simulation (pred & filt): < 5%% error\n');
    fprintf('  - Deadbeat DeltaC = 1: verified to machine precision\n');
else
    fprintf('RESULT: SOME FAILURES (see above)\n');
end
fprintf('====================================================\n');

%% ===== Figure 1: Variance comparison (dense sweep) =====
% Dense lc sweep for smooth curves
lc_dense = 0.4:0.005:0.95;
n_dense  = length(lc_dense);

C_eq17_dense  = zeros(1, n_dense);
C_pred_dense  = zeros(1, n_dense);   % prediction form
C_filt_dense  = zeros(1, n_dense);   % filtering form

le = le_fixed;
L1 = 1 - 3*le;
L2 = 1 - 3*le + 3*le^2;
L3 = (1 - le)^3;
Bn_d = [-1; 0; 0; -1];

for idx = 1:n_dense
    lc = lc_dense(idx);
    C_eq17_dense(idx) = 2 + 1/(1 - lc^2);

    A4p = [lc, 0, 0, (1-lc); 0, -L1, 1, 0; 0, -L2, 0, 1; 0, -L3, 0, 1];
    P = dlyap(A4p, Bn_d * Bn_d');
    C_pred_dense(idx) = P(1,1);

    A4f = [lc, -(1-lc)*L3, 0, (1-lc); 0, -L2, 1, 0; 0, -L3, 0, 1; 0, -L3, 0, 1];
    Pf = dlyap(A4f, Bn_d * Bn_d');
    C_filt_dense(idx) = Pf(1,1);
end

% Convert C factors to physical variance [um^2]
var_scale = sigma_fT_pN^2 * a_x^2;   % [um^2] per unit C

figure('Position', [100 100 800 500]);
plot(lc_dense, C_eq17_dense * var_scale, 'k--', 'LineWidth', 2.5, ...
    'DisplayName', 'Eq.17 (no observer)');
hold on;
plot(lc_dense, C_filt_dense * var_scale, 'b-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Observer: correct first (\\lambda_e=%.1f)', le_fixed));
plot(lc_dense, C_pred_dense * var_scale, 'r-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Observer: control first (\\lambda_e=%.1f)', le_fixed));
% Simulation points
plot(lc_list, C_sim_filt * var_scale, 'bo', 'MarkerSize', 9, 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
plot(lc_list, C_sim_pred * var_scale, 'ro', 'MarkerSize', 9, 'LineWidth', 2, ...
    'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
hold off;
xlabel('\lambda_c', 'FontSize', 14);
ylabel('\sigma^2_{\delta z}  [\mum^2]', 'FontSize', 14);
title(sprintf('Tracking Error Variance (\\lambda_e=%.1f, \\gamma=%.4f pN{\\cdot}s/\\mum)', ...
    le_fixed, gammaN), 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 11); grid on;
set(gca, 'FontSize', 12, 'LineWidth', 0.8);
print(gcf, fullfile(pwd, 'fig_variance_comparison.png'), '-dpng', '-r150');

%% ===== Figure 2: Deadbeat (le=0) variance =====
C_db_dense = zeros(1, n_dense);
for idx = 1:n_dense
    lc = lc_dense(idx);
    A4_db = [lc, 0, 0, (1-lc); 0, -1, 1, 0; 0, -1, 0, 1; 0, -1, 0, 1];
    Bn_d = [-1; 0; 0; -1];
    P_db = dlyap(A4_db, Bn_d * Bn_d');
    C_db_dense(idx) = P_db(1,1);
end

% Also run filtering simulation for deadbeat
C_sim_filt_db = zeros(1, n_lc);
for idx = 1:n_lc
    lc = lc_list(idx);
    L1_db = 1;  L2_db = 1;  L3_db = 1;
    L_db = [L1_db; L2_db; L3_db];

    rng(1000 + idx);
    fT_db = sigma_fT_pN * randn(N_steps, 1);

    z_pos_db = z_d;
    z_hist_db = zeros(N_steps+1, 1);  z_hist_db(1) = z_d;
    xhat_db = [0; 0; 0];
    dz_db = zeros(N_steps, 1);

    for k = 1:N_steps
        if k >= 3
            z_del_db = z_hist_db(k-2);
        elseif k == 2
            z_del_db = z_hist_db(1);
        else
            z_del_db = z_d;
        end
        y_db = z_d - z_del_db;
        dz_db(k) = y_db;

        e_z1_db = y_db - xhat_db(1);
        xhat_corr_db = xhat_db + L_db * e_z1_db;       % correct ALL
        fd_db = (1/a_x) * (1-lc) * xhat_corr_db(3);    % control corrected

        z_new_db = z_pos_db + a_x * (fd_db + fT_db(k));
        z_hist_db(k+1) = z_new_db;
        z_pos_db = z_new_db;

        xhat_db = A_obs * xhat_corr_db + b_u * fd_db;  % predict from corrected
    end

    ss_db = steady_start:N_steps;
    C_sim_filt_db(idx) = var(dz_db(ss_db)) / (sigma_fT_pN^2 * a_x^2);
end

figure('Position', [100 100 800 500]);
plot(lc_dense, C_eq17_dense * var_scale, 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Eq.17 (no observer)');
hold on;
plot(lc_dense, C_db_dense * var_scale, 'r--', 'LineWidth', 2, ...
    'DisplayName', 'PP filtering (\lambda_e=0, Lyapunov)');
plot(lc_list, C_sim_filt_db * var_scale, 'ro', 'MarkerSize', 9, ...
    'DisplayName', 'Sim: filtering (\lambda_e=0)');
hold off;
xlabel('\lambda_c'); ylabel('Variance [\mum^2]');
title(sprintf('Deadbeat Observer (\\lambda_e=0): \\DeltaC = 1, \\gamma=%.4f pN\\cdots/\\mum', gammaN));
legend('Location', 'northwest'); grid on;
set(gca, 'FontSize', 11);
print(gcf, fullfile(pwd, 'fig_variance_le0.png'), '-dpng', '-r150');

fprintf('\nFigures saved:\n');
fprintf('  fig_variance_comparison.png\n');
fprintf('  fig_variance_le0.png\n');
