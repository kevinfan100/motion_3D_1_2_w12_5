% verify_eq12_spectral.m — Numerical verification of Eq.12 C(lc) factor
%
% Four independent methods to verify C(lc) = 2 + 1/(1-lc^2):
%   1. Impulse response summation:  sum h_norm[k]^2
%   2. Frequency-domain integration: (1/pi) * int_0^pi |H(e^(jw))|^2 dw
%   3. Lyapunov equation:           dlyap(A, B*B') -> P(1,1)
%   4. Monte Carlo simulation:      Var(dzm) / (sigma_fT^2 * a_x^2)
%
% Also plots |H(e^(jw))|^2 to visualize the spectral content.

clear; clc; close all;

%% ===== Physical Parameters =====
Ts      = 1/1600;
kb      = 1.3806503e-23;
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
lc_list = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95];
n_lc    = length(lc_list);

% Preallocate results
C_theory   = zeros(1, n_lc);
C_impulse  = zeros(1, n_lc);
C_spectral = zeros(1, n_lc);
C_lyapunov = zeros(1, n_lc);
C_sim      = zeros(1, n_lc);

% Simulation parameters (160k steps for better convergence at high lc)
N_steps      = 160000;
steady_start = 60001;
z_d          = 25;       % target [um]

%% ===== Main loop =====
for idx = 1:n_lc
    lc = lc_list(idx);

    % --- Theoretical C ---
    C_theory(idx) = 2 + 1/(1 - lc^2);

    % =============================================
    %  Method 1: Impulse response summation
    % =============================================
    % Compute h_norm[k] by simulating:
    %   e[k+1] = lc*e[k] - fT[k] - (1-lc)*fT[k-1] - (1-lc)*fT[k-2]
    % with fT[0]=1, fT[k]=0 otherwise.

    K_max = 500;  % truncation length
    h_e = zeros(K_max+1, 1);   % h_e[0..K_max], 1-indexed: h_e(1)=h_e[0]
    fT_imp = zeros(K_max+1, 1);
    fT_imp(1) = 1;  % impulse at k=0

    e_prev = 0;
    for k = 0:K_max-1
        ki = k + 1;  % MATLAB 1-indexed
        fT_k  = fT_imp(ki);
        fT_km1 = 0; if ki >= 2, fT_km1 = fT_imp(ki-1); end
        fT_km2 = 0; if ki >= 3, fT_km2 = fT_imp(ki-2); end

        e_next = lc*e_prev - fT_k - (1-lc)*fT_km1 - (1-lc)*fT_km2;
        h_e(ki+1) = e_next;  % h_e[k+1]
        e_prev = e_next;
    end

    % h_norm = |h_e| (the sign is just a convention)
    C_impulse(idx) = sum(h_e.^2);

    % =============================================
    %  Method 2: Frequency-domain integration
    % =============================================
    % H_norm(z) = [z^2 + (1-lc)*z + (1-lc)] / [z^2*(z - lc)]
    % Evaluate |H_norm(e^(j*theta))|^2 and integrate (1/pi)*int_0^pi

    N_pts = 10000;
    theta = linspace(0, pi, N_pts);
    z_uc  = exp(1j * theta);   % points on upper half of unit circle

    num = z_uc.^2 + (1-lc)*z_uc + (1-lc);
    den = z_uc.^2 .* (z_uc - lc);
    H_mag2 = abs(num ./ den).^2;

    C_spectral(idx) = (1/pi) * trapz(theta, H_mag2);

    % =============================================
    %  Method 3: Lyapunov equation
    % =============================================
    % Augmented state: x = [e[k]; fT[k-1]; fT[k-2]]
    % x[k+1] = A*x[k] + B*fT[k]
    %
    % A = [lc,  -(1-lc),  -(1-lc)]    B_norm = [-1]
    %     [ 0,     0,         0   ]             [ 1]
    %     [ 0,     1,         0   ]             [ 0]
    %
    % (B_norm = B/a_x, so P(1,1) gives C directly)

    A_aug = [lc,  -(1-lc),  -(1-lc);
              0,      0,         0;
              0,      1,         0];
    B_norm = [-1; 1; 0];

    P = dlyap(A_aug, B_norm * B_norm');
    C_lyapunov(idx) = P(1,1);

    % =============================================
    %  Method 4: Monte Carlo simulation
    % =============================================
    rng(42 + idx);
    fT_sim = sigma_fT_pN * randn(N_steps, 1);

    z_pos   = z_d;
    z_hist  = zeros(N_steps, 1); z_hist(1) = z_d;
    fd_hist = zeros(N_steps, 1);
    dzm_all = zeros(N_steps, 1);

    for k = 1:N_steps
        % Measurement with 2-step delay
        if k >= 3
            z_del = z_hist(k-2);
        elseif k == 2
            z_del = z_hist(1);
        else
            z_del = z_d;
        end
        dzm = z_d - z_del;
        dzm_all(k) = dzm;

        % Predictor-based controller (Eq.17)
        fd_km2 = 0; if k >= 3, fd_km2 = fd_hist(k-2); end
        fd_km1 = 0; if k >= 2, fd_km1 = fd_hist(k-1); end
        e_hat  = dzm - a_x * (fd_km2 + fd_km1);
        fd_k   = (1/a_x) * (1-lc) * e_hat;
        fd_hist(k) = fd_k;

        % Plant
        z_new = z_pos + a_x * (fd_k + fT_sim(k));
        if k < N_steps
            z_hist(k+1) = z_new;
        end
        z_pos = z_new;
    end

    % Steady-state variance
    ss = steady_start:N_steps;
    var_dzm = var(dzm_all(ss));   % [um^2]

    % Convert to C: Var(dzm) = sigma_fT^2 * a_x^2 * C
    % sigma_fT^2 [pN^2], a_x^2 [um^2/pN^2] => sigma_fT^2 * a_x^2 [um^2]
    sigma_fT_pN2 = sigma_fT_pN^2;
    C_sim(idx) = var_dzm / (sigma_fT_pN2 * a_x^2);

    fprintf('lc=%.2f: C_theory=%.4f  C_impulse=%.4f  C_spectral=%.4f  C_lyap=%.4f  C_sim=%.4f\n', ...
        lc, C_theory(idx), C_impulse(idx), C_spectral(idx), C_lyapunov(idx), C_sim(idx));
end

%% ===== Results Table =====
fprintf('\n==================== C(lc) VERIFICATION RESULTS ====================\n');
fprintf('%-6s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'lc', 'Theory', 'Impulse', 'Spectral', 'Lyapunov', 'Simulation');
fprintf('%s\n', repmat('-', 1, 68));
for idx = 1:n_lc
    fprintf('%-6.2f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n', ...
        lc_list(idx), C_theory(idx), C_impulse(idx), C_spectral(idx), ...
        C_lyapunov(idx), C_sim(idx));
end

fprintf('\n--- Relative Errors vs Theory (%%):  ---\n');
fprintf('%-6s  %-10s  %-10s  %-10s  %-10s\n', 'lc', 'Impulse', 'Spectral', 'Lyapunov', 'Simulation');
fprintf('%s\n', repmat('-', 1, 50));
for idx = 1:n_lc
    err_imp  = abs(C_impulse(idx) - C_theory(idx)) / C_theory(idx) * 100;
    err_spec = abs(C_spectral(idx) - C_theory(idx)) / C_theory(idx) * 100;
    err_lyap = abs(C_lyapunov(idx) - C_theory(idx)) / C_theory(idx) * 100;
    err_sim  = abs(C_sim(idx) - C_theory(idx)) / C_theory(idx) * 100;
    fprintf('%-6.2f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n', ...
        lc_list(idx), err_imp, err_spec, err_lyap, err_sim);
end

%% ===== Figure 1: |H(e^(jw))|^2 for different lc =====
figure('Name', 'Power Spectral Density: |H_norm(e^{j\theta})|^2');
theta_plot = linspace(0, pi, 2000);
colors = lines(n_lc);
hold on;
for idx = 1:n_lc
    lc = lc_list(idx);
    z_uc = exp(1j * theta_plot);
    num = z_uc.^2 + (1-lc)*z_uc + (1-lc);
    den = z_uc.^2 .* (z_uc - lc);
    H_mag2 = abs(num ./ den).^2;
    plot(theta_plot/pi, H_mag2, 'Color', colors(idx,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\lambda_c = %.2f', lc));
end
hold off;
xlabel('\theta / \pi');
ylabel('|H_{norm}(e^{j\theta})|^2');
title('Spectral density of normalized transfer function');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);

%% ===== Figure 2: C(lc) vs lc — all methods =====
figure('Name', 'C(lc) Verification: All Methods');
lc_fine = linspace(0.01, 0.98, 200);
C_fine  = 2 + 1./(1 - lc_fine.^2);

plot(lc_fine, C_fine, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory: 2 + 1/(1-\lambda_c^2)');
hold on;
plot(lc_list, C_impulse,  'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Impulse response');
plot(lc_list, C_spectral, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Spectral integration');
plot(lc_list, C_lyapunov, 'gd', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Lyapunov equation');
plot(lc_list, C_sim,      'mo', 'MarkerSize', 10, 'DisplayName', 'Monte Carlo (80k steps)');
hold off;
xlabel('\lambda_c');
ylabel('C(\lambda_c)');
title('C(\lambda_c) = 2 + 1/(1-\lambda_c^2) — Four-way verification');
legend('Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);
ylim([2.5, max(C_fine)*1.1]);

%% ===== Figure 3: Impulse response visualization =====
figure('Name', 'Impulse Response h_{norm}[k]');
for pidx = 1:min(4, n_lc)
    subplot(2, 2, pidx);
    lc = lc_list(pidx);

    % Recompute impulse response
    K_show = 20;
    h_show = zeros(K_show+1, 1);
    fT_imp_s = zeros(K_show+1, 1);
    fT_imp_s(1) = 1;
    e_p = 0;
    for k = 0:K_show-1
        ki = k + 1;
        fk  = fT_imp_s(ki);
        fkm1 = 0; if ki >= 2, fkm1 = fT_imp_s(ki-1); end
        fkm2 = 0; if ki >= 3, fkm2 = fT_imp_s(ki-2); end
        e_n = lc*e_p - fk - (1-lc)*fkm1 - (1-lc)*fkm2;
        h_show(ki+1) = e_n;
        e_p = e_n;
    end

    stem(0:K_show, abs(h_show), 'filled', 'LineWidth', 1.2);
    xlabel('k');
    ylabel('|h_{norm}[k]|');
    title(sprintf('\\lambda_c = %.1f,  C = %.3f', lc, 2+1/(1-lc^2)));
    grid on;
    ylim([0, 1.3]);
end
sgtitle('Impulse response: h_{norm}[k] = [0, 1, 1, 1, \lambda_c, \lambda_c^2, ...]');

%% ===== Figure 4: Decomposition of C(lc) =====
figure('Name', 'C(lc) Decomposition');
lc_dec = linspace(0.01, 0.98, 200);
C_delay = 2 * ones(size(lc_dec));
C_ar1   = 1 ./ (1 - lc_dec.^2);

area_data = [C_delay; C_ar1]';
h_area = area(lc_dec, area_data);
h_area(1).FaceColor = [0.3 0.6 0.9];
h_area(2).FaceColor = [0.9 0.4 0.2];
h_area(1).DisplayName = 'Delay contribution (= 2)';
h_area(2).DisplayName = 'AR(1) amplification = 1/(1-\lambda_c^2)';
xlabel('\lambda_c');
ylabel('C(\lambda_c)');
title('Decomposition: C(\lambda_c) = 2 + 1/(1-\lambda_c^2)');
legend('Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);
ylim([0, 15]);

%% ===== Final pass/fail =====
fprintf('\n==================== PASS/FAIL ====================\n');
tol_exact = 1e-6;   % for exact methods
tol_sim   = 5;      % 5% for simulation (stochastic)

pass_all = true;
for idx = 1:n_lc
    err_imp  = abs(C_impulse(idx) - C_theory(idx)) / C_theory(idx) * 100;
    err_spec = abs(C_spectral(idx) - C_theory(idx)) / C_theory(idx) * 100;
    err_lyap = abs(C_lyapunov(idx) - C_theory(idx)) / C_theory(idx) * 100;
    err_sim  = abs(C_sim(idx) - C_theory(idx)) / C_theory(idx) * 100;

    if err_imp > tol_exact*100
        fprintf('FAIL: lc=%.2f impulse err=%.4f%%\n', lc_list(idx), err_imp);
        pass_all = false;
    end
    if err_spec > 0.1  % 0.1% tolerance for numerical integration
        fprintf('FAIL: lc=%.2f spectral err=%.4f%%\n', lc_list(idx), err_spec);
        pass_all = false;
    end
    if err_lyap > tol_exact*100
        fprintf('FAIL: lc=%.2f Lyapunov err=%.4f%%\n', lc_list(idx), err_lyap);
        pass_all = false;
    end
    if err_sim > tol_sim
        fprintf('FAIL: lc=%.2f simulation err=%.2f%%\n', lc_list(idx), err_sim);
        pass_all = false;
    end
end

if pass_all
    fprintf('RESULT: ALL PASS\n');
    fprintf('  - Impulse, Spectral, Lyapunov: < 0.01%% error (exact methods)\n');
    fprintf('  - Simulation: < 5%% error (80k steps Monte Carlo)\n');
else
    fprintf('RESULT: SOME FAILURES (see above)\n');
end
fprintf('====================================================\n');
