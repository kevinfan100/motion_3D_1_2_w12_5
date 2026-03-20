% verify_Cdx.m — Verify C_dx(lc, a_var) formula against Monte Carlo simulation
%
% IIR filter definition (matching C_dx formula convention):
%   dzm_bar[k] = a_var * dzm_bar[k-1] + (1-a_var) * dzm[k]
%   a_var=0: no filter (= Eq.17), a_var->1: heavy filtering
%
% System:
%   Plant:      x[k+1] = x[k] + a_x*(fd[k] + fT[k])
%   Measurement: dzm[k] = x_d - x[k-2]              (2-step delay)
%   IIR:        dzm_bar[k] = a_var*dzm_bar[k-1] + (1-a_var)*dzm[k]
%   Controller: fd[k] = (1/a_x)*(1-lc)*dzm_bar[k] - (1-lc)*(fd[k-1]+fd[k-2])
%
% C_dx formula:
%   C_dx(lc, av) = 2*(1-av)*(1-lc) / [1-(1-av)*lc]
%                + (2/(2-av)) / [(1+lc)*(1-(1-av)*lc)]
%   When av=0: C_dx = 2 + 1/(1-lc^2) = C(lc)

clear; clc; close all;

%% ===== Physical Parameters =====
Ts      = 1/1600;
kb      = 1.3806503e-23;
T_temp  = 310.15;
gammaN  = 0.0425;
a_x     = Ts / gammaN;

sigma2_fT_SI = 4 * kb * T_temp * (gammaN * 1e-6) / Ts;
sigma_fT_pN  = sqrt(sigma2_fT_SI * 1e24);
base_var     = 4 * kb * T_temp * a_x * 1e18;

%% ===== Test Parameters =====
a_var_list    = [0, 0.005, 0.5];
lambda_c_list = [0.3, 0.6, 0.9];
n_av = length(a_var_list);
n_lc = length(lambda_c_list);

N_steps      = 80000;
steady_start = 40000;
z_d          = 25;

%% ===== C factor functions =====
C_func    = @(lc)     2 + 1 / (1 - lc^2);
C_dx_func = @(lc, av) 2*(1-av)*(1-lc) / (1-(1-av)*lc) ...
                     + (2/(2-av)) / ((1+lc)*(1-(1-av)*lc));

%% ===== Observer Parameters =====
% 3-state observer pole placement (same as verify_eq13_unified.m)
lambda_e = 0.3;
L1 = 1 - 3*lambda_e;
L2 = 1 - 3*lambda_e + 3*lambda_e^2;
L3 = (1 - lambda_e)^3;

%% ===== Preallocate =====
C_formula  = zeros(n_av, n_lc);
C_lyapunov = zeros(n_av, n_lc);
C_sim      = zeros(n_av, n_lc);
C_sim_dzr  = zeros(n_av, n_lc);
C_sim_obs  = zeros(n_av, n_lc);
C_sim_dzr_obs = zeros(n_av, n_lc);

%% ===== Main loop =====
for ai = 1:n_av
    av = a_var_list(ai);

    for idx = 1:n_lc
        lc = lambda_c_list(idx);

        % --- C_dx formula ---
        C_formula(ai, idx) = C_dx_func(lc, av);

        % --- 6-state Lyapunov (ground truth) ---
        p = 1 - av;
        A = [1,  0, -(1-lc)*p, -(1-lc)*av,  (1-lc),   (1-lc);
             1,  0,  0,         0,           0,         0;
             0,  1,  0,         0,           0,         0;
             0,  0,  p,         av,          0,         0;
             0,  0,  (1-lc)*p,  (1-lc)*av, -(1-lc),  -(1-lc);
             0,  0,  0,         0,           1,         0];
        B_vec = [-1; 0; 0; 0; 0; 0];
        P = dlyap(A, B_vec*B_vec');
        C_lyapunov(ai, idx) = P(1,1);

        % --- Monte Carlo simulation ---
        rng(42 + idx);
        fT = sigma_fT_pN * randn(N_steps, 1);

        z      = z_d;
        z_hist = z_d * ones(N_steps, 1);
        fd_hist = zeros(N_steps, 1);
        dzm_all = zeros(N_steps, 1);
        dzr_all = zeros(N_steps, 1);
        dzm_bar = 0;

        for k = 1:N_steps
            % --- Measurement: 2-step delay ---
            if k >= 3,     z_delayed = z_hist(k-2);
            elseif k == 2, z_delayed = z_hist(1);
            else,          z_delayed = z_d;
            end
            dzm = z_d - z_delayed;
            dzm_all(k) = dzm;

            % --- IIR filter ---
            % a_var=0: dzm_bar = dzm (no filter = Eq.17)
            % a_var>0: smoothing
            dzm_bar = av * dzm_bar + (1 - av) * dzm;
            dzr_all(k) = dzm - dzm_bar;   % high-pass residual

            % --- Controller (Eq.17 form with dzm_bar) ---
            if k >= 2, fd_km1 = fd_hist(k-1); else, fd_km1 = 0; end
            if k >= 3, fd_km2 = fd_hist(k-2); else, fd_km2 = 0; end
            fd_k = (1/a_x)*(1-lc)*dzm_bar - (1-lc)*(fd_km1 + fd_km2);
            fd_hist(k) = fd_k;

            % --- Plant ---
            z_new = z + a_x * (fd_k + fT(k));
            if k < N_steps, z_hist(k+1) = z_new; end
            z = z_new;
        end

        % --- Steady-state variance ---
        ss = steady_start:N_steps;
        var_dzm = var(dzm_all(ss));
        C_sim(ai, idx) = var_dzm / (sigma_fT_pN^2 * a_x^2);
        var_dzr = var(dzr_all(ss));
        C_sim_dzr(ai, idx) = var_dzr / (sigma_fT_pN^2 * a_x^2);

        % --- Observer+IIR Monte Carlo (Simulink-style) ---
        rng(42 + idx);
        fT_obs = sigma_fT_pN * randn(N_steps, 1);

        z_obs      = z_d;
        z_hist_obs = z_d * ones(N_steps, 1);
        dzm_all_obs = zeros(N_steps, 1);
        dzr_all_obs = zeros(N_steps, 1);
        dz1_hat = 0;  dz2_hat = 0;  dz3_hat = 0;
        dzm_bar_obs = 0;

        for k = 1:N_steps
            % Measurement: 2-step delay
            if k >= 3,     z_del_obs = z_hist_obs(k-2);
            elseif k == 2, z_del_obs = z_hist_obs(1);
            else,          z_del_obs = z_d;
            end
            dzm_obs = z_d - z_del_obs;
            dzm_all_obs(k) = dzm_obs;

            % IIR filter (same as Eq.17 path)
            dzm_bar_obs = av * dzm_bar_obs + (1 - av) * dzm_obs;
            dzr_all_obs(k) = dzm_obs - dzm_bar_obs;   % high-pass residual

            % Innovation uses dzm_bar (Simulink x-axis style)
            innov = dzm_bar_obs - dz1_hat;

            % Control law (observer-based, known a_x)
            fd_k_obs = (1/a_x) * (1 - lc) * dz3_hat;

            % 3-state observer update
            dz1_new = dz2_hat      + L1 * innov;
            dz2_new = dz3_hat      + L2 * innov;
            dz3_new = lc * dz3_hat + L3 * innov;
            dz1_hat = dz1_new;  dz2_hat = dz2_new;  dz3_hat = dz3_new;

            % Plant
            z_new_obs = z_obs + a_x * (fd_k_obs + fT_obs(k));
            if k < N_steps, z_hist_obs(k+1) = z_new_obs; end
            z_obs = z_new_obs;
        end

        var_dzm_obs = var(dzm_all_obs(ss));
        C_sim_obs(ai, idx) = var_dzm_obs / (sigma_fT_pN^2 * a_x^2);
        var_dzr_obs = var(dzr_all_obs(ss));
        C_sim_dzr_obs(ai, idx) = var_dzr_obs / (sigma_fT_pN^2 * a_x^2);
    end
end

%% ===== Results =====
fprintf('==================== C_dx VERIFICATION ====================\n');
fprintf('IIR: dzm_bar[k] = a_var * dzm_bar[k-1] + (1-a_var) * dzm[k]\n');
fprintf('base_var = %.4e um^2\n\n', base_var);

for ai = 1:n_av
    av = a_var_list(ai);
    if av == 0
        fprintf('--- a_var = 0  (no filter = Eq.17) ---\n');
    else
        fprintf('--- a_var = %.3f  (%.1f%% old + %.1f%% new) ---\n', ...
            av, av*100, (1-av)*100);
    end
    fprintf('  %-6s  %-10s  %-10s  %-14s  %-14s  %-14s  %-14s  %-14s  %-14s\n', ...
        'lc', 'C_dx', 'C_lyap', 'formula(um^2)', 'lyap(um^2)', ...
        'sim_dzm(um^2)', 'sim_dzr(um^2)', 'sim_obs(um^2)', 'dzr_obs(um^2)');
    fprintf('  %s\n', repmat('-', 1, 135));
    for idx = 1:n_lc
        lc = lambda_c_list(idx);
        cf = C_formula(ai, idx);
        cl = C_lyapunov(ai, idx);
        cs = C_sim(ai, idx);
        cr = C_sim_dzr(ai, idx);
        co = C_sim_obs(ai, idx);
        cor = C_sim_dzr_obs(ai, idx);
        vf = cf * base_var;
        vl = cl * base_var;
        vs = cs * base_var;
        vr = cr * base_var;
        vo = co * base_var;
        vor = cor * base_var;
        if av == 0
            dzr_str = 'N/A (dzr=0)   ';
            dzr_obs_str = 'N/A (dzr=0)   ';
        else
            dzr_str = sprintf('%-14.4e', vr);
            dzr_obs_str = sprintf('%-14.4e', vor);
        end
        fprintf('  %-6.1f  %-10.4f  %-10.4f  %-14.4e  %-14.4e  %-14.4e  %s  %-14.4e  %s\n', ...
            lc, cf, cl, vf, vl, vs, dzr_str, vo, dzr_obs_str);
    end
    % Error summary for dzr vs C_dx (skip av=0)
    if av > 0
        fprintf('  dzr vs C_dx errors:');
        for idx = 1:n_lc
            vf = C_formula(ai, idx) * base_var;
            vr = C_sim_dzr(ai, idx) * base_var;
            err_dzr = (vr - vf) / vf * 100;
            fprintf('  lc=%.1f: %+.2f%%', lambda_c_list(idx), err_dzr);
        end
        fprintf('\n');
    end
    fprintf('\n');
end

%% ===== Figure =====
fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
lc_dense = 0.1:0.005:0.95;

colors = [0 0.45 0.74;     % blue
          0.85 0.33 0.1;    % red
          0.47 0.67 0.19];  % green
markers = {'o', 's', '^'};
line_styles = {'-', '--', ':'};

fig1 = figure('Position', [50 50 950 600], 'Color', 'w');
hold on;

for ai = 1:n_av
    av = a_var_list(ai);

    % Formula curve
    C_f_dense = arrayfun(@(lc) C_dx_func(lc, av), lc_dense);
    sigma2_f = C_f_dense * base_var;
    if av == 0
        lbl = 'a_{var}=0  (no filter)';
    else
        lbl = sprintf('a_{var}=%.3f', av);
    end
    plot(lc_dense, sigma2_f, line_styles{ai}, 'Color', colors(ai,:), ...
        'LineWidth', 2.5, 'DisplayName', ['C_{dx}: ' lbl]);

    % Sim markers
    sigma2_sim_pts = C_sim(ai, :) * base_var;
    plot(lambda_c_list, sigma2_sim_pts, markers{ai}, ...
        'Color', colors(ai,:), 'MarkerSize', 11, ...
        'MarkerFaceColor', colors(ai,:), 'LineWidth', 1.5, ...
        'DisplayName', ['Sim: ' lbl]);
end
hold off;

xlabel('\lambda_c', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\sigma^2_{dzm}  (\mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
title('C_{dx}(\lambda_c, a_{var}) formula vs simulation', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;

saveas(fig1, fullfile(fig_dir, 'fig_Cdx_verify.png'));
fprintf('Figure saved to figures/fig_Cdx_verify.png\n');

%% ===== Figure 2: with Lyapunov ground truth =====
% Dense Lyapunov curves
C_lyap_dense = zeros(n_av, length(lc_dense));
for ai = 1:n_av
    av = a_var_list(ai);
    p_av = 1 - av;
    for j = 1:length(lc_dense)
        lc_j = lc_dense(j);
        A_j = [1,  0, -(1-lc_j)*p_av, -(1-lc_j)*av,  (1-lc_j),   (1-lc_j);
               1,  0,  0,              0,              0,           0;
               0,  1,  0,              0,              0,           0;
               0,  0,  p_av,           av,             0,           0;
               0,  0,  (1-lc_j)*p_av, (1-lc_j)*av,  -(1-lc_j),  -(1-lc_j);
               0,  0,  0,              0,              1,           0];
        B_j = [-1; 0; 0; 0; 0; 0];
        P_j = dlyap(A_j, B_j*B_j');
        C_lyap_dense(ai, j) = P_j(1,1);
    end
end

fig2 = figure('Position', [50 50 950 600], 'Color', 'w');
hold on;

for ai = 1:n_av
    av = a_var_list(ai);

    % Lyapunov (solid thick = ground truth)
    sigma2_lyap = C_lyap_dense(ai, :) * base_var;
    if av == 0
        lbl_l = 'Lyapunov: a_{var}=0  (= Eq.17)';
    else
        lbl_l = sprintf('Lyapunov: a_{var}=%.3f', av);
    end
    plot(lc_dense, sigma2_lyap, '-', 'Color', colors(ai,:), ...
        'LineWidth', 2.5, 'DisplayName', lbl_l);

    % Formula (dashed)
    C_f_dense = arrayfun(@(lc) C_dx_func(lc, av), lc_dense);
    sigma2_f = C_f_dense * base_var;
    plot(lc_dense, sigma2_f, '--', 'Color', colors(ai,:)*0.7, ...
        'LineWidth', 1.8, 'DisplayName', sprintf('Formula: a_{var}=%.3f', av));

    % Sim markers
    sigma2_sim_pts = C_sim(ai, :) * base_var;
    plot(lambda_c_list, sigma2_sim_pts, markers{ai}, ...
        'Color', colors(ai,:), 'MarkerSize', 11, ...
        'MarkerFaceColor', colors(ai,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Sim: a_{var}=%.3f', av));
end
hold off;

xlabel('\lambda_c', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\sigma^2_{dzm}  (\mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
title('C_{dx} : solid = Lyapunov (exact),  dashed = formula,  markers = sim', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;

saveas(fig2, fullfile(fig_dir, 'fig_Cdx_lyapunov.png'));
fprintf('Figure saved to figures/fig_Cdx_lyapunov.png\n');

%% ===== Figure 3: Observer+IIR simulation =====
fig3 = figure('Position', [50 50 950 600], 'Color', 'w');
hold on;

for ai = 1:n_av
    av = a_var_list(ai);

    % Lyapunov (solid = ground truth for Eq.17+IIR system)
    sigma2_lyap = C_lyap_dense(ai, :) * base_var;
    if av == 0
        lbl_l = 'Lyapunov(Eq.17): a_{var}=0';
    else
        lbl_l = sprintf('Lyapunov(Eq.17): a_{var}=%.3f', av);
    end
    plot(lc_dense, sigma2_lyap, '-', 'Color', colors(ai,:), ...
        'LineWidth', 2.5, 'DisplayName', lbl_l);

    % Eq.17 sim markers (filled)
    sigma2_sim_pts = C_sim(ai, :) * base_var;
    plot(lambda_c_list, sigma2_sim_pts, markers{ai}, ...
        'Color', colors(ai,:), 'MarkerSize', 10, ...
        'MarkerFaceColor', colors(ai,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Sim(Eq.17): a_{var}=%.3f', av));

    % Observer sim markers (open)
    sigma2_obs_pts = C_sim_obs(ai, :) * base_var;
    plot(lambda_c_list, sigma2_obs_pts, markers{ai}, ...
        'Color', colors(ai,:)*0.6, 'MarkerSize', 12, ...
        'MarkerFaceColor', 'none', 'LineWidth', 2, ...
        'DisplayName', sprintf('Sim(Obs+IIR): a_{var}=%.3f', av));
end
hold off;

xlabel('\lambda_c', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\sigma^2_{dzm}  (\mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
title('Observer+IIR vs Eq.17+IIR: filled=Eq.17, open=Observer', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;

saveas(fig3, fullfile(fig_dir, 'fig_Cdx_observer.png'));
fprintf('Figure saved to figures/fig_Cdx_observer.png\n');

%% ===== Figure 4: C_dx vs Var(dzr) — IIR high-pass residual =====
% Only plot av > 0 (av=0 gives dzr=0 trivially)
av_plot_idx = find(a_var_list > 0);
colors_dzr = [0.85 0.33 0.1;    % red
              0.47 0.67 0.19;    % green
              0 0.45 0.74];      % blue (spare)

fig4 = figure('Position', [50 50 950 600], 'Color', 'w');
hold on;

for ii = 1:length(av_plot_idx)
    ai = av_plot_idx(ii);
    av = a_var_list(ai);
    col = colors_dzr(ii, :);

    % C_dx formula curve
    C_f_dense = arrayfun(@(lc) C_dx_func(lc, av), lc_dense);
    sigma2_f = C_f_dense * base_var;
    plot(lc_dense, sigma2_f, '-', 'Color', col, ...
        'LineWidth', 2.5, 'DisplayName', sprintf('C_{dx} formula: a_{var}=%.3f', av));

    % Var(dzr) MC markers (Eq.17+IIR)
    sigma2_dzr_pts = C_sim_dzr(ai, :) * base_var;
    plot(lambda_c_list, sigma2_dzr_pts, 'o', ...
        'Color', col, 'MarkerSize', 11, ...
        'MarkerFaceColor', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Var(dzr) Eq.17: a_{var}=%.3f', av));

    % Var(dzr) MC markers (Observer+IIR)
    sigma2_dzr_obs_pts = C_sim_dzr_obs(ai, :) * base_var;
    plot(lambda_c_list, sigma2_dzr_obs_pts, 's', ...
        'Color', col*0.6, 'MarkerSize', 12, ...
        'MarkerFaceColor', 'none', 'LineWidth', 2, ...
        'DisplayName', sprintf('Var(dzr) Obs: a_{var}=%.3f', av));

    % Lyapunov reference (dashed)
    sigma2_lyap = C_lyap_dense(ai, :) * base_var;
    plot(lc_dense, sigma2_lyap, '--', 'Color', col*0.5, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Lyapunov(dzm): a_{var}=%.3f', av));
end
hold off;

xlabel('\lambda_c', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Tracking error variance  (\mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
title('C_{dx} formula vs Var(dzr) — IIR high-pass residual (a_{var}>0 only)', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;

saveas(fig4, fullfile(fig_dir, 'fig_Cdx_dzr.png'));
fprintf('Figure saved to figures/fig_Cdx_dzr.png\n');
