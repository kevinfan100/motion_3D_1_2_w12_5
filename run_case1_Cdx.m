% run_case1_Cdx.m — Simulink cross-verification of C_dx formula
%
% Sweeps Avar = [0.005, 0.5] x lc = [0.3, 0.6, 0.9] using the actual
% Simulink model. Compares Var(dz_k2) with C_dx formula, 6-state Lyapunov,
% and verify_Cdx.m Monte Carlo results.
%
% The Simulink z-axis EKF innovation uses raw dz[k-2] (no IIR).
% Avar only affects azm_k (IIR mobility estimate) indirectly.
%
% Usage: open Simulink desktop, then run this script in MATLAB.

clear; clc; close all;

%% ===== Physical Constants =====
Ts      = 1/1600;
kb      = 1.3806503e-23;
T_temp  = 310.15;
gammaN  = 0.0425;
a_x     = Ts / gammaN;

sigma2_fT_SI = 4 * kb * T_temp * (gammaN * 1e-6) / Ts;
sigma_fT_pN  = sqrt(sigma2_fT_SI * 1e24);
base_var     = 4 * kb * T_temp * a_x * 1e18;   % [um^2]

%% ===== Sweep Parameters =====
Avar_list     = [0.005, 0.5];
lambda_c_list = [0.3, 0.6, 0.9];
Am_scaling_val = 1;

n_av = length(Avar_list);
n_lc = length(lambda_c_list);

%% ===== C_dx formula =====
C_dx_func = @(lc, av) 2*(1-av)*(1-lc) / (1-(1-av)*lc) ...
                     + (2/(2-av)) / ((1+lc)*(1-(1-av)*lc));

%% ===== Preallocate =====
C_formula    = zeros(n_av, n_lc);
C_lyapunov   = zeros(n_av, n_lc);
C_sim_eq17   = zeros(n_av, n_lc);
C_sim_slx    = zeros(n_av, n_lc);

%% ===== C_dx formula + 6-state Lyapunov =====
for ai = 1:n_av
    av = Avar_list(ai);
    for idx = 1:n_lc
        lc = lambda_c_list(idx);

        % C_dx formula
        C_formula(ai, idx) = C_dx_func(lc, av);

        % 6-state Lyapunov (Eq.17+IIR ground truth)
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
    end
end

%% ===== Monte Carlo (Eq.17+IIR, same as verify_Cdx.m) =====
N_steps      = 80000;
steady_start = 40000;
z_d          = 25;

for ai = 1:n_av
    av = Avar_list(ai);
    for idx = 1:n_lc
        lc = lambda_c_list(idx);

        rng(42 + idx);
        fT = sigma_fT_pN * randn(N_steps, 1);

        z      = z_d;
        z_hist = z_d * ones(N_steps, 1);
        fd_hist = zeros(N_steps, 1);
        dzm_all = zeros(N_steps, 1);
        dzm_bar = 0;

        for k = 1:N_steps
            if k >= 3,     z_delayed = z_hist(k-2);
            elseif k == 2, z_delayed = z_hist(1);
            else,          z_delayed = z_d;
            end
            dzm = z_d - z_delayed;
            dzm_all(k) = dzm;

            dzm_bar = av * dzm_bar + (1 - av) * dzm;

            if k >= 2, fd_km1 = fd_hist(k-1); else, fd_km1 = 0; end
            if k >= 3, fd_km2 = fd_hist(k-2); else, fd_km2 = 0; end
            fd_k = (1/a_x)*(1-lc)*dzm_bar - (1-lc)*(fd_km1 + fd_km2);
            fd_hist(k) = fd_k;

            z_new = z + a_x * (fd_k + fT(k));
            if k < N_steps, z_hist(k+1) = z_new; end
            z = z_new;
        end

        ss = steady_start:N_steps;
        var_dzm = var(dzm_all(ss));
        C_sim_eq17(ai, idx) = var_dzm / base_var;
    end
end
fprintf('Monte Carlo (Eq.17+IIR) done.\n');

%% ===== Simulink Setup =====
model = 'motion_3D_1_2_w12_5';
load_system(model);
open_system(model);
set_param(model, 'StopTime', '50');

% Get Stateflow chart handle
rt = sfroot;
chart = rt.find('-isa', 'Stateflow.EMChart', 'Path', [model '/Parameters']);

% Define ParametersBus (required by model, same as run_case1.m)
clear elms;
elms(1)  = Simulink.BusElement; elms(1).Name='Ts';            elms(1).DataType='double'; elms(1).Dimensions=1;
elms(2)  = Simulink.BusElement; elms(2).Name='lamdaC';        elms(2).DataType='double'; elms(2).Dimensions=1;
elms(3)  = Simulink.BusElement; elms(3).Name='theta';         elms(3).DataType='double'; elms(3).Dimensions=1;
elms(4)  = Simulink.BusElement; elms(4).Name='phi';           elms(4).DataType='double'; elms(4).Dimensions=1;
elms(5)  = Simulink.BusElement; elms(5).Name='pz';            elms(5).DataType='double'; elms(5).Dimensions=1;
elms(6)  = Simulink.BusElement; elms(6).Name='R';             elms(6).DataType='double'; elms(6).Dimensions=1;
elms(7)  = Simulink.BusElement; elms(7).Name='kb';            elms(7).DataType='double'; elms(7).Dimensions=1;
elms(8)  = Simulink.BusElement; elms(8).Name='T';             elms(8).DataType='double'; elms(8).Dimensions=1;
elms(9)  = Simulink.BusElement; elms(9).Name='ax_normal';     elms(9).DataType='double'; elms(9).Dimensions=1;
elms(10) = Simulink.BusElement; elms(10).Name='az_normal';    elms(10).DataType='double'; elms(10).Dimensions=1;
elms(11) = Simulink.BusElement; elms(11).Name='gammaN';       elms(11).DataType='double'; elms(11).Dimensions=1;
elms(12) = Simulink.BusElement; elms(12).Name='Avar';         elms(12).DataType='double'; elms(12).Dimensions=1;
elms(13) = Simulink.BusElement; elms(13).Name='Avar2';        elms(13).DataType='double'; elms(13).Dimensions=1;
elms(14) = Simulink.BusElement; elms(14).Name='Avar22';       elms(14).DataType='double'; elms(14).Dimensions=1;
elms(15) = Simulink.BusElement; elms(15).Name='Avar3';        elms(15).DataType='double'; elms(15).Dimensions=1;
elms(16) = Simulink.BusElement; elms(16).Name='Am_scaling';   elms(16).DataType='double'; elms(16).Dimensions=1;
elms(17) = Simulink.BusElement; elms(17).Name='beta';         elms(17).DataType='double'; elms(17).Dimensions=1;
elms(18) = Simulink.BusElement; elms(18).Name='lamdaF';       elms(18).DataType='double'; elms(18).Dimensions=1;
elms(19) = Simulink.BusElement; elms(19).Name='Pfz_11';       elms(19).DataType='double'; elms(19).Dimensions=1;
elms(20) = Simulink.BusElement; elms(20).Name='Pfz_22';       elms(20).DataType='double'; elms(20).Dimensions=1;
elms(21) = Simulink.BusElement; elms(21).Name='Pfz_33';       elms(21).DataType='double'; elms(21).Dimensions=1;
elms(22) = Simulink.BusElement; elms(22).Name='Pfz_44';       elms(22).DataType='double'; elms(22).Dimensions=1;
elms(23) = Simulink.BusElement; elms(23).Name='Pfz_55';       elms(23).DataType='double'; elms(23).Dimensions=1;
elms(24) = Simulink.BusElement; elms(24).Name='Pfz_66';       elms(24).DataType='double'; elms(24).Dimensions=1;
elms(25) = Simulink.BusElement; elms(25).Name='Pfz_77';       elms(25).DataType='double'; elms(25).Dimensions=1;
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

%% ===== Simulink Sweep =====
for ai = 1:n_av
    av = Avar_list(ai);
    for idx = 1:n_lc
        lc = lambda_c_list(idx);

        % Set Stateflow parameters
        set_sf_param(chart, 'lamdaC', lc);
        set_sf_param(chart, 'Avar', av);
        set_sf_param(chart, 'Am_scaling', Am_scaling_val);

        fprintf('Running Simulink: Avar=%.3f, lc=%.1f ...\n', av, lc);
        try
            out = sim(model);

            % Extract steady-state Var(dz_k2)
            N = length(dz_k2);
            ss = round(N/2):N;
            var_dz = var(dz_k2(ss));
            C_sim_slx(ai, idx) = var_dz / base_var;
            fprintf('  Var(dz_k2) = %.4e um^2,  C_slx = %.4f\n', var_dz, C_sim_slx(ai, idx));
        catch ME
            C_sim_slx(ai, idx) = NaN;
            fprintf('  FAILED: %s\n', ME.message);
        end
    end
end
fprintf('Simulink sweep done.\n\n');

%% ===== Results Table =====
fprintf('==================== C_dx SIMULINK CROSS-VERIFICATION ====================\n');
fprintf('base_var = %.4e um^2\n', base_var);
fprintf('Am_scaling = %d (Simulink)\n\n', Am_scaling_val);

fprintf('  %-7s %-5s  %-10s %-10s %-14s %-14s %-10s\n', ...
    'Avar', 'lc', 'C_formula', 'C_lyap', 'C_sim_eq17', 'C_sim_Slx', 'slx_err%');
fprintf('  %s\n', repmat('-', 1, 82));

for ai = 1:n_av
    av = Avar_list(ai);
    for idx = 1:n_lc
        lc = lambda_c_list(idx);
        cf = C_formula(ai, idx);
        cl = C_lyapunov(ai, idx);
        ce = C_sim_eq17(ai, idx);
        cs = C_sim_slx(ai, idx);
        % Error: Simulink vs Lyapunov
        if isnan(cs)
            fprintf('  %-7.3f %-5.1f  %-10.4f %-10.4f %-14.4f %-14s %s\n', ...
                av, lc, cf, cl, ce, 'FAILED', 'N/A');
        else
            err_sl = (cs - cl) / cl * 100;
            fprintf('  %-7.3f %-5.1f  %-10.4f %-10.4f %-14.4f %-14.4f %+.1f%%\n', ...
                av, lc, cf, cl, ce, cs, err_sl);
        end
    end
end
fprintf('\n');

% Also print variance table
fprintf('  %-7s %-5s  %-14s %-14s %-14s %-14s\n', ...
    'Avar', 'lc', 'V_formula', 'V_lyap', 'V_sim_eq17', 'V_sim_Slx');
fprintf('  %s\n', repmat('-', 1, 76));

for ai = 1:n_av
    av = Avar_list(ai);
    for idx = 1:n_lc
        lc = lambda_c_list(idx);
        vf = C_formula(ai, idx) * base_var;
        vl = C_lyapunov(ai, idx) * base_var;
        ve = C_sim_eq17(ai, idx) * base_var;
        vs = C_sim_slx(ai, idx) * base_var;
        fprintf('  %-7.3f %-5.1f  %-14.4e %-14.4e %-14.4e %-14.4e\n', ...
            av, lc, vf, vl, ve, vs);
    end
end

%% ===== Figure: fig_Cdx_simulink.png =====
% Only C_dx formula curves vs Simulink markers
lc_dense = 0.1:0.005:0.95;

colors = [0 0.45 0.74;      % blue  (Avar=0.005)
          0.85 0.33 0.1];    % red   (Avar=0.5)

fig1 = figure('Position', [50 50 1000 600], 'Color', 'w');
hold on;

for ai = 1:n_av
    av = Avar_list(ai);
    lbl_av = sprintf('A_{var}=%.3f', av);

    % C_dx formula (solid curve)
    C_f_dense = arrayfun(@(lc) C_dx_func(lc, av), lc_dense);
    sigma2_f = C_f_dense * base_var;
    plot(lc_dense, sigma2_f, '-', 'Color', colors(ai,:), ...
        'LineWidth', 2.5, 'DisplayName', ['C_{dx} formula: ' lbl_av]);

    % Simulink (filled markers) — skip NaN
    sigma2_slx = C_sim_slx(ai, :) * base_var;
    valid = ~isnan(sigma2_slx);
    if any(valid)
        plot(lambda_c_list(valid), sigma2_slx(valid), 'p', ...
            'Color', colors(ai,:), 'MarkerSize', 14, ...
            'MarkerFaceColor', colors(ai,:), 'LineWidth', 2, ...
            'DisplayName', ['Simulink Var(dz_{k2}): ' lbl_av]);
    end
end
hold off;

xlabel('\lambda_c', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Tracking error variance  (\mum^2)', 'FontSize', 14, 'FontWeight', 'bold');
title('C_{dx} formula vs Simulink Var(dz_{k2})', ...
    'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
grid on;

saveas(fig1, fullfile('figures', 'fig_Cdx_simulink.png'));
fprintf('\nFigure saved: figures/fig_Cdx_simulink.png\n');

%% ===== Local function =====
function set_sf_param(chart, param_name, value)
    script = chart.Script;
    lines = strsplit(script, newline);
    target = ['Parameters.' param_name];
    for i = 1:numel(lines)
        if contains(lines{i}, target) && contains(lines{i}, '=')
            lines{i} = sprintf('    %s = %g;', target, value);
            break;
        end
    end
    chart.Script = strjoin(lines, newline);
end
