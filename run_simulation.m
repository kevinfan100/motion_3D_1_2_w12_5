% RUN_SIMULATION  Main entry point for running the motion control simulation.
%
% Replaces path_control.m with modular architecture.
% Usage:
%   run_simulation            % default parameters
%   Edit user_config.m to change parameters, then re-run.

clear; clc; close all;

%% ===== Add paths =====
addpath('config');
addpath('controller');
addpath('wall_effect');
addpath('thermal_force');
addpath('trajectory');

%% ===== Build parameters =====
% Optional: pass overrides as struct, e.g.:
%   params = calc_simulation_params(struct('lamdaC', 0.5));
params = calc_simulation_params();

Ts = params.common.Ts;

%% ===== Load and configure model =====
model = 'motion_3D_1_2_w12_5';
load_system(model);
open_system(model);
set_param(model, 'StopTime', num2str(user_config().StopTime));

%% ===== Run simulation =====
fprintf('Running simulation...\n');
out = sim(model);
fprintf('Simulation complete.\n');

%% ===== Extract data =====
x_d_data  = x_d;
y_d_data  = y_d;
z_d_data  = z_d;
p_x_data  = p(:,1);
p_y_data  = p(:,2);
p_z_data  = p(:,3);
fd_x_data = squeeze(fd(1,1,:));
fd_y_data = squeeze(fd(2,1,:));
fd_z_data = squeeze(fd(3,1,:));

N = size(x_d_data, 1);
t = (0:N-1)' * Ts;

%% ===== Plot: X axis =====
figure;
subplot(311)
plot(t, x_d_data, 'o', t, p_x_data, 'r'); grid on;
title('Path in x'); xlabel('Time (s)'); ylabel('p_x (um)');
subplot(312)
plot(t, x_d_data - p_x_data, 'g'); grid on;
title('Motion error in x'); xlabel('Time (s)'); ylabel('error_x (um)');
subplot(313)
plot(t, fd_x_data, 'r'); grid on; axis([0 max(t) -10 10]);
title('Control effort in x'); xlabel('Time (s)'); ylabel('fd_x (pN)');

%% ===== Plot: Y axis =====
figure;
subplot(311)
plot(t, y_d_data, 'o', t, p_y_data, 'r'); grid on;
title('Path in y'); xlabel('Time (s)'); ylabel('p_y (um)');
subplot(312)
plot(t, y_d_data - p_y_data, 'g'); grid on;
title('Motion error in y'); xlabel('Time (s)'); ylabel('error_y (um)');
subplot(313)
plot(t, fd_y_data, 'r'); grid on; axis([0 max(t) -10 10]);
title('Control effort in y'); xlabel('Time (s)'); ylabel('fd_y (pN)');

%% ===== Plot: Z axis =====
figure;
subplot(311)
plot(t, z_d_data, 'o', t, p_z_data, 'r'); grid on;
title('Path in z'); xlabel('Time (s)'); ylabel('p_z (um)');
subplot(312)
plot(t, dz_k2, 'r'); grid on;
title('Motion error in z'); xlabel('Time (s)'); ylabel('error_z (um)');
subplot(313)
plot(t, fd_z_data, 'r'); grid on; axis([0 max(t) -10 10]);
title('Control effort in z'); xlabel('Time (s)'); ylabel('fd_z (pN)');

%% ===== Plot: Z motion gain =====
figure;
plot(t, azm_k, 'o', t, az_hat_k, 'r', t, mgain_z, 'g'); grid on;
title('Motion gain in z'); xlabel('Time (s)'); ylabel('az (um/(pN s))');

figure;
plot(t, az_hat_k, 'r', t, mgain_z, 'g'); grid on;
title('Motion gain in z (EKF vs plant)'); xlabel('Time (s)'); ylabel('az (um/(pN s))');

%% ===== Plot: Estimator feedback Lz =====
figure;
subplot(231)
plot(t, squeeze(Lz(1,1,:)), 'o', t, squeeze(Lz(2,1,:)), 'r', t, squeeze(Lz(3,1,:)), 'g'); grid on;
title('Lz col1 rows 1-3'); xlabel('Time (s)'); ylabel('Lz1');
subplot(234)
plot(t, squeeze(Lz(1,2,:)), 'o', t, squeeze(Lz(2,2,:)), 'r', t, squeeze(Lz(3,2,:)), 'g'); grid on;
title('Lz col2 rows 1-3'); xlabel('Time (s)'); ylabel('Lz2');
subplot(232)
plot(t, squeeze(Lz(4,1,:)), 'o', t, squeeze(Lz(5,1,:)), 'r'); grid on;
title('Lz col1 rows 4-5'); xlabel('Time (s)'); ylabel('Lzd1');
subplot(235)
plot(t, squeeze(Lz(4,2,:)), 'o', t, squeeze(Lz(5,2,:)), 'r'); grid on;
title('Lz col2 rows 4-5'); xlabel('Time (s)'); ylabel('Lzd2');
subplot(233)
plot(t, squeeze(Lz(6,1,:)), 'o', t, squeeze(Lz(7,1,:)), 'r'); grid on;
title('Lz col1 rows 6-7'); xlabel('Time (s)'); ylabel('Laz1');
subplot(236)
plot(t, squeeze(Lz(6,2,:)), 'o', t, squeeze(Lz(7,2,:)), 'r'); grid on;
title('Lz col2 rows 6-7'); xlabel('Time (s)'); ylabel('Laz2');

%% ===== Plot: Forecast error covariance Pfz =====
figure;
subplot(311)
plot(t, squeeze(Pfz_k(1,1,:)), 'o', t, squeeze(Pfz_k(2,2,:)), 'r', t, squeeze(Pfz_k(3,3,:)), 'g'); grid on;
title('Pfz diagonal (dx)'); xlabel('Time (s)'); ylabel('Pfz_dx');
subplot(312)
plot(t, squeeze(Pfz_k(4,4,:)), 'o', t, squeeze(Pfz_k(5,5,:)), 'r'); grid on;
title('Pfz diagonal (xd)'); xlabel('Time (s)'); ylabel('Pfz_xd');
subplot(313)
plot(t, squeeze(Pfz_k(6,6,:)), 'o', t, squeeze(Pfz_k(7,7,:)), 'r'); grid on;
title('Pfz diagonal (az)'); xlabel('Time (s)'); ylabel('Pfz_az');

%% ===== Plot: Estimated states =====
figure;
subplot(311)
plot(t, dz1_hat_k, 'r'); grid on;
title('dz1\_hat\_k'); xlabel('Time (s)'); ylabel('dz1 (um)');
subplot(312)
plot(t, dz2_hat_k, 'r'); grid on;
title('dz2\_hat\_k'); xlabel('Time (s)'); ylabel('dz2 (um)');
subplot(313)
plot(t, dz3_hat_k, 'r'); grid on;
title('dz3\_hat\_k'); xlabel('Time (s)'); ylabel('dz3 (um)');

figure;
plot(t, dz1_hat_k, 'b', t, dz2_hat_k, 'r', t, dz3_hat_k, 'g'); grid on;
title('Estimated states'); xlabel('Time (s)'); ylabel('dz (um)');

%% ===== Plot: 3D trajectory =====
figure;
plot3(p_x_data, p_y_data, p_z_data, '.-'); hold on;
plot3(x_d_data, y_d_data, z_d_data, 'r--');
grid on; axis equal;
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)');
title('3D Trajectory');
legend('actual', 'desired');
view(3); rotate3d on;
