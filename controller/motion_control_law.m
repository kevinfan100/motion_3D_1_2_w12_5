function [fd, dzm_bar_k, dz1_hat_k, dz2_hat_k, dz3_hat_k, dz_k2, ...
          dzr_k, azm_k, az1_hat_k, mgain_z, Lz, Pfz_k] = ...
          motion_control_law(t, xd_k, yd_k, zd_k, p_k, mgain, params)
% MOTION_CONTROL_LAW  IIR variance filter + Eq.13 + 7-state EKF + control.
%
%   [fd, ...] = motion_control_law(t, xd_k, yd_k, zd_k, p_k, mgain, params)
%
%   Inputs:
%     t      - Current time [s]
%     xd_k   - Desired x position [um]
%     yd_k   - Desired y position [um]
%     zd_k   - Desired z position [um]
%     p_k    - 2-step delayed measured position [3x1] [um]
%     mgain  - Motion gain from wall effect [mgain_x; mgain_y; mgain_z] [um/pN]
%     params - Flat Parameters struct (ParametersBus)
%
%   Outputs:
%     fd         - Control force [3x1] [pN]
%     dzm_bar_k  - Deterministic component of z tracking error
%     dz1_hat_k  - EKF estimated state 1
%     dz2_hat_k  - EKF estimated state 2
%     dz3_hat_k  - EKF estimated state 3
%     dz_k2      - z tracking error (2-step delayed)
%     dzr_k      - Random component of z tracking error
%     azm_k      - Estimated motion gain (IIR, scaled)
%     az1_hat_k  - EKF estimated motion gain
%     mgain_z    - z-axis motion gain (pass-through)
%     Lz         - EKF Kalman gain [7x2]
%     Pfz_k      - EKF forecast error covariance [7x7]

% ====================================================================
% PERSISTENT STATE VARIABLES
% ====================================================================
persistent xd_k1 xd_k2 yd_k1 yd_k2 zd_k1 zd_k2
persistent fd_X fd_Y fd_Z

persistent dx1_hat_k1 dx2_hat_k1 dx3_hat_k1
persistent dz1_hat_k1 dz2_hat_k1 dz3_hat_k1
persistent xD_hat_k1 dxD_hat_k1
persistent zD1_hat_k1 zD2_hat_k1
persistent ax_hat_k1 ax_hat_k2
persistent az1_hat_k1 az1_hat_k2
persistent dax_hat_k1 az2_hat_k1
persistent Pfx_k1 Pfz_k1
persistent fd_X_k1 fd_X_k2 fd_Z_k1 fd_Z_k2

persistent dxm_bar_k1 dxr_bar_k1 dxr2_bar_k1
persistent dzm_bar_k1 dzr_bar_k1 dzr2_bar_k1
persistent dxrr_bar_k1 dzrr_bar_k1
persistent dxm_k1 dzm_k1 azm_k1

persistent dx_k3 dy_k3 dz_k3

% ====================================================================
% INITIALIZATION
% ====================================================================
if isempty(xd_k1)
    xd_k1 = 0; xd_k2 = 0;
    yd_k1 = 0; yd_k2 = 0;
    zd_k1 = 0; zd_k2 = 0;
    fd_X = 0; fd_Y = 0; fd_Z = 0;

    dx1_hat_k1 = 0; dx2_hat_k1 = 0; dx3_hat_k1 = 0;
    dz1_hat_k1 = 0; dz2_hat_k1 = 0; dz3_hat_k1 = 0;
    xD_hat_k1 = 0; dxD_hat_k1 = 0;
    zD1_hat_k1 = 0; zD2_hat_k1 = 0;

    ax_hat_k1  = 0.012;
    ax_hat_k2  = 0;
    az1_hat_k1 = 0.012 / (1 + params.beta);
    az1_hat_k2 = 0;
    dax_hat_k1 = 0;
    az2_hat_k1 = 0;

    Pfx_k1 = 0.0147 * 4 * 1.3806503e-23 * (273.15+37) * diag(ones(1,7)) * 1/4 * 1e-2;

    Pfz_k1 = diag([params.Pfz_11, params.Pfz_22, params.Pfz_33, ...
                    params.Pfz_44, params.Pfz_55, params.Pfz_66, params.Pfz_77]);

    dxm_bar_k1  = 0; dzm_bar_k1  = 0;
    dxr_bar_k1  = 0; dzr_bar_k1  = 0;
    dxrr_bar_k1 = 0; dzrr_bar_k1 = 0;
    dxr2_bar_k1 = 0.9e-15;
    dzr2_bar_k1 = (0.012e6) / params.Am_scaling * 4 * params.kb * params.T ...
                  * (2 + 1/(1 - params.lamdaC^2)) * (1e6)^2;

    fd_X_k1 = 0; fd_X_k2 = 0;
    fd_Z_k1 = 0; fd_Z_k2 = 0;

    dxm_k1 = 0; dzm_k1 = 0;
    azm_k1 = 0.012 / params.Am_scaling;

    dx_k3 = 0; dy_k3 = 0; dz_k3 = 0;
end

% ====================================================================
% EXTRACT PARAMETERS
% ====================================================================
Ts     = params.Ts;
lamdaC = params.lamdaC;
kb     = params.kb;
T      = params.T;

Avar   = params.Avar;
Avar2  = params.Avar2;
Avar22 = params.Avar22;
Avar3  = params.Avar3;

mgain_x = mgain(1);
mgain_y = mgain(2);
mgain_z_val = mgain(3);

% ====================================================================
% MEASUREMENT WITH NOISE
% ====================================================================
Measurement_Variance_x_y = (0.7e-9)^2;      % [m^2]
Measurement_Variance_z   = (2.3e-9)^2;      % [m^2]

px_k = p_k(1) + sqrt(Measurement_Variance_x_y) * randn(1,1) * 1e6;   % [um]
py_k = p_k(2) + sqrt(Measurement_Variance_x_y) * randn(1,1) * 1e6;   % [um]
pz_k = p_k(3) + sqrt(Measurement_Variance_z)   * randn(1,1) * 1e6;   % [um]

% ====================================================================
% TRACKING ERROR (2-step delayed)
% ====================================================================
dx_k2 = xd_k2 - px_k;
dy_k2 = yd_k2 - py_k;
dz_k2 = zd_k2 - pz_k;

dxm_k = dx_k2;
dzm_k = dz_k2;

% ====================================================================
% IIR VARIANCE FILTER — X AXIS (Eqs. 9-10)
% ====================================================================
dxm_bar_k = Avar * dxm_k1 + (1-Avar) * dxm_bar_k1;
dxr_k     = dxm_k - dxm_bar_k;
dxr_bar_k = Avar2 * dxr_k + (1-Avar2) * dxr_bar_k1;
dxrr_bar_k = Avar22 * dxr_k + (1-Avar22) * dxrr_bar_k1;
dxr2_bar_k = Avar3 * dxr_k^2 + (1-Avar3) * dxr2_bar_k1;
sigma2_dxr = dxr2_bar_k - dxrr_bar_k^2;
if sigma2_dxr < 0, sigma2_dxr = 0; end

% ====================================================================
% IIR VARIANCE FILTER — Z AXIS (Eqs. 9-10)
% ====================================================================
dzm_bar_k = Avar * dzm_k + (1-Avar) * dzm_bar_k1;
dzr_k     = dzm_k1 - dzm_bar_k;
dzr_k     = min(max(dzr_k, -6e-2), 6e-2);
dzr_bar_k = Avar2 * dzr_k + (1-Avar2) * dzr_bar_k1;
dzrr_bar_k = Avar22 * dzr_k + (1-Avar22) * dzrr_bar_k1;
dzr2_bar_k = Avar3 * dzr_k^2 + (1-Avar3) * dzr2_bar_k1;
sigma2_dzr = dzr2_bar_k - dzrr_bar_k^2;
if sigma2_dzr < 0, sigma2_dzr = 0; end

% ====================================================================
% EQ. 13: MOTION GAIN ESTIMATION
% ====================================================================
sigma2_nx = (0.7e-3)^2;     % [um^2]
sigma2_nz = (2.3e-3)^2;     % [um^2]

den   = 4 * kb * T * (2 + 1/(1 - lamdaC^2));   % [N*m]
num_x = sigma2_dxr - (2/(1 + lamdaC)) * sigma2_nx;
num_z = sigma2_dzr - (2/(1 + lamdaC)) * sigma2_nz;
if num_x < 0, num_x = 0; end
if num_z < 0, num_z = 0; end

axm_k_dyna = (num_x * 1e-12) / den * (1e6/1e12);   % [um/(pN*s)]
azm_k_dyna = (num_z * 1e-12) / den * (1e6/1e12);   % [um/(pN*s)]

axm_k = 2.5 * axm_k_dyna;
azm_k = params.Am_scaling * azm_k_dyna;

% ====================================================================
% 7-STATE EKF — X AXIS
% ====================================================================
lamdaF = params.lamdaF;

e_x1_k = dxm_bar_k - dx1_hat_k1;
e_ax_k = axm_k - ax_hat_k1;
esti_error_x = [e_x1_k; e_ax_k];
H = [1 0 0 0 0 0 0; 0 0 0 0 0 1 0];

r11 = (0.7e-3)^2 * Avar;
r22 = (1e-3 * ax_hat_k1)^2;
RR = diag([r11, r22]);

q33 = ax_hat_k1 * 4 * kb * T * 1e6 * 1e12;
tau_az = 100.0;
K_tau  = 1 - exp(-Ts/tau_az);
betaQ  = K_tau^2;
q66 = betaQ * r22;
q77 = 10 * betaQ * r22;
Q_x = diag([0, 0, q33, 0, 0, q66, q77]);

Fe_x = [0 1 0  0 0      0         0;
        0 0 1  0 0      0         0;
        0 0 1 -1 0  -fd_X_k1     0;
        0 0 0  1 1      0         0;
        0 0 0  0 1      0         0;
        0 0 0  0 0      1         1;
        0 0 0  0 0      0         1];

S = H * Pfx_k1 * H.' + RR;
S = 0.5*(S + S.') + 1e-20*eye(size(S));
Lx = (Pfx_k1 * H.') / S;
Px = (eye(7) - Lx*H) * Pfx_k1;
Px = 0.5*(Px + Px.');
Pfx_k = Fe_x * Px * Fe_x.' + Q_x;
Pfx_k = 0.5*(Pfx_k + Pfx_k.');

inj_x = Lx * esti_error_x;
dx1_hat_k = dx2_hat_k1 + inj_x(1);
dx2_hat_k = dx3_hat_k1 + inj_x(2);
dx3_hat_k = lamdaC * dx3_hat_k1 + inj_x(3);
xD_hat_k  = xD_hat_k1 + dxD_hat_k1 + inj_x(4);
dxD_hat_k = dxD_hat_k1 + inj_x(5);
ax_hat_k  = ax_hat_k1 + dax_hat_k1 + inj_x(6);
dax_hat_k = dax_hat_k1 + inj_x(7);

% ====================================================================
% 7-STATE EKF — Z AXIS (with beta w1-w2)
% ====================================================================
e_z1_k = dz_k2 - dz1_hat_k1;
e_az_k = azm_k1 - az1_hat_k1;
esti_error_z = [e_z1_k; e_az_k];

beta_val = params.beta;

r11z = (az1_hat_k1 * 1e6) * 4 * kb * T * 1e12 * params.rz11_scaling;
r22z = (params.rz22_scaling * az1_hat_k1)^2;
RR_z = diag([r11z, r22z]);

q11z = (az1_hat_k1*1e6)*4*kb*T*1e12 * params.qz11_scaling;
q22z = (az1_hat_k1*1e6)*4*kb*T*1e12 * params.qz22_scaling;
q33z = (az1_hat_k1*1e6)*4*kb*T*1e12 * params.qz33_scaling;
q44z = (az1_hat_k1*1e6)*4*kb*T*1e12 * params.qz44_scaling;
q55z = (az1_hat_k1*1e6)*4*kb*T*1e12 * params.qz55_scaling;
q66z = params.qz66_scaling * r22z;
q77z = params.qz77_scaling * r22z;
Q_z = diag([q11z, q22z, q33z, q44z, q55z, q66z, q77z]);

Fe_z = [0 1 0  0          0         0          0;
        0 0 1  0          0         0          0;
        0 0 1 -1          0     -fd_Z_k1       0;
        0 0 0  1+beta_val -beta_val  0          0;
        0 0 0  1          0         0          0;
        0 0 0  0          0     1+beta_val  -beta_val;
        0 0 0  0          0         1          0];

S_z = H * Pfz_k1 * H.' + RR_z;
S_z = 0.5*(S_z + S_z.') + 1e-20*eye(size(S_z));
Lz = (Pfz_k1 * H.') / S_z;
Pz = 1/lamdaF * (eye(7) - Lz*H) * Pfz_k1;
Pz = 0.5*(Pz + Pz.');
Pfz_k = Fe_z * Pz * Fe_z.' + Q_z;
Pfz_k = 0.5*(Pfz_k + Pfz_k.');

inj_z = Lz * esti_error_z;
dz1_hat_k = dz2_hat_k1 + inj_z(1);
dz2_hat_k = dz3_hat_k1 + inj_z(2);
dz3_hat_k = lamdaC * dz3_hat_k1 + inj_z(3);
zD1_hat_k = (1+beta_val)*zD1_hat_k1 - beta_val*zD2_hat_k1 + inj_z(4);
zD2_hat_k = zD1_hat_k1 + inj_z(5);
az1_hat_k = (1+beta_val)*az1_hat_k1 - beta_val*az2_hat_k1 + inj_z(6);
az2_hat_k = az1_hat_k1 + inj_z(7);

% ====================================================================
% CONTROL FORCE
% ====================================================================
% z-axis: EKF-based control
fd_Z = 1/az1_hat_k * (zd_k - zd_k1 + (1-lamdaC)*dz3_hat_k - zD1_hat_k1);

% x,y-axis: conventional control with wall-corrected gain
fd_X = (1/mgain_x) * (xd_k - xd_k1 + (1-lamdaC)*dx_k3);
fd_Y = (1/mgain_y) * (yd_k - yd_k1 + (1-lamdaC)*dy_k3);

fd = [fd_X; fd_Y; fd_Z];
mgain_z = mgain_z_val;

% ====================================================================
% STATE UPDATES
% ====================================================================
xd_k2 = xd_k1;  xd_k1 = xd_k;
yd_k2 = yd_k1;  yd_k1 = yd_k;
zd_k2 = zd_k1;  zd_k1 = zd_k;

dx1_hat_k1 = dx1_hat_k; dx2_hat_k1 = dx2_hat_k; dx3_hat_k1 = dx3_hat_k;
dz1_hat_k1 = dz1_hat_k; dz2_hat_k1 = dz2_hat_k; dz3_hat_k1 = dz3_hat_k;
xD_hat_k1  = xD_hat_k;  dxD_hat_k1 = dxD_hat_k;
zD1_hat_k1 = zD1_hat_k; zD2_hat_k1 = zD2_hat_k;
ax_hat_k2  = ax_hat_k1; ax_hat_k1  = ax_hat_k;
az1_hat_k2 = az1_hat_k1; az1_hat_k1 = az1_hat_k;
dax_hat_k1 = dax_hat_k; az2_hat_k1 = az2_hat_k;
Pfx_k1 = Pfx_k; Pfz_k1 = Pfz_k;

dxm_bar_k1  = dxm_bar_k;  dzm_bar_k1  = dzm_bar_k;
dxr_bar_k1  = dxr_bar_k;  dzr_bar_k1  = dzr_bar_k;
dxrr_bar_k1 = dxrr_bar_k; dzrr_bar_k1 = dzrr_bar_k;
dxr2_bar_k1 = dxr2_bar_k; dzr2_bar_k1 = dzr2_bar_k;

fd_X_k2 = fd_X_k1; fd_X_k1 = fd_X;
fd_Z_k2 = fd_Z_k1; fd_Z_k1 = fd_Z;

dx_k3 = dx_k2;
dy_k3 = dy_k2;
dz_k3 = dz_k2;

dxm_k1 = dxm_k;
dzm_k1 = dzm_k;
azm_k1 = azm_k;

end
