function ctrl = calc_ctrl_params(cfg, constants)
% CALC_CTRL_PARAMS  Derived controller parameters.
%
%   ctrl = calc_ctrl_params(cfg, constants)
%
%   Inputs:
%     cfg       - user_config struct
%     constants - physical_constants struct
%
%   Outputs:
%     ctrl - struct with all controller tuning parameters

    ctrl.lamdaC     = cfg.lamdaC;
    ctrl.Avar       = cfg.Avar;
    ctrl.Avar2      = cfg.Avar2;
    ctrl.Avar22     = cfg.Avar22;
    ctrl.Avar3      = cfg.Avar3;
    ctrl.Am_scaling = cfg.Am_scaling;
    ctrl.beta       = cfg.beta;
    ctrl.lamdaF     = cfg.lamdaF;

    % Initial Pfz diagonal
    ctrl.Pfz_11 = cfg.Pfz_11;
    ctrl.Pfz_22 = cfg.Pfz_22;
    ctrl.Pfz_33 = cfg.Pfz_33;
    ctrl.Pfz_44 = cfg.Pfz_44;
    ctrl.Pfz_55 = cfg.Pfz_55;
    ctrl.Pfz_66 = cfg.Pfz_66;
    ctrl.Pfz_77 = cfg.Pfz_77;

    % EKF R-matrix scaling
    ctrl.rz11_scaling = cfg.rz11_scaling;
    ctrl.rz22_scaling = cfg.rz22_scaling;

    % EKF Q-matrix scaling
    ctrl.qz11_scaling = cfg.qz11_scaling;
    ctrl.qz22_scaling = cfg.qz22_scaling;
    ctrl.qz33_scaling = cfg.qz33_scaling;
    ctrl.qz44_scaling = cfg.qz44_scaling;
    ctrl.qz55_scaling = cfg.qz55_scaling;
    ctrl.qz66_scaling = cfg.qz66_scaling;
    ctrl.qz77_scaling = cfg.qz77_scaling;

    % Measurement noise variances
    ctrl.meas_noise_var_xy = cfg.meas_noise_var_xy;
    ctrl.meas_noise_var_z  = cfg.meas_noise_var_z;

    % Pre-computed constants
    ctrl.ax_normal = cfg.ax_normal;
    ctrl.az_normal = cfg.az_normal;
end
