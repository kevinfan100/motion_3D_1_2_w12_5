function cfg = user_config()
% USER_CONFIG  All user-adjustable parameters for the simulation.
%
%   cfg = user_config()
%
%   Returns a struct with fields grouped by subsystem.

    % --- Controller bandwidth ---
    cfg.lamdaC = 0.7;

    % --- Wall geometry ---
    cfg.theta = pi/2;          % inclined plane angle [rad]
    cfg.phi   = 0;             % inclined plane angle [rad]
    cfg.pz    = -15;           % wall offset [um]

    % --- Nominal motion gain ---
    cfg.ax_normal = 0.0147;    % [um/(pN*s)]
    cfg.az_normal = 0.0147;    % [um/(pN*s)]

    % --- IIR filter for deterministic component (Eqs. 9-10) ---
    cfg.Avar   = 0.45;         % IIR pole for mean tracking error
    cfg.Avar2  = 0.005;        % IIR pole for mean of random component
    cfg.Avar22 = 0.05;         % IIR pole for mean square (numerator)
    cfg.Avar3  = 0.05;         % IIR pole for mean square (denominator)
    cfg.Am_scaling = 10;       % Eq.13 azm scaling factor

    % --- w1-w2 estimator ---
    cfg.beta = 0.5;

    % --- Kalman filter for z ---
    cfg.lamdaF = 1;            % forgetting factor (1 = no forgetting)

    % Initial forecast error covariance Pfz diagonal
    cfg.Pfz_11 = 0;                    % [um^2]
    cfg.Pfz_22 = 0;                    % [um^2]
    cfg.Pfz_33 = 1e-4;                 % [um^2]
    cfg.Pfz_44 = 1e-4;                 % [um^2]
    cfg.Pfz_55 = 0;                    % [um^2]
    cfg.Pfz_66 = 10*(0.0147)^2;        % [um/(pN*s)]^2
    cfg.Pfz_77 = 0;                    % [um/(pN*s)]^2

    % EKF R-matrix scaling
    cfg.rz11_scaling = 1;              % base: (az1_hat*1e6)*4*kb*T*(1e12)
    cfg.rz22_scaling = 1e-1;           % base: (az1_hat)^2

    % EKF Q-matrix scaling
    cfg.qz11_scaling = 0;
    cfg.qz22_scaling = 0;
    cfg.qz33_scaling = 1e4;
    cfg.qz44_scaling = 1e2;
    cfg.qz55_scaling = 0;
    cfg.qz66_scaling = 8e-5;
    cfg.qz77_scaling = 0;

    % --- Thermal force ---
    cfg.thermal_enable = 1;            % 1 = on, 0 = off

    % --- Measurement noise ---
    cfg.meas_noise_var_xy = (0.7e-3)^2;   % [um^2]
    cfg.meas_noise_var_z  = (2.3e-3)^2;   % [um^2]

    % --- Simulation ---
    cfg.StopTime = 3.5;               % [s]
end
