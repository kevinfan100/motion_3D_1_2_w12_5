function [p_d, del_pd] = trajectory_generator(t, params)
% TRAJECTORY_GENERATOR  Generate desired trajectory.
%
%   [p_d, del_pd] = trajectory_generator(t, params)
%
%   Inputs:
%     t      - Current time [s]
%     params - Parameter struct (unused for now, reserved for future)
%
%   Outputs:
%     p_d    - Desired position [x_d; y_d; z_d] [um]
%     del_pd - Desired velocity (unused, set to zeros) [3x1]

    rx = 0;       % [um]
    ry = 0;       % [um]
    rz = 6.0;     % [um]

    theta_x = 2*pi*5*t;     % 5 Hz
    theta_y = 2*pi*5*t;     % 5 Hz

    x_d = rx * sin(theta_x);
    y_d = ry * sin(theta_y);

    if t < 1.5
        z_d = 0;
    elseif t < 2.75
        theta_z = 2*pi*0.4*(t - 1.5);   % 0.4 Hz
        z_d = -rz + rz * cos(theta_z);
    else
        z_d = -2*rz;
    end

    p_d    = [x_d; y_d; z_d];
    del_pd = [0; 0; 0];
end
