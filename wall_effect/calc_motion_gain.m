function mgain = calc_motion_gain(p, params)
% CALC_MOTION_GAIN  Compute wall-corrected motion gain vector.
%
%   mgain = calc_motion_gain(p, params)
%
%   Inputs:
%     p      - Position [3x1] [um]
%     params - Parameter struct (needs wall, common)
%
%   Outputs:
%     mgain  - Motion gain [mgain_x; mgain_y; mgain_z] [um/pN]

    w_hat  = params.wall.w_hat;
    pz     = params.wall.pz;
    R      = params.wall.R;
    gammaN = params.common.gammaN;
    Ts     = params.common.Ts;

    % Gap distance
    h     = dot(p, w_hat) - pz;
    h_bar = h / R;

    % Correction functions
    [c_par, c_perp] = calc_correction_functions(h_bar);

    % Motion gains
    gamma_x = gammaN * c_par;
    gamma_y = gammaN * c_par;
    gamma_z = gammaN * c_perp;

    mgain = [Ts/gamma_x; Ts/gamma_y; Ts/gamma_z];
end
