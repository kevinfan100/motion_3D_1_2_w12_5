function [Gamma_inv, h_bar] = calc_gamma_inv(p, params)
% CALC_GAMMA_INV  3x3 inverse mobility matrix with wall effect.
%
%   [Gamma_inv, h_bar] = calc_gamma_inv(p, params)
%
%   Inputs:
%     p      - Position [3x1] [um]
%     params - Parameter struct (needs wall.w_hat, wall.pz, wall.R, wall.W,
%              common.gammaN)
%
%   Outputs:
%     Gamma_inv - 3x3 inverse drag matrix [um/(pN*s)]
%     h_bar     - Normalized gap distance h/R

    w_hat  = params.wall.w_hat;
    pz     = params.wall.pz;
    R      = params.wall.R;
    gammaN = params.common.gammaN;
    W      = params.wall.W;

    % Gap distance
    h     = dot(p, w_hat) - pz;
    h_bar = h / R;

    % Correction functions
    [c_par, c_perp] = calc_correction_functions(h_bar);

    % Inverse mobility matrix
    I3    = eye(3);
    coeff = 1 / (gammaN * c_par);
    alpha = (c_perp - c_par) / c_perp;
    Gamma_inv = coeff * (I3 - alpha * W);
end
