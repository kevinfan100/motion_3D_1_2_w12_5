function f_th = calc_thermal_force(p, params)
% CALC_THERMAL_FORCE  Brownian force with wall-effect correction.
%
%   f_th = calc_thermal_force(p, params)
%
%   Inputs:
%     p      - Position [3x1] [um]
%     params - Parameter struct (needs wall, common, thermal)
%
%   Outputs:
%     f_th   - Thermal force [3x1] [pN]

    w_hat  = params.wall.w_hat;
    u_hat  = params.wall.u_hat;
    v_hat  = params.wall.v_hat;
    pz     = params.wall.pz;
    R      = params.wall.R;
    gammaN = params.common.gammaN;
    kb     = params.common.kb;
    T      = params.common.T;
    Ts     = params.common.Ts;

    % Gap distance and correction factors
    h     = dot(p, w_hat) - pz;
    h_bar = h / R;
    [c_par, c_perp] = calc_correction_functions(h_bar);

    % Direction-dependent correction vector
    C = c_par * (u_hat + v_hat) + c_perp * w_hat;

    % Variance: 4*kb*T*gammaN_SI / Ts * |C|
    variance_coeff = 4 * kb * T * (gammaN * 1e-6) / Ts;   % [N^2]  (MKS)
    Variance = variance_coeff * abs(C);

    % Generate random thermal force
    f_th = sqrt(Variance) .* randn(3, 1) * 1e12;            % [pN]
end
