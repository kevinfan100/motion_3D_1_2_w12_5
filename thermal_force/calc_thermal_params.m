function th = calc_thermal_params(cfg, constants)
% CALC_THERMAL_PARAMS  Pre-compute thermal force parameters.
%
%   th = calc_thermal_params(cfg, constants)
%
%   Inputs:
%     cfg       - user_config struct
%     constants - physical_constants struct
%
%   Outputs:
%     th - struct with fields:
%       enable          - 1/0 thermal force on/off
%       variance_coeff  - 4*kb*T*gammaN_SI/Ts [N^2]

    th.enable = cfg.thermal_enable;
    th.variance_coeff = 4 * constants.kb * constants.T ...
                        * (constants.gammaN * 1e-6) / constants.Ts;
end
