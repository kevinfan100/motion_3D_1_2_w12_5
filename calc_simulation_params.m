function params = calc_simulation_params(cfg_overrides)
% CALC_SIMULATION_PARAMS  Build hierarchical parameter struct + Bus objects.
%
%   params = calc_simulation_params()
%   params = calc_simulation_params(cfg_overrides)
%
%   Inputs:
%     cfg_overrides - (optional) struct with fields to override in user_config
%
%   Outputs:
%     params - Hierarchical struct:
%       .common   - physical constants (R, gammaN, Ts, kb, T)
%       .wall     - wall geometry (theta, phi, pz, w_hat, u_hat, v_hat, W, R)
%       .ctrl     - controller parameters
%       .thermal  - thermal force parameters
%
%   Side effect: Creates ParametersBus in base workspace.

    % --- Load defaults ---
    constants = physical_constants();
    cfg = user_config();

    % --- Apply overrides ---
    if nargin > 0 && ~isempty(cfg_overrides)
        fnames = fieldnames(cfg_overrides);
        for i = 1:numel(fnames)
            cfg.(fnames{i}) = cfg_overrides.(fnames{i});
        end
    end

    % --- Build hierarchical struct ---
    params.common  = constants;
    params.wall    = calc_wall_params(cfg, constants);
    params.ctrl    = calc_ctrl_params(cfg, constants);
    params.thermal = calc_thermal_params(cfg, constants);

    % --- Also build flat ParametersBus for backward compatibility ---
    build_parameters_bus(cfg, constants);
end

function build_parameters_bus(cfg, c)
% Build the flat 34-element ParametersBus for the Simulink model.
    fields = { ...
        'Ts',           c.Ts;
        'lamdaC',       cfg.lamdaC;
        'theta',        cfg.theta;
        'phi',          cfg.phi;
        'pz',           cfg.pz;
        'R',            c.R;
        'kb',           c.kb;
        'T',            c.T;
        'ax_normal',    cfg.ax_normal;
        'az_normal',    cfg.az_normal;
        'gammaN',       c.gammaN;
        'Avar',         cfg.Avar;
        'Avar2',        cfg.Avar2;
        'Avar22',       cfg.Avar22;
        'Avar3',        cfg.Avar3;
        'Am_scaling',   cfg.Am_scaling;
        'beta',         cfg.beta;
        'lamdaF',       cfg.lamdaF;
        'Pfz_11',       cfg.Pfz_11;
        'Pfz_22',       cfg.Pfz_22;
        'Pfz_33',       cfg.Pfz_33;
        'Pfz_44',       cfg.Pfz_44;
        'Pfz_55',       cfg.Pfz_55;
        'Pfz_66',       cfg.Pfz_66;
        'Pfz_77',       cfg.Pfz_77;
        'rz11_scaling', cfg.rz11_scaling;
        'rz22_scaling', cfg.rz22_scaling;
        'qz11_scaling', cfg.qz11_scaling;
        'qz22_scaling', cfg.qz22_scaling;
        'qz33_scaling', cfg.qz33_scaling;
        'qz44_scaling', cfg.qz44_scaling;
        'qz55_scaling', cfg.qz55_scaling;
        'qz66_scaling', cfg.qz66_scaling;
        'qz77_scaling', cfg.qz77_scaling;
    };

    n = size(fields, 1);
    clear elms;
    for i = 1:n
        elms(i) = Simulink.BusElement; %#ok<AGROW>
        elms(i).Name = fields{i, 1};
        elms(i).DataType = 'double';
        elms(i).Dimensions = 1;
    end

    ParametersBus = Simulink.Bus;
    ParametersBus.Elements = elms;
    assignin('base', 'ParametersBus', ParametersBus);
end
