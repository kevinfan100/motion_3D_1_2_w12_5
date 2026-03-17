function wall = calc_wall_params(cfg, constants)
% CALC_WALL_PARAMS  Compute wall geometry vectors from config.
%
%   wall = calc_wall_params(cfg, constants)
%
%   Inputs:
%     cfg       - user_config struct (needs theta, phi, pz)
%     constants - physical_constants struct (needs R)
%
%   Outputs:
%     wall - struct with fields:
%       theta, phi, pz, R, w_hat[3x1], u_hat[3x1], v_hat[3x1], W[3x3]

    theta = cfg.theta;
    phi   = cfg.phi;
    pz    = cfg.pz;

    % Plane normal
    wall.w_hat = [cos(theta)*sin(phi);
                  sin(theta)*sin(phi);
                  cos(phi)];

    % Plane tangent vectors
    wall.u_hat = [-cos(theta)*cos(phi);
                  -sin(theta)*cos(phi);
                  sin(phi)];
    wall.v_hat = [sin(theta);
                  -cos(theta);
                  0];

    % Outer-product matrix
    wall.W = wall.w_hat * wall.w_hat.';

    % Store scalar params
    wall.theta = theta;
    wall.phi   = phi;
    wall.pz    = pz;
    wall.R     = constants.R;
end
