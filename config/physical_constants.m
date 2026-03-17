function c = physical_constants()
% PHYSICAL_CONSTANTS  Fixed physical constants for the microprobe system.
%
%   c = physical_constants()
%
%   Returns a struct with:
%     R      - Probe radius [um]
%     gammaN - Nominal Stokes drag coefficient [pN*s/um]
%     Ts     - Sampling period [s]  (1/1600)
%     kb     - Boltzmann constant [J/K]
%     T      - Temperature [K]  (37 C = 310.15 K)

    c.R      = 2.25;                   % [um]
    c.gammaN = 0.0425;                 % [pN*s/um]
    c.Ts     = 1/1600;                 % [s]
    c.kb     = 1.3806503e-23;          % [J/K]
    c.T      = 273.15 + 37;            % [K]
end
