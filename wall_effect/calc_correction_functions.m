function [c_par, c_perp] = calc_correction_functions(h_bar)
% CALC_CORRECTION_FUNCTIONS  Wall drag correction factors.
%
%   [c_par, c_perp] = calc_correction_functions(h_bar)
%
%   Inputs:
%     h_bar  - Normalized gap distance h/R (scalar)
%
%   Outputs:
%     c_par  - Parallel correction factor (>= 1)
%     c_perp - Perpendicular correction factor (>= 1)

    invh = 1 ./ h_bar;
    invh3 = invh.^3;
    invh4 = invh.^4;
    invh5 = invh.^5;

    % c_parallel(h_bar)
    denom_par = 1 ...
                - (9/16)   * invh ...
                + (1/8)    * invh3 ...
                - (45/256) * invh4 ...
                - (1/16)   * invh5;
    c_par = 1 ./ denom_par;

    % c_perp(h_bar)
    invh11 = invh.^11;
    invh12 = invh.^12;
    denom_perp = 1 ...
                 - (9/8)    * invh ...
                 + (1/2)    * invh3 ...
                 - (57/100) * invh4 ...
                 + (1/5)    * invh5 ...
                 + (7/200)  * invh11 ...
                 - (1/25)   * invh12;
    c_perp = 1 ./ denom_perp;
end
