function phase_offset = apply_phase_errors(virtual_pos_i, t, theta, beta, kappa)
% APPLY_PHASE_ERRORS Generates residual phase offset for a single element
%
% Inputs:
%   virtual_pos_i - scalar position of this virtual element [m]
%   t             - time vector [1 x n_samples]
%   theta         - angle of arrival [rad]
%   beta          - normalized CFO (scalar, per unit)
%   kappa         - TO in seconds (scalar, per unit)
%
% Output:
%   phase_offset  - 1 x n_samples vector of phase error (radians)

% Constants
fn = 78e9;
lambda = 3e8 / fn;
c = 3e8;

% Default zero offset
phase_offset = zeros(1, length(t));

% Apply only to right-side elements
if virtual_pos_i > 0
    % Generate element-specific localization error
    delta_pos = 1e-3 * randn();  % e.g., 1 mm standard deviation

    % Phase components
    phi_CFO = 2 * pi * fn * beta * t;                           % time-varying
    phi_TO = -2 * pi * fn * kappa;                              % constant
    phi_geom = -2 * pi * fn / c * delta_pos * sin(theta);       % constant

    phase_offset = phi_CFO + phi_TO + phi_geom;
end
end
