function phase_offset = apply_phase_errors(virtual_pos_i, t, theta, beta, kappa)
% Applies clock frequency offset and time offset errors to a signal
% APPLY_PHASE_ERRORS Generates residual phase offset for a single element
%
% Inputs:
%   virtual_pos_i - scalar position of this virtual element [m]
%   t             - time vector [1 x n_samples]
%   theta         - angle of arrival [rad]
%   beta          - normalized CFO (scalar, per unit)
%   kappa         - TO in seconds (scalar, per unit)
%   delta_pos     - optional positional error for this element [m]
%
% Output:
%   phase_offset  - 1 x n_samples vector of phase error (radians)

% Constants
P = uwb_params();
fn = P.fc;
lambda = P.lambda;
c = P.c;

% Default zero offset
phase_offset = zeros(1, length(t));

% Apply Localization error 
delta_pos = 1e-3 * randn();  % 1 mm std
% Phase error due to small position mismatch
phi_geom = -2 * pi * fn / c * delta_pos * sin(theta);  % constant offset

phase_offset = phase_offset + phi_geom;


% Apply CFO and TO only to right-side elements
if virtual_pos_i > 0    
    % Phase components
    phi_CFO = 2 * pi * fn * beta * t;                           % time-varying
    phi_TO = -2 * pi * fn * kappa;                              % constant
    % Sum CFO and TO contributions for this element

    phase_offset = phi_CFO + phi_TO;
end
end
