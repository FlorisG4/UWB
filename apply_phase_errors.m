function total_phase = apply_phase_errors(phase_shift, virtual_pos, t, varargin)
% COMPUTE_PHASE_TAGLIAFERRI Adds Tagliaferri-style CFO, TO, and position error
% Inputs:
%   phase_shift     - Ideal geometric phase [n_el x 1 or n_el x n_samples]
%   virtual_pos     - Virtual positions [n_el x 1] (m)
%   t               - Time vector [1 x n_samples] (s)
%
% Name-Value Optional Inputs:
%   'fn'            - Carrier frequency [Hz] (default 78e9)
%   'beta'          - Normalized CFO per element (Hz/Hz, i.e., Δf/fn) [n_el x 1]
%   'kappa'         - TO (timing offset) per element [s] [n_el x 1]
%   'delta_pos'     - Position errors per element [m] [n_el x 1]
%   'theta'         - Angle of arrival [rad] (default 0)
%
% Output:
%   total_phase     - Phase shift including errors [n_el x n_samples]

% === Parse Inputs ===
p = inputParser;
addParameter(p, 'fn', 78e9);
addParameter(p, 'beta', zeros(size(virtual_pos)));
addParameter(p, 'kappa', zeros(size(virtual_pos)));
addParameter(p, 'delta_pos', zeros(size(virtual_pos)));
addParameter(p, 'theta', 0);  % angle of arrival
parse(p, varargin{:});

fn = p.Results.fn;
beta = p.Results.beta(:);
kappa = p.Results.kappa(:);
delta_pos = p.Results.delta_pos(:);
theta = p.Results.theta;
lambda = 3e8 / fn;
c = 3e8;

% === Matrix Expansion ===
n_el = length(virtual_pos);
n_samples = length(t);
t_mat = repmat(t, n_el, 1);

% === Tagliaferri 2024 Phase Terms ===
% 1. Carrier Frequency Offset (CFO): α_nm = 2π f_n β_nm t
phi_CFO = 2 * pi * fn * beta .* t_mat;

% 2. Time Offset (TO): 2π f_n * kappa
phi_TO = -2 * pi * fn * repmat(kappa, 1, n_samples);

% 3. Localization Error: Δpos * sin(θ)
phi_geom = -2 * pi * fn / c * (delta_pos .* sin(theta));
phi_geom = repmat(phi_geom, 1, n_samples);

% === Total Phase ===
total_phase = phase_shift + phi_CFO + phi_TO + phi_geom;

end
