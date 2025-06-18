function P = uwb_params()
%UWB_PARAMS  Central definition of radar constants
%   P = UWB_PARAMS() returns a struct containing commonly used
%   parameters for the UWB radar simulations.

P.c       = 3e8;        % Speed of light [m/s]
P.fc      = 78e9;       % Carrier frequency [Hz]
P.Fs      = 10e6;       % Sampling frequency [Hz]
P.T_chirp = 25.6e-6;    % Chirp duration [s]
P.BW      = 250e6;      % Bandwidth [Hz]

% Derived quantities
P.lambda  = P.c / P.fc; % Wavelength [m]
end