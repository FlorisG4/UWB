
% =========================================================================
% UWB MIMO Radar Localization Simulation
% 
% =========================================================================

clear; clc;

%% === Simulation Parameters ===
c = 3e8;                  % Speed of light [m/s]
fc = 78e9;                % Carrier frequency [Hz]
lambda = c / fc;          % Wavelength [m]
Fs = 10e6;                % Sampling frequency [Hz]
T_chirp = 25.6e-6;        % Chirp duration [s]
BW = 250e6;                % Bandwidload th [Hz]
SNR_dB = 20;              % Signal-to-noise ratio
n_chirps = 16;

%% === MIMO Radar Setup ===
N_tx = 3;
N_rx = 4;
d = lambda / 2; %Inter-array spacing in (m)
Baseline = 0.05 ; %Baseline in (m)

% ===  Array Positions === %
[tx_pos,rx_pos,virtual_pos] = create_array(N_tx,N_rx,d,Baseline);

% ===  Visualize the array === %
plot_mimo_array(tx_pos, rx_pos,virtual_pos);

% ===  Compute ff distance for Baseline === %
D= max(virtual_pos) - min(virtual_pos); %aperture size of whole array in m
R_ff= (2*D^2)/lambda;
fprintf('Far-field distance: %.2f meters\n', R_ff);


%% === Target Setup ===
targets = [1500, deg2rad(10);   % [Range (m), Angle (rad)]
           ];

if any(targets(:,1) < R_ff)
    error('One or more targets are in the near field — computation not valid.');
end

n_targets = size(targets,1);

%% === Signal Generation ===
t = 0:1/Fs:T_chirp-1/Fs;
n_samples = length(t);
rx_signals = zeros(length(virtual_pos), n_samples);

for el = 1:length(virtual_pos)
    for tgt = 1:n_targets
        R = targets(tgt, 1);
        theta = targets(tgt, 2);

        fb = 2 * BW * R / (c * T_chirp);  % Beat frequency

        phase_shift = 2 * pi * virtual_pos(el) * sin(theta) / lambda; %Phase shift induced by position of element & AOA

        % delta_f_c = f_tx - f_rx;                                        % The carrier offset introduced by different osicllators
        % phase_shift_CFO = 2 * pi * delta_f_c * t;                         % CFO-induced phase rotation over time
        % phase_shift_TO = -2 * pi * fb * delta_t;                        % TO leads to a beat frequency phase offset
        % phase_shift_geom = 2 * pi / lambda * delta_pos(el) * sin(theta); % localization error per element
        % 
        % phase_shift = phase_shift + phase_shift_CFO + phase_shift_TO + phase_shift_geom;


        signal = exp(1j * (2*pi*fb*t + phase_shift));

        rx_signals(el,:) = rx_signals(el,:) + signal;
    end
end

 %% === Add Noise ===
signal_power = var(rx_signals(:));
noise_power = signal_power / 10^(SNR_dB/10);
noise = sqrt(noise_power) * randn(size(rx_signals));
% rx_signals = rx_signals + noise;

%% === FFT and Localization Processing ===
% This is where you'd insert range FFT, angle FFT, and optionally Block-FOCUSS
%% === Range FFT ===
N_fft = 512;  % Zero-padding improves frequency resolution
range_fft = fft(rx_signals, N_fft, 2);  % Along time axis (fast-time)
range_fft = range_fft(:, 1:N_fft/2);    % Keep positive frequencies

% Corresponding range bins
range_axis = ((0:N_fft/2-1) * Fs / N_fft) * (c * T_chirp / (2 * BW));




%% === Angle Steering Matrix ===
angle_grid = deg2rad(-30:0.5:30);  % Search angles from -30° to 30°
steering_matrix = zeros(length(virtual_pos), length(angle_grid));

for a = 1:length(angle_grid)
    steering_matrix(:,a) = exp(1j * 2 * pi * virtual_pos * sin(angle_grid(a)) / lambda);
end

%% === Angle Estimation at Strongest Range Bin ===
range_power = sum(abs(range_fft).^2, 1);
[~, max_bin] = max(range_power);

signal_vector = range_fft(:, max_bin);
spectrum = abs(steering_matrix' * signal_vector).^2;

% Normalize and plot
spectrum = spectrum / max(spectrum);
figure;
plot(rad2deg(angle_grid), 10*log10(spectrum));
xlabel('Angle (deg)');
ylabel('Power (dB)');
title('Angle Spectrum at Strongest Range Bin');
grid on;

figure;
plot(0:N_fft/2-1, 10*log10(range_power));
xlabel('Range (m)');
ylabel('Power (dB)');
title('Range Profile');
grid on;


disp('Simulation complete.');
