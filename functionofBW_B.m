% =========================================================================
% UWB MIMO Radar Localization Simulation
% 
% =========================================================================

clear; clc;
%% === Define parameter range ===
baseline_values = linspace(0.05, 0.5 , 10);     % in meters
bw_values = linspace(50e6, 500e6, 10);         % in Hz

rmse_matrix = zeros(length(baseline_values), length(bw_values));

% === Loop over BW and Baseline combinations ===
for i = 1:length(baseline_values)
    for j = 1:length(bw_values)
        Baseline = baseline_values(i);
        BW = bw_values(j);

%% === Simulation Parameters ===
P = uwb_params();
c = P.c;
fc = P.fc;
lambda = P.lambda;
Fs = P.Fs;
T_chirp = P.T_chirp;
SNR_dB = 20;              % Signal-to-noise ratio
n_chirps = 16;

%% === MIMO Radar Setup ===
N_tx = 3;
N_rx = 4;
d = lambda / 2; %Inter-array spacing in (m)

% ===  Array Positions === %
[tx_pos,rx_pos,virtual_pos,~] = create_array(N_tx,N_rx,d,Baseline);

virtual_pos_x_sep = virtual_pos(:,[1,3]);
virtual_pos_x = [virtual_pos_x_sep(:,1); virtual_pos_x_sep(:,2)];  

n_radars = 2;
n_elements_per_radar = length(virtual_pos_x)/n_radars;
radar_indices = reshape(1:length(virtual_pos_x), n_elements_per_radar, n_radars);

% === Initialize storage ===
I = cell(n_radars, n_radars);  % only monostatic for now


%% === Target Setup ===
targets = [40, deg2rad(10);   % [Range (m), Angle (rad)]
            50,deg2rad(-10);
             30, deg2rad(0)];


n_targets = size(targets,1);

%% === Signal Generation ===
t = 0:1/Fs:T_chirp-1/Fs;
n_samples = length(t);
rx_signals = zeros(length(virtual_pos_x), n_samples);

% Parameters for true offsets
beta_true = 100 * randn();  % CFO frequency offset [Hz]
kappa_true = 1e-9 * randn();     % TO parameter [s
beta_values = [];

for el = 1:length(virtual_pos_x)
    rx_signals_el = zeros(1, n_samples);  % for accumulation over targets

    for tgt = 1:n_targets
        R = targets(tgt, 1);
        theta = targets(tgt, 2);

        fb = 2 * BW * R / (c * T_chirp);  % Beat frequency

        % === Ideal signal (no errors)
        phase_shift_ideal = 2 * pi * virtual_pos_x(el) * sin(theta) / lambda;
        signal_ideal = exp(1j * (2*pi*fb*t + phase_shift_ideal));

        % === Corrupted signal
        % phase_offsets = apply_phase_errors(virtual_pos_x(el), t, theta, beta_true/fc, kappa_true); % NOTE: normalized beta
        phase_offsets = 0;
        phase_shift_corrup = phase_shift_ideal + phase_offsets;
        signal_corrup = exp(1j * (2*pi*fb*t + phase_shift_corrup));

        % === Accumulate signals
        rx_signals_el = rx_signals_el + signal_corrup;

        % === Estimate CFO (only for one target and right-side elements)
        if el > 12 && tgt == 1
            phi_residual = unwrap(angle(signal_corrup ./ signal_ideal));
            p = polyfit(t, phi_residual, 1);
            beta_tmp = p(1) / (2 * pi * fc);
            beta_values(end+1) = beta_tmp;
        end
    end

    % Store final signal after all targets added
    rx_signals(el,:) = rx_signals_el;
end

% === Final estimate across all right-side elements
beta_est_norm = mean(beta_values, 'omitnan');
beta_est = beta_est_norm * fc;  % frequency offset in Hz


 %% === Add Noise ===
signal_power = var(rx_signals(:));
noise_power = signal_power / 10^(SNR_dB/10);
noise = sqrt(noise_power) * randn(size(rx_signals));
rx_signals = rx_signals + noise;

%% === Range FFT ===
N_fft =1024;  % Zero-padding improves frequency resolution
range_fft = fft(rx_signals, N_fft, 2);  % Along time axis (fast-time)
range_fft = range_fft(:, 1:N_fft/2);    % Keep positive frequencies

% Corresponding range bins
range_axis = ((0:N_fft/2-1) * Fs / N_fft) * (c * T_chirp / (2 * BW));

%% === Range-Angle Map ===
N_ffta = 1024;
angle_axis = asind(linspace(-1,1,N_ffta));
RAOA = fftshift(fft(range_fft, N_ffta, 1),1);
%% === Localization approach with offset present ===
%Method 1 simple peak detection
[est_ranges, est_angles, peak_values] = localize_targets_peak(RAOA, range_axis, angle_axis, n_targets);

[rmse, FAR, P_det] = evaluation(targets, est_ranges, est_angles, BW, Baseline);

rmse_matrix(i, j) = rmse;

    end
end

% === Gemiddelde RMSE over baselines ===
mean_rmse_per_bw = mean(rmse_matrix, 1);

% === Plot ===
figure;
plot(bw_values/1e6, mean_rmse_per_bw, 'o-','LineWidth',1.5);
xlabel('Bandwidth (MHz)');
ylabel('Mean RMSE');
title('RMSE vs Bandwidth (averaged over Baseline)');
grid on;

% === Gemiddelde RMSE over BW ===
mean_rmse_per_baseline = mean(rmse_matrix, 2);

% === Plot ===
figure;
plot(baseline_values, mean_rmse_per_baseline, 's-','LineWidth',1.5);
xlabel('Baseline (m)');
ylabel('Mean RMSE');
title('RMSE vs Baseline (averaged over Bandwidth)');
grid on;