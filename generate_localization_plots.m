function [] = generate_localization_plots(rx_signals,SNR_dB,Fs,T_chirp,c,BW,n_targets)
% Generate range-angle map plots and estimate peak locations
 % === Add Noise ===
signal_power = var(rx_signals(:));
noise_power = signal_power / 10^(SNR_dB/10);
noise = sqrt(noise_power) * randn(size(rx_signals));
rx_signals = rx_signals + noise;

% === Range FFT ===
N_fft =1024;  % Zero-padding improves frequency resolution
range_fft = fft(rx_signals, N_fft, 2);  % Along time axis (fast-time)
range_fft = range_fft(:, 1:N_fft/2);    % Keep positive frequencies

% Corresponding range bins
range_axis = ((0:N_fft/2-1) * Fs / N_fft) * (c * T_chirp / (2 * BW));

% === Range-Angle Map ===
N_ffta = 1024;
angle_axis = asind(linspace(-1,1,N_ffta));
RAOA = fftshift(fft(range_fft, N_ffta, 1),1);

%Method 1 simple peak detection
[est_ranges, est_angles, peak_values] = localize_targets_peak(RAOA, range_axis, angle_axis, n_targets);

%Localization plot
figure;
pcolor(angle_axis, range_axis, 20*log10(abs(RAOA.')/max(abs(RAOA(:)))));
hold on
for l = 1 : length(est_angles)
plot(est_angles(l), est_ranges(l) , 'm.', 'MarkerSize', 20);
end
shading flat; colormap jet; caxis([-40 0]);
xlabel('Angle (deg)'); ylabel('Range (m)');
title('Range Angle Map with pre-sync Localization');
hold off
