
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
BW = 250e6;                % Bandwidth [Hz]
SNR_dB = 20;              % Signal-to-noise ratio
n_chirps = 16;

%% === MIMO Radar Setup ===
N_tx = 3;
N_rx = 4;
d = lambda / 2; %Inter-array spacing in (m)
Baseline = 0.05 ; %Baseline in (m)

% ===  Array Positions === %
[tx_pos,rx_pos,virtual_pos,~] = create_array(N_tx,N_rx,d,Baseline);

virtual_pos_x_sep = virtual_pos(:,[1,3]);
virtual_pos_x = [virtual_pos_x_sep(:,1); virtual_pos_x_sep(:,2)]; 

n_radars = 2;
n_elements_per_radar = length(virtual_pos_x)/n_radars;
radar_indices = reshape(1:length(virtual_pos_x), n_elements_per_radar, n_radars);

% === Initialize storage ===
I = cell(n_radars, n_radars);  % only monostatic for now

% ===  Visualize the array === %
% plot_mimo_array(tx_pos, rx_pos,virtual_pos);

% ===  Compute minimum ff distance for Baseline === %
D= max(virtual_pos_x) - min(virtual_pos_x); %aperture size of whole array in m
R_ff= (2*D^2)/lambda;
% fprintf('Far-field distance: %.2f meters\n', R_ff);

% === Maximum unambiguous range ===
R_max = Fs * c * T_chirp / (4 * BW);
% fprintf('Max unambiguous range: %.2f meters\n', R_max);

fprintf('Pick a target in the range: (%.2f - %.2f) meters\n', R_ff, R_max);
%% === Target Setup ===
targets = [40, deg2rad(0);   % [Range (m), Angle (rad)]
           50,deg2rad(0)];
            % 30, deg2rad(0)];

if any(targets(:,1) < R_ff)
    error('One or more targets are in the near field — computation not valid.');
end

if any(targets(:,1) > R_max)
    error('Target range exceeds maximum unambiguous range of %.2f m', R_max);
end
n_targets = size(targets,1);

%% === Signal Generation ===
t = 0:1/Fs:T_chirp-1/Fs;
n_samples = length(t);
rx_signals = zeros(length(virtual_pos_x), n_samples);

%Parameters for offsets
beta = 100 * randn(); %  CFO parameter
kappa = 1e-9 * randn();  % TO parameter
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
        phase_offsets = apply_phase_errors(virtual_pos_x(el), t, theta, beta/fc, kappa); % NOTE: normalized beta
        phase_shift_corrup = phase_shift_ideal + phase_offsets;
        signal_corrup = exp(1j * (2*pi*fb*t + phase_shift_corrup));

        % === Accumulate signals
        rx_signals_el = rx_signals_el + signal_corrup;

        % === Estimate CFO (only for one target and right-side elements)
        if el > 12 && tgt == 1
            phi_residual = unwrap(angle(signal_corrup ./ signal_ideal));
            p = polyfit(t, phi_residual, 1);
            beta_est = p(1) / (2 * pi * fc);
            beta_values(end+1) = beta_est;
        end
    end

    % Store final signal after all targets added
    rx_signals(el,:) = rx_signals_el;
end

% === Final estimate across all right-side elements
beta_est = mean(beta_values, 'omitnan')*fc;
fprintf("Estimated CFO beta: %.3e (Δf ≈ %.2f Hz)\n", beta_est, beta_est);


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

%%Localization plot
% figure;
% pcolor(angle_axis, range_axis, 20*log10(abs(RAOA.')/max(abs(RAOA(:)))));
% hold on
% for l = 1 : length(est_angles)
% plot(est_angles(l), est_ranges(l) , 'm.', 'MarkerSize', 20);
% end
% shading flat; colormap jet; caxis([-40 0]);
% xlabel('Angle (deg)'); ylabel('Range (m)');
% title('Range Angle Map with pre-sync Localization');
% hold off


%% Individual radar images
for r = 1:n_radars
    idx = radar_indices(:, r);  % indices of virtual channels for radar r
    radar_data = rx_signals(idx, :);  % select only those rows

    % === Range FFT ===
    range_fft_r = fft(radar_data, N_fft, 2);
    range_fft_r = range_fft_r(:, 1:N_fft/2);

    % === Angle FFT ===
    RAOA_r = fftshift(fft(range_fft_r, N_ffta, 1), 1);

    % === Store image ===
    I{r, r} = RAOA_r;
end

%plot all monostatic images 
% figure;
% for r = 1:n_radars
%     subplot(1,n_radars,r);
%     imagesc(angle_axis, range_axis, 20*log10(abs(I{r,r}')./max(abs(I{r,r}(:)))));
%     title(sprintf('Radar %d Monostatic Image', r));
%     xlabel('Angle (deg)'); ylabel('Range (m)');
%     caxis([-40 0]); colormap jet; shading flat;
% end

%% Step 2: Coarse synchronization

% Step 2.1
% Find best correlation between I_ref and different shifts and rotations
[rx_pos_CS, tx_pos_CS] = estimate_sensor_transformations(rx_pos,tx_pos,I);

% % Step 2.2
% % ===Reconstructed Virtual Array Positions ===
% virtual_pos_CS = cell(1, n_radars);  % use cells to hold each radar's virtual pos
% for k = 1:n_radars
%     virtual_pos_CS{k} = [];  % initialize each radar's list
% 
%     for i = 1:N_tx
%         for j = 1:N_rx
%             pos = tx_pos_CS{k}(i,:) + rx_pos_CS{k}(j,:);
%             virtual_pos_CS{k}(end+1,:) = pos;
%         end
%     end
% end
% 
% virtual_pos_all = [virtual_pos(:, 1:2); virtual_pos(:, 3:4)];
% virtual_pos_CS_all = [virtual_pos_CS{1}; virtual_pos_CS{2}];
% 
% %Questionble! Dont use target location??
% %Calculate original distance and corrected distance
% u_hat = [cos(theta), sin(theta)];         % 1×2
% r_orig = virtual_pos_all * u_hat.';           % [N × 1]
% r_corr = virtual_pos_CS_all * u_hat.';    % [N × 1]
% 
% %Calculate phase shift due to mismatch
% delta_phi = 2*pi*(r_corr - r_orig)/lambda;
% rx_signals_corrected = rx_signals .* exp(-1j * delta_phi);
% 
% % generate_localization_plots(rx_signals_corrected,SNR_dB,Fs,T_chirp,c,BW,n_targets)

%Step 2.3

for n = 1:n_radars
    for m = 1:n_radars
        tx = tx_pos_CS{n};  % TX positions for radar n
        rx = rx_pos_CS{m};  % RX positions for radar m

        virtual_pos = [];
        for i = 1:size(tx,1)
            for j = 1:size(rx,1)
                virtual_pos(end+1,:) = tx(i,:) + rx(j,:);
            end
        end

        rx_signals = zeros(size(virtual_pos,1), n_samples);

        for el = 1:size(virtual_pos,1)
            signal_el = zeros(1, n_samples);
            for tgt = 1:n_targets
                R = targets(tgt, 1);
                theta = targets(tgt, 2);
                fb = 2 * BW * R / (c * T_chirp);
                u_hat = [cos(theta), sin(theta)];
                proj = dot(virtual_pos(el,:), u_hat);
                phi = 2 * pi * proj / lambda;

                signal_el = signal_el + exp(1j * (2*pi*fb*t + phi));
            end
            rx_signals(el,:) = signal_el;
        end

        % Add noise (optional)
        signal_power = var(rx_signals(:));
        noise_power = signal_power / 10^(SNR_dB/10);
        noise = sqrt(noise_power) * randn(size(rx_signals));
        rx_signals_all{n,m} = rx_signals + noise;
    end
end

I_CS = generate_bistatic_images(rx_signals_all,N_fft,N_ffta);

% Angle axis for linear array along x
angle_axis = asind(linspace(-1, 1, N_ffta));
angle_axis = angle_axis;  % or +90 depending on the direction



figure;
imagesc(angle_axis, range_axis, 20*log10(abs(I_CS{1,2}')./max(abs(I_CS{1,2}(:)))));
xlabel('Angle (deg)');
ylabel('Range (m)');
title('Bistatic Image I_{1,2}');
colormap jet; colorbar; caxis([-40 0]);
set(gca, 'YDir', 'normal');  % range = 0 at bottom, increasing upward


figure;
imagesc(angle_axis, range_axis, 20*log10(abs(I_CS{2,1}')./max(abs(I_CS{1,2}(:)))));
xlabel('Angle (deg)');
ylabel('Range (m)');
title('Bistatic Image I_{1,2}');
colormap jet; colorbar; caxis([-40 0]);
set(gca, 'YDir', 'normal');  % range = 0 at bottom, increasing upward







% %% === Reconstructed Signal Generation ===
% t = 0:1/Fs:T_chirp-1/Fs;
% n_samples = length(t);
% rx_signals = zeros(length(virtual_pos_CS), n_samples);
% 
% for el = 1:length(virtual_pos_CS)
%     rx_signals_el = zeros(1, n_samples);  % for accumulation over targets
% 
%     for tgt = 1:n_targets
%         R = targets(tgt, 1);
%         theta = targets(tgt, 2);
% 
%         fb = 2 * BW * R / (c * T_chirp);  % Beat frequency
% 
%      % === Corrected signal (no artificial errors)
%     phase_shift_corrected = 2 * pi * virtual_pos_CS{n}(el,1) * sin(theta) / lambda;
%     signal_corrected = exp(1j * (2*pi*fb*t + phase_shift_corrected));
%     rx_signals_el = rx_signals_el + signal_corrected;
% 
%     end
% 
%     % Store final signal after all targets added
%     rx_signals(el,:) = rx_signals_el;
% 
% end
% 
