
% =========================================================================
% UWB MIMO Radar Localization Simulation
% 
% =========================================================================

% Main script to simulate a dual-radar MIMO localization scenario
clear; clc;

%% === Simulation Parameters ===
P = uwb_params();
c = P.c;
fc = P.fc;
lambda = P.lambda;
Fs = P.Fs;
T_chirp = P.T_chirp;
BW = P.BW;
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

% % ===  Visualize the array === %
% plot_mimo_array(sort([tx_pos{1}(:,1), tx_pos{2}(:,1)]), sort([rx_pos{1}(:,1),rx_pos{2}(:,1)]),virtual_pos_x);

% ===  Compute minimum ff distance for Baseline === %
D= max(virtual_pos_x) - min(virtual_pos_x); %aperture size of whole array in m
R_ff= (2*D^2)/lambda;
% fprintf('Far-field distance: %.2f meters\n', R_ff);

% === Maximum unambiguous range ===
R_max = Fs * c * T_chirp / (4 * BW);
% fprintf('Max unambiguous range: %.2f meters\n', R_max);

fprintf('Pick a target in the range: (%.2f - %.2f) meters\n', R_ff, R_max);
%% === Target Setup ===
targets = [40, deg2rad(10);   % [Range (m), Angle (rad)]
            50,deg2rad(-10);
             30, deg2rad(0)];

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
        phase_offsets = apply_phase_errors(virtual_pos_x(el), t, theta, beta_true/fc, kappa_true); % NOTE: normalized beta
        % phase_offsets = 0;
        phase_shift_corrup = phase_shift_ideal + phase_offsets;
        signal_corrup = exp(1j * (2*pi*fb*t + phase_shift_corrup));

        % === Accumulate signals
        rx_signals_el = rx_signals_el + signal_corrup;

        % === Estimate CFO (only for one target and right-side elements)
        % Derive CFO by measuring residual phase slope
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
fprintf("Estimated CFO beta: %.3e (Δf ≈ %.2f Hz)\n", beta_est_norm, beta_est);

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
% [est_ranges, est_angles, peak_values] = localize_targets_peak(RAOA, range_axis, angle_axis, n_targets);
% 
% evaluation(targets,est_ranges,est_angles)
% 
% 
% %Localization plot
% figure;
% pcolor(angle_axis, range_axis, 20*log10(abs(RAOA.')/max(abs(RAOA(:)))));
% hold on
% for l = 1 : length(est_angles)
% plot(est_angles(l), est_ranges(l) , 'm.', 'MarkerSize', 18);
% end
% shading flat; colormap jet; caxis([-40 0]);
% xlabel('Angle (deg)'); ylabel('Range (m)');
% title('Pre-sync Localization with offsets');
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
[rx_pos_CS, tx_pos_CS, virtual_pos_CS] = estimate_sensor_transformations(rx_pos,tx_pos,virtual_pos,I);


%Step 2.2 and 2.3
% % ===Reconstructed Virtual Array Positions ===
for n = 1:n_radars
    for m = 1:n_radars
        tx = tx_pos_CS{n};
        rx = rx_pos_CS{m};

        virtual_pos = [];
        for i = 1:size(tx,1)
            for j = 1:size(rx,1)
                virtual_pos(end+1,:) = tx(i,:) + rx(j,:);
            end
        end

        rx_block = zeros(size(virtual_pos,1), n_samples);

        for el = 1:size(virtual_pos,1)
            signal_el_offsets = zeros(1, n_samples);

            for tgt = 1:n_targets
                R = targets(tgt, 1);
                theta = targets(tgt, 2);
                fb = 2 * BW * R / (c * T_chirp);
                u_hat = [sin(theta), cos(theta)];
                proj = dot(virtual_pos(el,:), u_hat);
                phi = 2 * pi * proj / lambda;
                phase_offsets = apply_phase_errors(virtual_pos(el,1), t, theta, beta_true/fc, kappa_true);

                signal_corr  = exp(1j * (2*pi*fb*t + phi + phase_offsets));
                signal_el_offsets = signal_el_offsets + signal_corr;
            end

            rx_block(el,:) = signal_el_offsets;
        end

        % Add noise
        signal_power = var(rx_block(:));
        noise_power = signal_power / 10^(SNR_dB/10);
        noise = sqrt(noise_power) * randn(size(rx_block));

        rx_signals_all{n,m} = rx_block + noise;
    end
end

I_CS = generate_bistatic_images(rx_signals_all, N_fft, N_ffta);

% Pre Coarse sync fusion
I_fused = fuse_bistatic_images(I_CS, I, 'coherent');

% figure;
% imagesc(angle_axis, range_axis, 20*log10(abs(I_fused')./max(abs(I_fused(:)))));
% xlabel('Angle (deg)');
% ylabel('Range (m)');
% title('Fused Image - Presync');
% colormap jet; colorbar; caxis([-40 0]);
% set(gca, 'YDir', 'normal');

[est_ranges, est_angles, peak_values] = localize_targets_peak(I_fused, range_axis, angle_axis, n_targets);
evaluation(targets,est_ranges,est_angles)
%Localization plot
figure;
pcolor(angle_axis, range_axis, 20*log10(abs(I_fused.')/max(abs(I_fused(:)))));
hold on
for l = 1 : length(est_angles)
plot(est_angles(l), est_ranges(l) , 'm.', 'MarkerSize', 18);
end
shading flat; colormap jet; caxis([-40 0]);
xlabel('Angle (deg)'); ylabel('Range (m)');
title('Pre-Coarse sync fusion plot');
hold off

% figure;
% imagesc(angle_axis, range_axis, 20*log10(abs(I_CS{1,2}')./max(abs(I_CS{1,2}(:)))));
% xlabel('Angle (deg)');
% ylabel('Range (m)');
% title('Bistatic Image I_{1,2}');
% colormap jet; colorbar; caxis([-40 0]);
% set(gca, 'YDir', 'normal');  % range = 0 at bottom, increasing upward
% 
% figure;
% imagesc(angle_axis, range_axis, 20*log10(abs(I_CS{2,1}')./max(abs(I_CS{2,1}(:)))));
% xlabel('Angle (deg)');
% ylabel('Range (m)');
% title('Bistatic Image I_{2,1}');
% colormap jet; colorbar; caxis([-40 0]);
% set(gca, 'YDir', 'normal');  


%Kappa estimation
I_ref = I{1,1};
I_nm = I_CS{1,2};
I_mn = I_CS{2,1};

N_size = 5; % Local region around each pixel to determine if its a local max
dB_thr = -1; % Threshold relative to peak

[kappa_est_nm, tof_nm] = estimate_kappa(I_nm, I_ref, tx_pos_CS, beta_est_norm,dB_thr, N_size); 
[kappa_est_mn, tof_mn] = estimate_kappa(I_mn, I_ref, tx_pos_CS,beta_est_norm, dB_thr, N_size);

fprintf("Actual TO kappa is: %.3e \n", kappa_true);

fprintf("Estimated TO kappa_nm: %.3e \n", kappa_est_nm );
fprintf("Estimated TO kappa_mn: %.3e \n", kappa_est_mn );



% Subtracting the errors
t_corr_nm = (1)*tof_nm + kappa_est_nm;
t_corr_mn = (1+beta_est_norm)*tof_mn + kappa_est_mn;


phase_corr_nm = exp(+j*2*pi*fc*(1)*tof_nm) ;
phase_corr_mn = exp(+j*2*pi*fc*(1+beta_est_norm)*tof_mn) ;

% Physical units in meters
x_min = -90; x_max = 90;     % lateral axis
y_min =  0; y_max = 76.80;     % range axis (downfield)
N_x = 1024;  % number of pixels in x
N_y = 512;  % number of pixels in y
x_vals = linspace(x_min, x_max, N_x);
y_vals = linspace(y_min, y_max, N_y);
[image_grid_x, image_grid_y] = meshgrid(x_vals, y_vals);  % size: [N_y × N_x]


[Ny, Nx] = size(image_grid_x);  % number of pixels
I_nm = zeros(Nx, Ny);           % final image
I_mn = zeros(Nx, Ny);


% Caluclating n_m
kappa = kappa_est_nm;
tx_pos = tx_pos_CS{1,1};          % TX positions (3×2)
rx_pos = rx_pos_CS{1,2};          % RX positions (4×2)
virtual_pos = virtual_pos_CS{1,2};  % (12×2)


for yi = 1:Ny  % range (rows)
    for xi = 1:Nx  % angle (columns)
        x = [image_grid_x(yi, xi), image_grid_y(yi, xi)];  % [x y]

        % === Compute TOF for each virtual element ===
        for el = 1:size(virtual_pos,1)
            tx = tx_pos(mod(el-1, N_tx) + 1, :);
            rx = rx_pos(floor((el-1)/N_tx) + 1, :);

            tof = (norm(x - tx) + norm(x - rx)) / c;
            t_corr = (1) * tof + kappa;

            % Interpolate signal at corrected time
            s = interp1(t, rx_signals_all{1,2}(el,:), t_corr, 'linear', 0);

            % Phase correction
            phase = exp(1j * 2 * pi * fc * (1) * tof);

            % Accumulate contribution
            I_nm(xi, yi) = I_nm(xi, yi) + s * phase;
        end
    end
end
I_CS{1,2} = I_nm;

%Caluclating m_n
beta = beta_est_norm;
kappa = kappa_est_mn;
tx_pos = tx_pos_CS{1,2};          % TX positions (3×2)
rx_pos = rx_pos_CS{1,1};          % RX positions (4×2)
virtual_pos = virtual_pos_CS{2,1};  % (12×2)


for yi = 1:Ny  % range (rows)
    for xi = 1:Nx  % angle (columns)
        x = [image_grid_x(yi, xi), image_grid_y(yi, xi)];  % [x y]

        %=== Compute TOF for each virtual element ===
        for el = 1:size(virtual_pos,1)
            tx = tx_pos(mod(el-1, N_tx) + 1, :);
            rx = rx_pos(floor((el-1)/N_tx) + 1, :);

            tof = (norm(x - tx) + norm(x - rx)) / c;
            t_corr = (1 + beta) * tof + kappa;

            %Interpolate signal at corrected time
            s = interp1(t, rx_signals_all{2,1}(el,:), t_corr, 'linear', 0);

            %Phase correction
            phase = exp(1j * 2 * pi * fc * (1 + beta) * tof);

            %Accumulate contribution
            I_mn(xi, yi) = I_mn(xi, yi) + s * phase;
        end
    end
end
I_CS{2,1} = I_mn;

% Beta substraction from {2,2}
    I_corr = I_CS{2,2};
    for xi = 1:Nx
        for yi = 1:Ny
            x_pix = [image_grid_x(yi, xi), image_grid_y(yi, xi)];
            tx = tx_pos_CS{2}(ceil(end/2), :);  % assume central TX (or average)
            tau = 2 * norm(x_pix - tx) / c;
            phase_corr = exp(-1j * 2*pi * fc * beta * tau);
            I_corr(xi, yi) = I_corr(xi, yi) * phase_corr;
        end
    end
    I_CS{2,2} = I_corr;



% % === Step 2 Last Bistatic image fusion ===
I_fused = fuse_bistatic_images(I_CS, I, 'coherent');


[est_ranges, est_angles, peak_values] = localize_targets_peak(I_fused, range_axis, angle_axis, n_targets);
evaluation(targets,est_ranges,est_angles)
%Localization plot
figure;
pcolor(angle_axis, range_axis, 20*log10(abs(I_fused.')/max(abs(I_fused(:)))));
hold on
for l = 1 : length(est_angles)
plot(est_angles(l), est_ranges(l) , 'm.', 'MarkerSize', 18);
end
shading flat; colormap jet; caxis([-40 0]);
xlabel('Angle (deg)'); ylabel('Range (m)');
title('Post-Coarse sync fusion plot');
hold off

% figure;
% imagesc(angle_axis, range_axis, 20*log10(abs(I_fused')./max(abs(I_fused(:)))));
% xlabel('Angle (deg)');
% ylabel('Range (m)');
% title('Fused Image');
% colormap jet; colorbar; caxis([-40 0]);
% set(gca, 'YDir', 'normal');
% 
% 
% %% Step 3 
% % === Section III.C/D: Residual Phase Error Estimation and Correction ===
% 
% % Stap 1: Detecteer peaks in beide fine-synced beelden
% dB_threshold = -10;
% neigh_size = 5;
% [peaks_12, ~] = find_peaks_2D(abs(I_CS{1,2}), dB_threshold, neigh_size);
% [peaks_21, ~] = find_peaks_2D(abs(I_CS{2,1}), dB_threshold, neigh_size);
% 
% % Stap 2: Vind de meest overlappende peak
% min_dist = Inf;
% for i = 1:size(peaks_12,1)
%     for j = 1:size(peaks_21,1)
%         dist = norm(peaks_12(i,:) - peaks_21(j,:));
%         if dist < min_dist
%             min_dist = dist;
%             idx_12 = i;
%             idx_21 = j;
%         end
%     end
% end
% 
% % Stap 3: Extract complex values at peak locations
% val_12 = I_CS{1,2}(peaks_12(idx_12,2), peaks_12(idx_12,1));
% val_21 = I_CS{2,1}(peaks_21(idx_21,2), peaks_21(idx_21,1));
% 
% 
% % Stap 4: Bereken residuele faseverschil
% delta_phi = angle(val_12) - angle(val_21);  % moet in [-pi, pi]
% 
% % Stap 5: Bereken kappa (en eventueel delta_beta)
% % We nemen aan: τₙ en τₘ ≈ tof_nm en tof_mn
% delta_tof = tof_nm - tof_mn;
% 
% % === Gecombineerde foutmodel:
% % Δϕ ≈ 2π f_c (β_est * Δτ + κ_n - κ_m)
% % ⇒ κ_n - κ_m ≈ (Δϕ - 2π f_c β_est * Δτ) / (2π f_c)
% kappa_diff_est = (delta_phi - 2*pi*fc*beta_est*delta_tof) / (2*pi*fc);
% 
% fprintf('Residual κ_n - κ_m estimate: %.3e seconds\n', kappa_diff_est);
% 
% % ⇒ Je kunt deze waarde gebruiken om je eerdere kappa-schattingen te verfijnen:
% kappa_est_nm = kappa_est_nm + kappa_diff_est/2;
% kappa_est_mn = kappa_est_mn - kappa_diff_est/2;




%% calculating fusion with actual kappa and beta
% Caluclating n_m
kappa = kappa_true;
tx_pos = tx_pos_CS{1,1};          % TX positions (3×2)
rx_pos = rx_pos_CS{1,2};          % RX positions (4×2)
virtual_pos = virtual_pos_CS{1,2};  % (12×2)


for yi = 1:Ny  % range (rows)
    for xi = 1:Nx  % angle (columns)
        x = [image_grid_x(yi, xi), image_grid_y(yi, xi)];  % [x y]

        % === Compute TOF for each virtual element ===
        for el = 1:size(virtual_pos,1)
            tx = tx_pos(mod(el-1, N_tx) + 1, :);
            rx = rx_pos(floor((el-1)/N_tx) + 1, :);

            tof = (norm(x - tx) + norm(x - rx)) / c;
            t_corr = (1) * tof + kappa;

            % Interpolate signal at corrected time
            s = interp1(t, rx_signals_all{1,2}(el,:), t_corr, 'linear', 0);

            % Phase correction
            phase = exp(1j * 2 * pi * fc * (1) * tof);

            % Accumulate contribution
            I_nm(xi, yi) = I_nm(xi, yi) + s * phase;
        end
    end
end
I_CS{1,2} = I_nm;

%Caluclating m_n
beta = beta_est_norm;
kappa = kappa_true;
tx_pos = tx_pos_CS{1,2};          % TX positions (3×2)
rx_pos = rx_pos_CS{1,1};          % RX positions (4×2)
virtual_pos = virtual_pos_CS{2,1};  % (12×2)


for yi = 1:Ny  % range (rows)
    for xi = 1:Nx  % angle (columns)
        x = [image_grid_x(yi, xi), image_grid_y(yi, xi)];  % [x y]

        %=== Compute TOF for each virtual element ===
        for el = 1:size(virtual_pos,1)
            tx = tx_pos(mod(el-1, N_tx) + 1, :);
            rx = rx_pos(floor((el-1)/N_tx) + 1, :);

            tof = (norm(x - tx) + norm(x - rx)) / c;
            t_corr = (1 + beta) * tof + kappa;

            %Interpolate signal at corrected time
            s = interp1(t, rx_signals_all{2,1}(el,:), t_corr, 'linear', 0);

            %Phase correction
            phase = exp(1j * 2 * pi * fc * (1 + beta) * tof);

            %Accumulate contribution
            I_mn(xi, yi) = I_mn(xi, yi) + s * phase;
        end
    end
end
I_CS{2,1} = I_mn;

I_fused = fuse_bistatic_images(I_CS, I, 'coherent');

[est_ranges, est_angles, peak_values] = localize_targets_peak(I_fused, range_axis, angle_axis, n_targets);
evaluation(targets,est_ranges,est_angles)
%Localization plot
figure;
pcolor(angle_axis, range_axis, 20*log10(abs(I_fused.')/max(abs(I_fused(:)))));
hold on
for l = 1 : length(est_angles)
plot(est_angles(l), est_ranges(l) , 'm.', 'MarkerSize', 18);
end
shading flat; colormap jet; caxis([-40 0]);
xlabel('Angle (deg)'); ylabel('Range (m)');
title('Ideal sync fusion plot');
hold off

