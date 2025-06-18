function [kappa_est, tof_nm] = estimate_kappa(I_nm,I_ref, tx_pos_CS, beta_est, dB_thr, N_size)

P = uwb_params();
c = P.c;
fc = P.fc;
Fs = P.Fs;
T_chirp = P.T_chirp;
BW = P.BW;

N_fft =1024;
N_ffta = 1024;
range_axis = ((0:N_fft/2-1) * Fs / N_fft) * (c * T_chirp / (2 * BW));
angle_axis = asind(linspace(-1,1,N_ffta));



[peaks_nm, vals_nm] = find_peaks_2D(I_nm, dB_thr, N_size);
[peaks_ref, vals_ref] = find_peaks_2D(I_ref, dB_thr, N_size);

% Match closest pair (Euclidean distance)
min_dist = Inf;
for i = 1:size(peaks_nm,1)
    for j = 1:size(peaks_ref,1)
        dist = norm(peaks_nm(i,:) - peaks_ref(j,:));
        if dist < min_dist
            min_dist = dist;
            idx_nm = i;
            idx_ref = j;
        end
    end
end

coord_nm  = peaks_nm(idx_nm, :);    % [x, y] in I_nm
coord_ref = peaks_ref(idx_ref, :); % [x, y] in I_ref

range_nm  = range_axis(coord_nm(:,1));
angle_nm  = angle_axis(coord_nm(:,2));

range_ref = range_axis(coord_ref(:,1));
angle_ref = angle_axis(coord_ref(:,2));

u_nm  = [cosd(angle_nm), sind(angle_nm)];
u_ref = [cosd(angle_ref), sind(angle_ref)];

x_nm  = range_nm  * u_nm;
x_ref = range_ref * u_ref;

s_n = tx_pos_CS{1,1}(2,:);
s_m = tx_pos_CS{1,2}(2,:);

tof_nm  = (norm(x_nm - s_n) + norm(s_m - x_nm)) / c;
tof_ref = (norm(x_ref - s_n) + norm(s_m - x_ref)) / c;

alpha =1e3; % kappa is still a lot of nonsense fix later

kappa_est =  (1 / (1 + beta_est)) * (tof_nm - tof_ref);