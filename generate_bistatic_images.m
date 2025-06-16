function I_CS = generate_bistatic_images(rx_signals_all, N_fft, N_ffta)
% Generate bistatic images from actual received signals
% Inputs:
%   rx_signals_all  - cell{n,m} of [N_virtual × time] beat signals
%   virtual_pos_CS  - cell{n,m} of [x,y] positions (optional, for consistency)
%   N_fft, N_ffta   - FFT sizes for range and angle
% Output:
%   I_CS{n,m}       - FFT image for radar pair n,m

N = size(rx_signals_all, 1);
I_CS = cell(N, N);

for n = 1:N
    for m = 1:N
        rx_signals = rx_signals_all{n,m};  % [N_virtual × N_samples]

        % === Range FFT
        range_fft = fft(rx_signals, N_fft, 2);
        range_fft = range_fft(:, 1:N_fft/2);
        range_fft = flipud(range_fft);

        % === Angle FFT
        RAOA = fftshift(fft(range_fft, N_ffta, 1), 1);

        % === Store image
        I_CS{n,m} = RAOA;
    end
end
end
