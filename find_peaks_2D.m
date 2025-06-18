function [peak_coords, peak_vals] = find_peaks_2D(I, threshold_dB, neighborhood)
% Utility helper used for simple peak detection in range-angle maps
% Find local peaks in 2D image above a threshold (in dB)
% Returns coordinates and peak magnitudes
%
% I: input complex image
% threshold_dB: relative dB threshold from max (e.g. -10)
% neighborhood: size of local neighborhood (e.g. 3 or 5)

if nargin < 2
    threshold_dB = -10;
end
if nargin < 3
    neighborhood = 3;
end

mag = 20*log10(abs(I) / max(abs(I(:))));
BW = imregionalmax(mag);

% Apply threshold
BW = BW & (mag > threshold_dB);

% Find coordinates
[y, x] = find(BW);
peak_coords = [x, y];  % (col, row)
peak_vals = mag(sub2ind(size(I), y, x));
end

