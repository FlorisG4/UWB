function [est_ranges, est_angles, peak_values] = localize_targets_peak(RAOA, range_axis, angle_axis, n_targets, displayResults)
% LOCALIZE_TARGETS_PEAK - Estimate [range, angle] of multiple targets from RAOA map
% % Uses simple peak picking on the RA map to estimate target locations
% Inputs:
%   RAOA_mag       - Magnitude of the range-angle map [angle_bins x range_bins]
%   range_axis     - Vector of range values corresponding to rows of RAOA
%   angle_axis     - Vector of angle values corresponding to columns of RAOA
%   n_targets      - Number of targets to detect
%   displayResults - (optional) if true, prints estimated results to console
%
% Outputs:
%   est_ranges     - Estimated range values [n_targets x 1]
%   est_angles     - Estimated angle values [n_targets x 1]
%   peak_values    - Peak magnitude at each estimated location

    if nargin < 5
        displayResults = true;  % default behavior
    end
    RAOA_mag = abs(RAOA);
    % Find strongest peaks
    [flat_peaks, flat_locs] = findpeaks(RAOA_mag(:), ...
        'SortStr', 'descend', 'NPeaks', n_targets);

    [angle_idx, range_idx] = ind2sub(size(RAOA_mag), flat_locs);

    est_ranges = range_axis(range_idx);
    est_angles = angle_axis(angle_idx);
    peak_values = flat_peaks(:);

    % Optional printing
    if displayResults
        fprintf('\n=== Estimated Target Locations ===\n');
        fprintf('  %-8s %-10s %-10s %-10s\n', 'Target', 'Range (m)', 'Angle (deg)', 'Peak Value');
        fprintf('  ---------------------------------------------\n');

        for i = 1:n_targets
            fprintf('  %-8d %-10.2f %-10.2f %-10.2f\n', ...
                i, est_ranges(i), est_angles(i), peak_values(i));
        end

        fprintf('  ---------------------------------------------\n\n');
    end
end
