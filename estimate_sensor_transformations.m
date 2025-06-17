function [rx_pos_CS, tx_pos_CS, virtual_pos_CS] = estimate_sensor_transformations(rx_pos, tx_pos, virtual_pos, I)
% Coarse synchronization: align monostatic images via rigid-body transform
% Input: I - cell array where I{n,n} is monostatic image of radar n
% Output:
%   rotations: estimated rotation angles [rad]
%   translations: estimated shifts [dx, dy] in pixels

N = size(I, 1);           % Number of radars
I_ref = abs(I{1,1});      % Reference image
rotations = zeros(N,1);   
translations = zeros(N,2);

% --- Parameters for search ---
theta_range = deg2rad(-5:1:5);     % Search over small rotation angles
shift_range = -5:5;                % Search over pixel shifts (can refine later)
[dx_grid, dy_grid] = meshgrid(shift_range, shift_range);

for n = 2:N
    I_n = abs(I{n,n});
    best_corr = -Inf;
    best_theta = 0;
    best_shift = [0,0];
    
    for theta = theta_range
        % Rotate image
        I_rot = imrotate(I_n, rad2deg(theta), 'crop');  % use 'crop' to keep size
        
        for idx = 1:numel(dx_grid)
            dx = dx_grid(idx);
            dy = dy_grid(idx);
            
            % Shift image
            I_shifted = circshift(I_rot, [dy, dx]);
            
            % Compute similarity
            corr_val = sum(sum(I_ref .* I_shifted));
            
            if corr_val > best_corr
                best_corr = corr_val;
                best_theta = theta;
                best_shift = [dx, dy];
            end
        end
    end
    
    rotations(n) = best_theta;
    translations(n,:) = best_shift;
end

% radar 1 is reference
rotations(1) = 0;
translations(1,:) = [0 0];
ref_center = 0;

tx_pos_CS{1} = tx_pos{1};
rx_pos_CS{1} = rx_pos{1};
virtual_pos_CS{1} = virtual_pos(:,1:2);

for n = 2:N
    R = [cos(rotations(n)), -sin(rotations(n));
         sin(rotations(n)),  cos(rotations(n))];

    % --- TX elements ---
    for k = 1:size(tx_pos{n},1)
        tx_pos_CS{n}(k,:) = (R * (tx_pos{n}(k,:) - ref_center)')' + translations(n,:);
    end

    % --- RX elements ---
    for k = 1:size(rx_pos{n},1)
        rx_pos_CS{n}(k,:) = (R * (rx_pos{n}(k,:) - ref_center)')' + translations(n,:);
    end
end
   % Initialize full 2×2 virtual array cell
virtual_pos_CS = cell(N, N);  % N = number of radars

for n = 1:N          % TX radar
    for m = 1:N      % RX radar
        virtual_pos = [];
        for i = 1:size(tx_pos_CS{n}, 1)
            for j = 1:size(rx_pos_CS{m}, 1)
                virtual_pos(end+1, :) = tx_pos_CS{n}(i,:) + rx_pos_CS{m}(j,:);
            end
        end
        virtual_pos_CS{n, m} = virtual_pos;
    end
end

