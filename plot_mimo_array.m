function plot_mimo_array(tx_pos, rx_pos,virtual_pos)

% === Plot ===
figure; hold on; grid on; axis equal;

% Plot Tx in green
scatter(tx_pos, zeros(size(tx_pos)), 100, 'g', 'filled', 'DisplayName', 'Tx');

% Plot Rx in red
scatter(rx_pos, -0.002*ones(size(rx_pos)), 100, 'r', 'filled', 'DisplayName', 'Rx');

% Plot Virtual in blue
scatter(virtual_pos, 0.002*ones(size(virtual_pos)), 80, 'b', 'filled', 'DisplayName', 'Virtual');


% Label axes and legend
xlabel('x position (meters)');
% ylabel('y offset (for clarity)');
title('Tx, Rx, and Virtual Antenna Positions');
ylim([-0.01 0.04]);
xlim([min(virtual_pos) max(virtual_pos)]);
legend('Location', 'northwest');