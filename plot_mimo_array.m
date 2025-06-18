function plot_mimo_array(tx_pos, rx_pos, virtual_pos)
% Simple helper to visualize antenna geometry
    figure; hold on; grid on; axis equal;

    % Plot and save handles
    hTx = scatter(tx_pos, zeros(size(tx_pos)), 100, 'g', 'filled');
    hRx = scatter(rx_pos, -0.002*ones(size(rx_pos)), 100, 'r', 'filled');
    hVirt = scatter(virtual_pos, 0.002*ones(size(virtual_pos)), 80, 'b', 'filled');

    % Labels and limits
    xlabel('x position (meters)');
    title('Tx, Rx, and Virtual Antenna Positions');
    ylim([-0.01 0.04]);
    xlim([min(virtual_pos) max(virtual_pos)]);

    % Use only the handles in the legend
    legend([hTx(1), hRx(1), hVirt(1)], {'Tx', 'Rx', 'Virtual'}, 'Location', 'northwest');
end
