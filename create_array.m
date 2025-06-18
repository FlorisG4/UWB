function [tx_pos,rx_pos,virtual_pos, virtual_pos_cells] = create_array(N_tx,N_rx,d,Baseline)
%% Creates transmit, receive and virtual arrays for a dual-radar MIMO setup

% Uniform linear arrays for TX and RX 
tx_pos = ((0:N_tx-1) - (N_tx-1)/2) * 4*d;
rx_pos = ((0:N_rx-1) - (N_rx-1)/2) * d;


y_tx = zeros(size(tx_pos));       % y-coordinates = 0 (line array)
y_rx = zeros(size(rx_pos));


% === Virtual Array Positions ===
% Build virtual array by combining every TX with every RX
virtual_pos = [];
for i = 1:N_tx
    for j = 1:N_rx
        virtual_pos(end+1) = tx_pos(i) + rx_pos(j);
    end
end

% Center around zero
virtual_pos = virtual_pos - mean(virtual_pos);  %  automatic centering
virtual_pos = sort(virtual_pos);

% Split array around origin to model two separate radars

% Left Radar (Radar 1)
tx_pos_l = tx_pos - Baseline/2;
rx_pos_l = rx_pos - Baseline/2;
virtual_pos_l = virtual_pos - Baseline/2;

% Right Radar (Radar 2)
tx_pos_r = tx_pos + Baseline/2;
rx_pos_r = rx_pos + Baseline/2;
virtual_pos_r = virtual_pos + Baseline/2;

%Add y coordinates (Which are all zeros)
tx_pos_l = [tx_pos_l(:), zeros(size(tx_pos_l(:)))];
tx_pos_r = [tx_pos_r(:), zeros(size(tx_pos_r(:)))];
rx_pos_l = [rx_pos_l(:), zeros(size(rx_pos_l(:)))];
rx_pos_r = [rx_pos_r(:), zeros(size(rx_pos_r(:)))];
virtual_pos_l = [virtual_pos_l(:), zeros(size(virtual_pos_l(:)))];
virtual_pos_r = [virtual_pos_r(:), zeros(size(virtual_pos_r(:)))];

% === Store per radar in cell arrays ===
tx_pos = cell(2,1);
rx_pos = cell(2,1);
virtual_pos_cells = cell(2,1);

tx_pos{1} = tx_pos_l;
rx_pos{1} = rx_pos_l;
virtual_pos_cells{1} = virtual_pos_l;

tx_pos{2} = tx_pos_r;
rx_pos{2} = rx_pos_r;
virtual_pos_cells{2} = virtual_pos_r;

% Optional combined arrays (if you still want them for global processing)
virtual_pos = sort([virtual_pos_l, virtual_pos_r]);


