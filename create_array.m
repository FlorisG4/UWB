function [tx_pos,rx_pos,virtual_pos] = create_array(N_tx,N_rx,d,Baseline)

tx_pos = ((0:N_tx-1) - (N_tx-1)/2) * 4*d;
rx_pos = ((0:N_rx-1) - (N_rx-1)/2) * d;


y_tx = zeros(size(tx_pos));       % y-coordinates = 0 (line array)
y_rx = zeros(size(rx_pos));


% === Virtual Array Positions ===
virtual_pos = [];
for i = 1:N_tx
    for j = 1:N_rx
        virtual_pos(end+1) = tx_pos(i) + rx_pos(j);
    end
end

% Center around zero
virtual_pos = virtual_pos - mean(virtual_pos);  % ✅ automatic centering
virtual_pos = sort(virtual_pos);


%Left Array
tx_pos_l = tx_pos - Baseline/2;
rx_pos_l = rx_pos - Baseline/2;
virtual_pos_l = virtual_pos -Baseline/2; 

%Right Array
tx_pos_r = tx_pos + Baseline/2;
rx_pos_r = rx_pos + Baseline/2;
virtual_pos_r = virtual_pos + Baseline/2;

%Combine them
tx_pos = [tx_pos_l, tx_pos_r];
rx_pos = [rx_pos_l, rx_pos_r];
virtual_pos = sort([virtual_pos_l, virtual_pos_r]);