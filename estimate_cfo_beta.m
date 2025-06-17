function beta_est = estimate_cfo_beta(rx_signals, virtual_pos, t, est_range, est_angle)
% ESTIMATE_CFO_BETA - Estimate CFO (beta) from residual phase drift
%
% Output:
%   beta_est - estimated normalized CFO for the right-side unit

    fn = 78e9;
    lambda = 3e8 / fn;
    c = 3e8;
    BW = 250e6;
    T_chirp = 25.6e-6;

    fb = 2 * BW * est_range / (c * T_chirp);  % Beat frequency
    s_ideal = exp(1j * (2 * pi * fb * t));    % Ideal beat signal (no phase shift)

    is_right = virtual_pos > 0;
    phase_slopes = [];

    for el = find(is_right)'
        % Add geometric phase
        geom_phase = 2 * pi / lambda * virtual_pos(el) * sin(est_angle);
        s_ref = s_ideal .* exp(1j * geom_phase);

        % Residual phase (due to CFO)
        phase_diff = angle(rx_signals(el,:) ./ s_ref);
        phase_diff = unwrap(phase_diff);

        % Fit a line to phase_diff vs time → slope = dφ/dt
        p = polyfit(t, phase_diff, 1);
        dphi_dt = p(1);

        % Estimate beta
        beta = dphi_dt / (2 * pi * fn);
        phase_slopes(end+1) = beta;
    end

    % Average over all right-side elements
    beta_est = mean(phase_slopes);
end
