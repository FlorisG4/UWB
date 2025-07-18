function beta_est = estimate_cfo_beta(rx_signals, virtual_pos, t, est_range, est_angle)
% Estimates clock frequency offset based on phase drift across right-side
% array elements
% ESTIMATE_CFO_BETA - Estimate CFO (beta) from residual phase drift
%
% Output:
%   beta_est - estimated normalized CFO for the right-side unit

    P = uwb_params();
    fn = P.fc;
    lambda = P.lambda;
    c = P.c;
    BW = P.BW;
    T_chirp = P.T_chirp;    

    % Compute ideal beat signal for assumed target range
    fb = 2 * BW * est_range / (c * T_chirp);  % Beat frequency
    s_ideal = exp(1j * (2 * pi * fb * t));    % Ideal beat signal (no phase shift)
    % Consider only elements on the right side of the array
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
        % Linear fit to estimate phase drift over time
        p = polyfit(t, phase_diff, 1);
        dphi_dt = p(1);

        % Estimate beta
        beta = dphi_dt / (2 * pi * fn);
        phase_slopes(end+1) = beta;
    end

    % Average over all right-side elements
    beta_est = mean(phase_slopes);
end
