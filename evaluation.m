function evaluation(targets, est_ranges, est_angles)


true_positions = sort(targets);
est_positions = sort([est_ranges' , deg2rad(est_angles')]);
errors = est_positions - true_positions;
rmse = sqrt(mean(sum(errors.^2, 2)));

%FAR
tol = 0.5; % distance threshold for a valid detection
false_alarms = 0;

for i = 1:size(est_positions, 1)
    dists = vecnorm(true_positions - est_positions(i,:), 2, 2);
    if min(dists) > tol
        false_alarms = false_alarms + 1;
    end
end

FAR = false_alarms / size(est_positions, 1);

%Detection probability
detections = 0;

for i = 1:size(true_positions, 1)
    dists = vecnorm(est_positions - true_positions(i,:), 2, 2);
    if min(dists) < tol
        detections = detections + 1;
    end
end

P_det = detections / size(true_positions, 1);

% === Evaluation Summary ===
fprintf('========================================\n');
fprintf('        Radar Evaluation Metrics         \n');
fprintf('========================================\n');
fprintf('Root Mean Square Error (RMSE):    %.3f meters\n', rmse);
fprintf('False Alarm Rate (FAR):           %.2f %%\n', FAR * 100);
fprintf('Detection Probability (P_det):    %.2f %%\n', P_det * 100);
fprintf('========================================\n');
