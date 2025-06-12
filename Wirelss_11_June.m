clc; clear; close all;

%% === System Parameters ===
M = 64;
fc = 2.4e9;
c = 3e8;
lambda = c / fc;
d = lambda / 2;
Pt = 1;
v_human = 0.5;
max_doppler = 2 * v_human / lambda;
T_sweep = 10; fps = 30;
total_frames = T_sweep * fps;

%% === Human Setup ===
num_humans = 5;
true_human_angles = randi([0, 359], 1, num_humans);
true_human_rad = deg2rad(true_human_angles);
human_gains = 0.8 + 0.4 * rand(1, num_humans);

angles_deg = 0:359;
doppler_accum = zeros(1, length(angles_deg));

%% === Phase 1: Doppler-Based Human Detection ===
for t = 1:total_frames
    for theta = 1:length(angles_deg)
        theta_rad = deg2rad(angles_deg(theta));
        a_rx = exp(1j * 2 * pi * d * (0:M-1)' * sin(theta_rad));
        doppler_signal = 0;

        for i = 1:num_humans
            phi = true_human_rad(i);
            a_obj = exp(1j * 2 * pi * d * (0:M-1)' * sin(phi));
            gain = human_gains(i) * (a_rx' * a_obj) / norm(a_obj);
            rel_angle = cos(theta_rad - phi);
            doppler_shift = (2 * v_human * rel_angle) / lambda;
            doppler_weight = abs(doppler_shift) / max_doppler;
            doppler_signal = doppler_signal + doppler_weight * abs(gain)^2;
        end

        doppler_accum(theta) = doppler_accum(theta) + doppler_signal;
    end
end

doppler_power = doppler_accum / total_frames;
noise_floor = median(doppler_power);
threshold = noise_floor + 0.10 * (max(doppler_power) - noise_floor);
[peak_vals, locs] = findpeaks(doppler_power, 'MinPeakHeight', threshold);
detected_angles = mod(angles_deg(locs), 360);

%% === MUSIC Refinement (Optional) ===
% Assume MUSIC returns refined peaks (placeholder - insert if needed)
% For now, skip and use doppler-based peaks

%% === Merge Close Detections (Within 10°) ===
merged_angles = [];
used = false(size(detected_angles));

for i = 1:length(detected_angles)
    if used(i), continue; end
    group = detected_angles(abs(detected_angles - detected_angles(i)) <= 10);
    group = [group, detected_angles(abs(detected_angles - detected_angles(i) + 360) <= 10)];
    group = mod(unique(group), 360);
    merged_angles(end+1) = mean(group);
    used(ismember(detected_angles, group)) = true;
end
merged_angles = mod(round(merged_angles), 360);

%% === Phase 1: Log ===
fprintf('\n=== Phase 1: Human Detection Log ===\n');
fprintf('True Human Angles    : %s\n', mat2str(true_human_angles));
fprintf('Detected Human Angles: %s\n', mat2str(detected_angles));
fprintf('Merged Detections    : %s\n', mat2str(merged_angles));
fprintf('=============================\n');

%% === Phase 2: Power Allocation ===
SNR_min_db = 10;
SNR_min = 10^(SNR_min_db/10);
theta_comm = deg2rad(90);
a_comm = exp(1j * 2 * pi * d * (0:M-1)' * sin(theta_comm));
gain_90 = abs(a_comm' * a_comm)^2 / norm(a_comm)^2;
P_comm = SNR_min / gain_90;
P_comm = min(Pt, P_comm);
P_tracking = Pt - P_comm;

% Confidence = doppler peak heights
confidences = peak_vals ./ sum(peak_vals);
merged_conf = zeros(size(merged_angles));

for i = 1:length(merged_angles)
    % Merge corresponding original confidence
    close_idx = abs(detected_angles - merged_angles(i)) <= 10;
    merged_conf(i) = sum(confidences(close_idx));
end
merged_conf = merged_conf / sum(merged_conf); % Normalize weights
P_alloc = P_tracking * merged_conf;

%% === Phase 2: Log ===
fprintf('\n=== Phase 2: Power Allocation ===\n');
fprintf('Total Power Budget     : %.2f W\n', Pt);
fprintf('Power for Communication (90°): %.4f W\n', P_comm);
fprintf('Remaining Power for Tracking  : %.4f W\n\n', P_tracking);
fprintf('Tracking Beam Powers (weighted by merged confidence):\n');
for i = 1:length(merged_angles)
    fprintf('Angle %4d° → Power: %.4f W (Weight: %.2f%%)\n', ...
        merged_angles(i), P_alloc(i), 100*merged_conf(i));
end

%% === Phase 2: Plot ===
figure;
bar(categorical(string(merged_angles)), P_alloc, 'FaceColor', [0.2 0.6 0.8]);
ylabel('Power Allocated (W)');
title('Tracking Beam Power Allocation (Phase 2)');
grid on;

%% === Phase 3: Tracking Using Kalman Filter (Simple Simulation) ===
% Simulate object at detected angles with slight angular velocity
num_frames = 20;
track_history = zeros(length(merged_angles), num_frames);
Q = 1;  % Process noise
R = 2;  % Measurement noise
for i = 1:length(merged_angles)
    x = merged_angles(i);
    v = 1 + randn*0.2;
    P = 1;
    for t = 1:num_frames
        % Prediction
        x = x + v;
        P = P + Q;
        % Measurement (simulate noisy detection)
        z = x + randn()*sqrt(R);
        % Kalman update
        K = P / (P + R);
        x = x + K*(z - x);
        P = (1 - K)*P;
        track_history(i, t) = mod(x, 360);
    end
end

%% === Phase 3: Plot Tracking ===
figure;
hold on;
for i = 1:size(track_history,1)
    plot(1:num_frames, track_history(i,:), '-o');
end
xlabel('Frame');
ylabel('Angle (°)');
title('Kalman Filter-based Human Tracking (Phase 3)');
legend(arrayfun(@(x) sprintf('Obj@%d°', x), merged_angles, 'UniformOutput', false));
grid on;

figure;
pax = polaraxes;
hold(pax, 'on');
for i = 1:size(track_history, 1)
    r = ones(1, num_frames) * (i + 1); % radial offset to separate lines
    theta = deg2rad(track_history(i,:));
    polarplot(pax, theta, r, '-o', 'LineWidth', 1.5);
end
title('Polar Trajectory of Tracked Humans (Phase 3)');
rlim([0 size(track_history, 1) + 2]);
legend(arrayfun(@(x) sprintf('Obj@%d°', x), merged_angles, 'UniformOutput', false), 'Location', 'bestoutside');

% === Phase 3: Heatmap of Tracking Densities ===
angle_bins = 0:359;
heatmap_matrix = zeros(num_frames, length(angle_bins));

for t = 1:num_frames
    for i = 1:size(track_history, 1)
        angle_idx = mod(round(track_history(i,t)), 360) + 1;
        heatmap_matrix(t, angle_idx) = heatmap_matrix(t, angle_idx) + 1;
    end
end

figure;
imagesc(angle_bins, 1:num_frames, heatmap_matrix);
xlabel('Angle (°)');
ylabel('Frame');
title('Tracking Density Heatmap (Phase 3)');
colormap hot; colorbar;