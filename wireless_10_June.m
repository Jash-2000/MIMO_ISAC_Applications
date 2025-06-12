clc; clear;

% === System Parameters ===
M = 64;
fc = 2.4e9;
c = 3e8;
lambda = c / fc;
d = lambda / 2;
Pt = 1;

% === Doppler and Time Parameters ===
v_human = 1;  % m/s walking speed
max_doppler = 2 * v_human / lambda;  % Max Doppler shift
T_sweep = 5;   % Sweep duration in seconds
fps = 30;      % Sweep resolution: 30 sweeps per second
total_frames = T_sweep * fps;

% === Beam Setup ===
N = 8;  % Number of beams
beam_dirs = linspace(0, 360, N+1); beam_dirs(end) = [];
beam_dirs_rad = deg2rad(beam_dirs);

% === Human Setup ===
num_humans = 3;
human_angles = randi([0, 359], 1, num_humans);
human_angles_rad = deg2rad(human_angles);
human_gains = 0.8 + 0.4 * rand(1, num_humans);  % Higher, consistent gain

% === Reflection Sweeping ===
angles_deg = 0:359;
received_power = zeros(N, length(angles_deg));
doppler_accum = zeros(1, length(angles_deg));

% === Simulate Temporal Sweeping ===
for t = 1:total_frames
    for b = 1:N
        for theta = 1:length(angles_deg)
            theta_rad = deg2rad(angles_deg(theta));
            a_rx = exp(1j * 2 * pi * d * (0:M-1)' * sin(theta_rad));
            
            reflection = 0;
            doppler_signal = 0;
            
            for i = 1:num_humans
                phi = human_angles_rad(i);
                a_obj = exp(1j * 2 * pi * d * (0:M-1)' * sin(phi));
                gain = human_gains(i) * (a_rx' * a_obj) / norm(a_obj);
                
                rel_angle = cos(theta_rad - phi);
                doppler_shift = (2 * v_human * rel_angle) / lambda;
                doppler_weight = abs(doppler_shift) / max_doppler;

                reflection = reflection + abs(gain)^2;
                doppler_signal = doppler_signal + doppler_weight * abs(gain)^2;
            end
            
            received_power(b, theta) = reflection;
            doppler_accum(theta) = doppler_accum(theta) + doppler_signal;
        end
    end
end

% === Normalize Doppler Power ===
doppler_power = doppler_accum / total_frames;

% === Doppler-based Detection ===
noise_floor = median(doppler_power);  % Adaptive noise baseline
threshold = noise_floor + 0.06* (max(doppler_power) - noise_floor);
[~, locs] = findpeaks(doppler_power, 'MinPeakHeight', threshold);
detected_angles = mod(angles_deg(locs), 360);

% === Output ===
fprintf('\n=== Human Detection Log ===\n');
fprintf('True Human Angles    : %s\n', mat2str(human_angles));
fprintf('Detected Human Angles: %s\n', mat2str(detected_angles));
fprintf('=============================\n');

% === Plot Results ===
figure;
polarplot(deg2rad(angles_deg), doppler_power, 'b', 'LineWidth', 1.4); hold on;
polarplot(deg2rad(human_angles), max(doppler_power) * ones(1, num_humans), 'ro', 'MarkerSize', 8, 'DisplayName', 'True Human');
polarplot(deg2rad(detected_angles), max(doppler_power) * ones(1, length(detected_angles)), 'gx', 'MarkerSize', 8, 'DisplayName', 'Detected');
title('Enhanced Doppler-Based Human Detection (5s Sweep)');
legend('Doppler Power', 'True Human', 'Detected Human');