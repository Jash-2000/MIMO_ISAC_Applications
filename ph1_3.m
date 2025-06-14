clc; clear; close all;

%% === SYSTEM PARAMETERS ===
fc = 2.4e9; c = 3e8;
lambda = c / fc; d = lambda / 2;
Pt = 1; M = 64;
v_human = 0.5; max_doppler = 2 * v_human / lambda;
T_sweep = 10; fps = 30; total_frames = T_sweep * fps;
angles_deg = 0:359;

%% === HUMAN SETUP ===
num_humans = 5;
true_human_angles = randi([0, 359], 1, num_humans);
true_human_rad = deg2rad(true_human_angles);
human_gains = 0.8 + 0.4 * rand(1, num_humans);
doppler_accum = zeros(1, length(angles_deg));

%% === PHASE 1: DOPPLER DETECTION ===
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

% Merge nearby detections
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

fprintf('\n=== PHASE 1: DETECTION LOG ===\n');
fprintf('True Human Angles    : %s\n', mat2str(true_human_angles));
fprintf('Detected Human Angles: %s\n', mat2str(detected_angles));
fprintf('Merged Detections    : %s\n', mat2str(merged_angles));

%% === PHASE 2: INITIAL POWER ALLOCATION ===
SNR_min = 10^(10/10);
theta_comm = deg2rad(90);
a_comm = exp(1j * 2 * pi * d * (0:M-1)' * sin(theta_comm));
gain_90 = abs(a_comm' * a_comm)^2 / norm(a_comm)^2;
P_comm_static = min(Pt, SNR_min / gain_90);
P_tracking_init = Pt - P_comm_static;

confidences = peak_vals ./ sum(peak_vals);
merged_conf = zeros(size(merged_angles));
for i = 1:length(merged_angles)
    close_idx = abs(detected_angles - merged_angles(i)) <= 10;
    merged_conf(i) = sum(confidences(close_idx));
end
merged_conf = merged_conf / sum(merged_conf);
P_alloc = P_tracking_init * merged_conf;

fprintf('\n=== PHASE 2: POWER ALLOCATION ===\n');
fprintf('Static Communication Power (90°): %.4f W\n', P_comm_static);
fprintf('Initial Tracking Power : %.4f W\n', P_tracking_init);
for i = 1:length(merged_angles)
    fprintf('Angle %3d° → Power: %.4f W (%.2f%%)\n', ...
        merged_angles(i), P_alloc(i), 100 * merged_conf(i));
end

figure;
bar(categorical(string(merged_angles)), P_alloc, 'FaceColor', [0.2 0.6 0.8]);
ylabel('Power (W)'); grid on;
title('Phase 2: Initial Power Allocation (Tracking)');

%% === PHASE 3: ADAPTIVE TRACKING + SNR-AWARE POWER ===
beamAz = merged_angles;
B = length(beamAz); T = 50;
Nt = 8; sigma2n = 1e-3;
steer = @(th) exp(1j*2*pi*d/lambda*(0:Nt-1).' * sind(th)) / sqrt(Nt);
objAng = repelem(beamAz', 1, T);
K = size(objAng,1); RCS = ones(1,K);

% Smooth SNR variation
t_vec = linspace(0, 2*pi, T);
SNR_dB_dynamic = 12.5 + 2.5 * sin(t_vec) + 1.0 * sin(2.5 * t_vec + pi/3);
SNR_dB_dynamic = min(max(SNR_dB_dynamic, 10), 15);

% Compute baseline comm power for 10 dB
P_comm_baseline = min(Pt, 10^(10/10) / gain_90);
fprintf('\n=== BASELINE COMMUNICATION POWER ===\n');
fprintf('Minimum SNR = 10 dB → Baseline Power = %.4f W\n', P_comm_baseline);

% Buffers
P = ones(B,1)*Pt/B;
P_hist = zeros(B,T); r_hist = zeros(B,T); motionFlag = zeros(B,T);
P_comm_log = zeros(1, T); P_sense_log = zeros(1, T);
smooth_r = zeros(B,1); betaSmooth = 0.7; alphaP = 0.8;

for t = 1:T
    SNR_min_t = 10^(SNR_dB_dynamic(t)/10);
    P_comm = min(Pt, SNR_min_t / gain_90);
    P_sense = Pt - P_comm;
    P_comm_log(t) = P_comm;
    P_sense_log(t) = P_sense;

    h = (randn(K,1) + 1j * randn(K,1)) / sqrt(2);
    G = zeros(B,K);
    for b = 1:B
        ab = steer(beamAz(b));
        for k = 1:K
            ak = steer(objAng(k,t));
            G(b,k) = abs(ab' * ak)^2;
        end
    end

    r = zeros(B,1);
    for b = 1:B
        signal = 0;
        for k = 1:K
            signal = signal + sqrt(P(b)) * sqrt(G(b,k)) * sqrt(RCS(k)) * h(k);
        end
        r(b) = abs(signal)^2 + sigma2n/2 * randn;
    end

    smooth_r = betaSmooth * smooth_r + (1 - betaSmooth) * r;
    if t > 1
        delta = abs(smooth_r - r_hist(:,t-1));
        thresh = 3 * sqrt(sigma2n);
        motionFlag(:,t) = delta > thresh;
    end

    snrEst = max(r - sigma2n, 0);
    if sum(snrEst)==0, snrEst = ones(B,1); end
    Pnew = P_sense * snrEst / sum(snrEst);
    P = alphaP * P + (1 - alphaP) * Pnew;
    P_hist(:,t) = P;
    r_hist(:,t) = smooth_r;
end

%% === PLOTS ===
tAxis = 1:T;

figure;
plot(tAxis, 10*log10(r_hist'+1e-8),'LineWidth',1.2); grid on;
ylabel('Received Power (dB)'); title('Phase 3: Smoothed Return Envelopes');

figure;
stem(tAxis, motionFlag','filled'); grid on;
xlabel('Frame Index'); ylabel('Motion Flag');
yticks([0 1]); ylim([-0.2 1.2]);
title('Phase 3: Detected Movement Events');

% === Communication vs Sensing Power Plot
figure;
plot(tAxis, P_comm_log, 'r-', 'LineWidth', 2); hold on;
plot(tAxis, P_sense_log, 'b--', 'LineWidth', 2);
yline(P_comm_baseline, 'k--', 'Baseline Comm Power (10 dB)', 'LineWidth', 1.5, ...
       'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle');
xlabel('Time Frame'); ylabel('Power (W)');
title('Communication vs Tracking Power over Time');
legend('Communication Power', 'Tracking Power', 'Location', 'best'); grid on;

% === Detailed Per-Beam Power Plot
figure('Name','Per-Beam Power Allocation');
hold on;
for b = 1:B
    plot(tAxis, P_hist(b,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Sensing Beam %d°', beamAz(b)));
end
plot(tAxis, P_comm_log, 'k--', 'LineWidth', 2.5, ...
     'DisplayName', 'Communication Beam (90°)');
yline(P_comm_baseline, 'k:', 'Baseline Comm Power', 'LineWidth', 1.2);
xlabel('Time Frame'); ylabel('Power (W)');
title('Individual Beam Power Allocation Over Time');
legend('Location', 'eastoutside'); grid on;
