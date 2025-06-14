clc; clear;

% === System Parameters ===
M = 64;
fc = 2.4e9;
c = 3e8;
lambda = c / fc;
d = lambda / 2;
Pt = 1;

% === SNR Threshold ===
SNR_min_db = 10;
SNR_min = 10^(SNR_min_db / 10);

% === Simulation Parameters ===
num_paths = 40;
num_runs = 1000;
angles_deg = 0:359;
N_values = zeros(num_runs, 1);
rng(1);
for run = 1:num_runs
    % === Generate Rayleigh multipath channel ===
    h_total = zeros(M,1);
    for p = 1:num_paths
        phi = 2 * pi * rand();
        gain = (randn() + 1j * randn()) / sqrt(2);
        a_phi = exp(1j * 2 * pi * d * (0:M-1)' * sin(phi));
        h_total = h_total + gain * a_phi;
    end

    % === Compute gain over 360° sweep ===
    gain_values = zeros(size(angles_deg));
    for idx = 1:length(angles_deg)
        theta_rad = deg2rad(angles_deg(idx));
        a_theta = exp(1j * 2 * pi * d * (0:M-1)' * sin(theta_rad));
        h_proj = (h_total' * a_theta) / norm(a_theta);
        gain_values(idx) = abs(h_proj)^2;
    end

    % === Focus on user region: 45°–135° ===
    user_sector_idx = (angles_deg >= 45 & angles_deg <= 135);
    user_sector_gains = gain_values(user_sector_idx);
    min_gain = min(user_sector_gains);

    % === Calculate number of beams ===
    N_values(run) = floor(min_gain / SNR_min);
end

% === Output ===
avg_N = mean(N_values);
fprintf('===== Average over %d Monte Carlo Runs =====\n', num_runs);
fprintf('Average N (beams possible) : %.2f\n', avg_N);
fprintf('Max N observed             : %d\n', max(N_values));
fprintf('Min N observed             : %d\n', min(N_values));
fprintf('===========================================\n');

% Optional histogram
figure;
histogram(N_values, 'BinMethod', 'integers');
xlabel('Number of Beams (N)');
ylabel('Frequency');
title('Distribution of N over 100 Monte Carlo Runs');
grid on;
