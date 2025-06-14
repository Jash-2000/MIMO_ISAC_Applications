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

% === Monte Carlo Parameters ===
num_paths = 40;
num_runs = 1000;
angles_deg = 0:359;
N_values = zeros(num_runs, 1);
rng(1);

for run = 1:num_runs
    h_total = zeros(M,1);
    for p = 1:num_paths
        phi = 2 * pi * rand();
        gain = (randn() + 1j * randn()) / sqrt(2);
        a_phi = exp(1j * 2 * pi * d * (0:M-1)' * sin(phi));
        h_total = h_total + gain * a_phi;
    end

    gain_values = zeros(size(angles_deg));
    for idx = 1:length(angles_deg)
        theta_rad = deg2rad(angles_deg(idx));
        a_theta = exp(1j * 2 * pi * d * (0:M-1)' * sin(theta_rad));
        h_proj = (h_total' * a_theta) / norm(a_theta);
        gain_values(idx) = abs(h_proj)^2;
    end

    % Focus on user region: 45°–135°
    user_sector_idx = (angles_deg >= 45 & angles_deg <= 135);
    min_gain = min(gain_values(user_sector_idx));

    N_values(run) = floor(min_gain / SNR_min);
end

% === Output Monte Carlo Summary ===
avg_N = mean(N_values);
fprintf('===== Monte Carlo Result over %d Runs =====\n', num_runs);
fprintf('Average N (Beams) : %.2f\n', avg_N);
fprintf('Max N             : %d\n', max(N_values));
fprintf('Min N             : %d\n', min(N_values));
fprintf('==========================================\n');

% === Plot Histogram of N ===
figure;
histogram(N_values, 'BinMethod', 'integers');
xlabel('Number of Beams (N)');
ylabel('Frequency');
title('Distribution of N over 1000 Monte Carlo Runs');
grid on;

% === Visualization of Power vs Time and Angle ===
% Assume average N beams rotate over 10 seconds
N = round(avg_N);
beam_angles_initial = linspace(0, 360, N+1); beam_angles_initial(end) = [];

T = 15;                        % total time (sec)
omega = 360 / T;              % rotation speed (deg/sec)
time = 0:0.1:T;               % time steps
angles = 0:1:359;             % sweep angles
k = 2*pi / lambda;
num_times = length(time);
num_angles = length(angles);
faded_power = zeros(num_times, num_angles);

for ti = 1:num_times
    t = time(ti);
    current_beams = mod(beam_angles_initial + omega * t, 360);

    for ai = 1:num_angles
        phi = angles(ai);
        total_gain = 0;

        for bi = 1:length(current_beams)
            theta_b = current_beams(bi);
            phi_rad = deg2rad(phi);
            theta_b_rad = deg2rad(theta_b);
            n = 0:M-1;
            phase_shift = k * n * d * (cos(phi_rad) - cos(theta_b_rad));
            field = sum(exp(1j * phase_shift)) / sqrt(M);
            total_gain = total_gain + abs(field)^2;
        end

        h = (randn() + 1j * randn()) / sqrt(2);  % Rayleigh fading
        faded_power(ti, ai) = total_gain * abs(h)^2;
    end
end

% === Visualization ===
figure;
imagesc(angles, time, 10*log10(faded_power + eps));
xlabel('Angle (degrees)');
ylabel('Time (s)');
title(sprintf('Signal Power (dB) over Time and Angle [N = %.0f Beams]', avg_N));
colormap jet; colorbar;
caxis([-20 20]);  % adjust for better contrast
