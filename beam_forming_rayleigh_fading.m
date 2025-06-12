clc; clear;

%% Parameters
M = 64;            % Transmit antennas
K = 5;             % Users
fc = 28e9;         % Carrier frequency
c = 3e8; lambda = c/fc;
d = lambda/2;      % Antenna spacing
T = 1;             % Time step
num_steps = 100;   % Duration

% Communication users' angles
angles_comm = [-10, -5, 0, 5, 10];
theta_comm_rad = deg2rad(angles_comm);

% Object movement
theta_obj_true = linspace(60, -20, num_steps);  % Object trajectory
theta_obj_est = zeros(1, num_steps);
theta_obj_est(1) = theta_obj_true(1) + randn;

% Kalman-style smoothing factor
alpha_track = 0.85;

% Steering vector function
a_theta = @(theta_deg) exp(1j*2*pi*d*(0:M-1)'*sin(deg2rad(theta_deg))/lambda) / sqrt(M);

% Large-scale fading (static)
user_dists = 50 + (200-50)*rand(K,1);
shadow = 8*randn(K,1);
PL_dB = -128.1 - 36.7*log10(user_dists/1000) + shadow;
xi = 10.^(PL_dB/10);  % path loss gain

%% Preallocate
beam_pattern_time = zeros(361, num_steps);
theta_scan = -90:90;

% Loop over time
for t = 1:num_steps
    %% Channel fading (Rayleigh, time-varying)
    H_t = zeros(K, M);
    for k = 1:K
        h_k = (randn(1,M) + 1j*randn(1,M))/sqrt(2);  % Rayleigh fading
        a_k = a_theta(rad2deg(theta_comm_rad(k)))';
        H_t(k,:) = sqrt(xi(k)) * h_k .* a_k;        % Combine fading + AoA
    end

    %% MRT beamformer for communication users
    T_comm = H_t' * inv(H_t*H_t' + 1e-6*eye(K));
    T_comm = T_comm / norm(T_comm, 'fro');  % Normalize

    %% Mean direction for communication
    mean_comm_angle = mean(angles_comm);
    a_comm = a_theta(mean_comm_angle);

    %% Movement detection and tracking
    if t > 1
        aoa_measured = theta_obj_true(t) + 1.5*randn;
        theta_obj_est(t) = alpha_track * theta_obj_est(t-1) + ...
                          (1 - alpha_track) * aoa_measured;
    else
        theta_obj_est(t) = theta_obj_true(t);
    end

    %% Dual beam (comm + sensing)
    a_obj = a_theta(theta_obj_est(t));
    alpha = sqrt(0.7); beta = sqrt(0.3);
    w_dual = alpha*a_comm + beta*a_obj;

    %% Compute beam pattern
    for i = 1:length(theta_scan)
        a_scan = a_theta(theta_scan(i));
        beam_pattern_time(i, t) = abs(a_scan' * w_dual)^2;
    end
end

%% Plot beam pattern evolution
figure;
imagesc(1:num_steps, theta_scan, 10*log10(beam_pattern_time));
xlabel('Time Step'); ylabel('Angle (deg)');
title('Beam Pattern Over Time with Fading and Tracking');
colorbar; colormap jet;
axis xy;