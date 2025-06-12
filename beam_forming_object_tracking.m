% Dynamically track movement over time using an adaptive beam steering mechanism based on:
% % Time-varying angle of arrival (AoA) (simulated object moving across space).
% % Periodic feedback from the receiver to the transmitter.
% % A simple Kalman-like filter to smooth AoA estimation (optional but improves realism).
% % A continuously updated dual-beam pattern.

clc; clear;

%% Parameters
M = 64;          % Transmit antennas
K = 5;           % Number of users
fc = 28e9;       % Carrier frequency
c = 3e8;
lambda = c/fc;
d = lambda/2;    % Antenna spacing

angles_comm = [-10, -5, 0, 5, 10]; % User angles in degrees
theta_comm_rad = deg2rad(angles_comm);

% Time settings
T = 1;           % Time step (s)
num_steps = 100; % Simulation duration

% Object movement
theta_obj_true = linspace(60, -20, num_steps);  % Moving from 60° to -20°
theta_obj_est = zeros(1, num_steps);            % Receiver's estimate
theta_obj_est(1) = theta_obj_true(1) + randn;   % Initial noisy estimate

% Optional: Kalman-like smoothing
alpha_track = 0.85;

% Steering vector
a_theta = @(theta_deg) exp(1j*2*pi*d*(0:M-1)'*sin(deg2rad(theta_deg))/lambda) / sqrt(M);

%% Channel matrix (fixed)
H = zeros(K, M);
for k = 1:K
    H(k,:) = a_theta(rad2deg(theta_comm_rad(k)))';
end

% Communication beam (fixed)
mean_comm_angle = mean(theta_comm_rad);
a_comm = a_theta(rad2deg(mean_comm_angle));

%% Visualization setup
theta_scan = -90:0.5:90;
beam_pattern_time = zeros(length(theta_scan), num_steps);

%% Main Simulation Loop
for t = 2:num_steps
    % Receiver estimates AoA with noise
    aoa_measured = theta_obj_true(t) + 1.5*randn;
    
    % Apply tracking filter (e.g., low-pass filter)
    theta_obj_est(t) = alpha_track * theta_obj_est(t-1) + (1 - alpha_track) * aoa_measured;
    
    % Construct sensing beam
    a_obj = a_theta(theta_obj_est(t));
    
    % Power allocation
    alpha = sqrt(0.7); beta = sqrt(0.3);
    w_dual = alpha*a_comm + beta*a_obj;
    
    % Beam pattern at this time
    for i = 1:length(theta_scan)
        a_scan = a_theta(theta_scan(i));
        beam_pattern_time(i, t) = abs(a_scan' * w_dual)^2;
    end
end

%% Plotting: Dynamic Beam Pattern as Heatmap
figure;
imagesc(1:num_steps, theta_scan, 10*log10(beam_pattern_time));
xlabel('Time Step');
ylabel('Angle (deg)');
title('Dynamic Dual-Beam Pattern Over Time');
colorbar; colormap jet;
axis xy;
