%Initial MRT Beamforming to communication users.
%A reflection at +50° is detected.
%TX then splits its energy between: mean comm beam and sensing beam toward the object.

clc; clear;

%% Parameters
M = 64;         % Number of TX antennas
K = 5;          % Number of users (Rx)
fc = 28e9;      % Carrier frequency (28 GHz)
c = 3e8;        % Speed of light
lambda = c/fc;
d = lambda/2;   % Antenna spacing

angles_comm = [-10, -5, 0, 5, 10]; % User angles in degrees
theta_comm_rad = deg2rad(angles_comm);

% Generate channel matrix H (Line of Sight only for simplicity)
H = zeros(K, M);
for k = 1:K
    a_tx = exp(j*2*pi*d*(0:M-1)'*sin(theta_comm_rad(k))/lambda) / sqrt(M);
    H(k,:) = a_tx';
end

%% Initial Beamforming (Communication Only)
% Perform MRT (maximum ratio transmission)
T_comm = H' * inv(H*H' + 1e-6*eye(K));  % Precoder
T_comm = T_comm / norm(T_comm,'fro');   % Normalize

% Visualize beam pattern
theta_scan = -90:0.5:90;
a_scan = @(theta) exp(1j*2*pi*d*(0:M-1)'*sin(deg2rad(theta))/lambda) / sqrt(M);

beam_pattern_comm = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a = a_scan(theta_scan(i));
    beam_pattern_comm(i) = norm(a' * T_comm,2)^2;
end

%% Simulate movement detection (e.g., reflection at +50 degrees)
movement_detected = true;
theta_obj = 50; % Detected object AoA [deg]
theta_obj_rad = deg2rad(theta_obj);

%% Adaptive Beamforming: Dual Beam
% 1 beam to communication users, 1 beam to detected object

% Communication beam (average direction)
mean_comm_angle = mean(theta_comm_rad);
a_comm = exp(1j*2*pi*d*(0:M-1)'*sin(mean_comm_angle)/lambda) / sqrt(M);

% Sensing beam (toward object)
a_obj = exp(1j*2*pi*d*(0:M-1)'*sin(theta_obj_rad)/lambda) / sqrt(M);

% Power allocation (e.g., 70% comm, 30% sensing)
alpha = sqrt(0.7);
beta = sqrt(0.3);
w_dual = alpha*a_comm + beta*a_obj;

%% Visualize updated beam pattern
beam_pattern_dual = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a = a_scan(theta_scan(i));
    beam_pattern_dual(i) = abs(a' * w_dual)^2;
end

%% Plot results
figure;
plot(theta_scan, 10*log10(beam_pattern_comm), 'b--', 'LineWidth', 2); hold on;
plot(theta_scan, 10*log10(beam_pattern_dual), 'r-', 'LineWidth', 2);
xlabel('Angle (deg)'); ylabel('Beam Gain (dB)');
legend('Initial Beam (Comm Only)', 'Dual Beam (Comm + Sensing)');
title('Adaptive Beam Steering in ISAC');
grid on;
