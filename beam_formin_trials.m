clc ; close all; clear all;

M = 64;                     % Number of antennas
fc = 28e9;                  % Carrier frequency
c = 3e8; lambda = c/fc;     
d = lambda/2;               % Spacing
theta = 30;                 % Beam angle in degrees
k = 2*pi/lambda;

% Beamforming (steering) vector
w = exp(j * k * d * (0:M-1)' * sin(deg2rad(theta))) / sqrt(M);


theta_scan = -90:0.5:90;
beam_pattern = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a = exp(j * k * d * (0:M-1)' * sin(deg2rad(theta_scan(i)))) / sqrt(M);
    beam_pattern(i) = abs(a' * w).^2;
end

plot(theta_scan, 10*log10(beam_pattern));
xlabel('Angle (deg)'); ylabel('Gain (dB)');
title('Beam Pattern');
grid on;
