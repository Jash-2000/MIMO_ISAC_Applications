clear all ; close all ; clc ; 

% Parameters
N = 64;                      % number of ULA elements
fc = 2.4e9;                  % carrier frequency (Hz)
c = 3e8;                     % speed of light (m/s)
lambda = c/fc;               % wavelength (m)
d = lambda/2;                % element spacing (half-wavelength)
beam_angles_initial = 0:45:360;   % initial central angles of 5 beams (degrees)
T = 10;                      % total rotation period (seconds)
omega = 360/T;               % rotation speed (degrees per second)
time = 0:0.1:T;              % time samples from 0 to 10 s
angles = 0:1:359;            % angles in degrees to evaluate (0-359)
lambda = c / fc;        % wavelength (m)
d      = lambda / 2;     % half-wavelength spacing
k      = 2*pi / lambda;  % <-- wave-number (rad/m)

% Preallocate results
num_times = length(time);
num_angles = length(angles);
signal_power = zeros(num_times, num_angles);   % will hold received power (without fading)
faded_power = zeros(num_times, num_angles);    % received power including Rayleigh fading

% Simulation loop
for ti = 1:num_times
    t = time(ti);
    % Compute beam steering angles at time t
    current_beams = mod( beam_angles_initial + omega*t, 360 );
    % Loop over angles to compute array gain
    for ai = 1:num_angles
        phi = angles(ai);
        total_gain = 0;  % linear power gain (array factor)
        % Sum contributions from each of the 5 beams
        for bi = 1:length(current_beams)
            theta_b = current_beams(bi);
            % Calculate array factor for beam bi at angle phi:
            % AF = sum_{n=0}^{N-1} exp(j * k * n * d * (cos(phi) - cos(theta_b)) )
            % (phi and theta_b converted to radians for cosine)
            phi_rad = deg2rad(phi);
            theta_b_rad = deg2rad(theta_b);
            % Sum over elements n=0...N-1:
            element_phases = 0:N-1;
            % Compute phase difference for each element
            phase_diff = k * element_phases * d * (cos(phi_rad) - cos(theta_b_rad));
            % Element contributions and beam sum
            field_sum = sum( exp(1j * phase_diff) );
            % Optionally normalize each beam's weight so each beam has unit power:
            field_sum = field_sum / sqrt(N);  
            total_gain = total_gain + abs(field_sum)^2;  % add beam's power
        end
        signal_power(ti, ai) = total_gain;
        % Apply Rayleigh fading (random complex gain with unit variance)
        % Generate one random complex coefficient for this angle/time
        h = (randn() + 1j*randn())/sqrt(2);      % Rayleigh fading coefficient
        faded_power(ti, ai) = total_gain * (abs(h)^2);
    end
end

% Visualization 1: Beam pattern at t = 0 (initial orientation)
figure;
gain_dB = 10*log10(signal_power(1, :) + eps);
plot(angles, gain_dB, 'LineWidth', 1.5);
xlabel('Angle (degrees)'); ylabel('Array gain (dB)');
title('Initial 6-Beam Pattern (t = 0)');
grid on;

% Visualization 2: Signal strength over time and angle
figure;
% Plot a color map of power vs angle vs time
imagesc(angles, time, 10*log10(faded_power + eps));
xlabel('Angle (degrees)'); ylabel('Time (s)');
title('Signal Power (dB) over Time and Angle');
colormap jet; colorbar; caxis([-20 20]);  % color axis limits for better contrast
