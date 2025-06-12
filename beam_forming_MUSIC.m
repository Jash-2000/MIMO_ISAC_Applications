%% beam_forming_MUSIC.m
% MUSIC AoA estimation + beamforming for a 64-element ULA
% -------------------------------------------------------

clear;  close all;  clc;

%% 1. Array & signal parameters
fc      = 2.4e9;          % carrier (Hz)
c       = 3e8;            % speed of light (m/s)
lambda  = c / fc;         % wavelength (m)
d       = lambda / 2;     % half-wavelength element spacing
k       = 2*pi / lambda;  % wave-number (rad/m)

N       = 64;             % # array elements
M       = 200;            % # snapshots
trueAng = [ 10  -25 ];    % test AoAs (deg)
S       = numel(trueAng); % # sources
SNR_dB  = 10;             % per-source SNR (dB)

%% 2. Snapshot synthesis:  x = A * s + n
angles_rad = deg2rad(trueAng);
A = exp(1j * k * d * (0:N-1).' .* sin(angles_rad));   % N×S steering matrix
sig = (randn(S,M) + 1j*randn(S,M)) / sqrt(2);         % i.i.d. complex Gauss
x_clean = A * sig;

sigma2 = 10^(-SNR_dB/10);                             % noise variance
noise  = sqrt(sigma2) * (randn(N,M) + 1j*randn(N,M)) / sqrt(2);
x = x_clean + noise;                                  % received snapshots

%% 3. MUSIC algorithm
R = (x * x') / M;                 % sample covariance (N×N)
[Evecs,Eval]  = eig(R);           % eigen-decomposition
[evals,idx]   = sort(real(diag(Eval)),'descend');
U_n = Evecs(:, idx(S+1:end));     % noise subspace (N×(N-S))

theta_scan = -90 : 0.1 : 90;      % search grid (deg)
music_spec = zeros(size(theta_scan));

for ii = 1:numel(theta_scan)
    a_th   = exp(1j * k * d * (0:N-1).' * sin(deg2rad(theta_scan(ii))));
    denom  = a_th' * (U_n * U_n') * a_th;
    music_spec(ii) = 1 ./ abs(denom);     % |…| ensures real-positive
end

% Convert to dB and normalise
music_spec_dB = 10 * log10( music_spec / max(music_spec) );

% --- Peak picking (find S largest peaks) -------------------------------
[~, peakIdx] = findpeaks(music_spec_dB, ...
                         'SortStr','descend', ...
                         'NPeaks',  S, ...
                         'MinPeakDistance', 5);   % avoid very close peaks
estAng = sort(theta_scan(peakIdx));

%% 4. Form beams towards estimated AoAs
w = exp(1j * k * d * (0:N-1).' .* sin(deg2rad(estAng)));   % N×S weights
w = w ./ vecnorm(w);                                       % unit-norm cols

phi = -90 : 0.1 : 90;                                      % evaluation grid
pattern = zeros(size(phi));
for jj = 1:S
    beam_field = w(:,jj)' * exp(1j * k * d * (0:N-1).' * sin(deg2rad(phi)));
    pattern = pattern + abs(beam_field).^2;                % power pattern
end
pattern_dB = 10 * log10( pattern / max(pattern) );

%% 5. Visualisation
figure;
plot(theta_scan, music_spec_dB, 'LineWidth', 1.3);
hold on;
stem(estAng, max(music_spec_dB)*ones(size(estAng)), 'r', 'filled');
hold off;
xlabel('\theta (deg)'); ylabel('MUSIC spectrum (dB)');
title('MUSIC Pseudospectrum'); grid on;

figure;
plot(phi, pattern_dB, 'LineWidth', 1.3);
xlabel('\phi (deg)'); ylabel('Beam pattern (dB)');
title('Composite Beam Pattern steered to MUSIC estimates'); grid on;

%% 6. Console output
fprintf('\nTrue AoAs (deg):      %s\n', mat2str(trueAng));
fprintf('Estimated AoAs (deg): %s\n', mat2str(estAng));
