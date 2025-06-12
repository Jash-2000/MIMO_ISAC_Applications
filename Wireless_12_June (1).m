% ------------------------------------------------------------
% Adaptive-power multi-beam demo with motion detection
% ------------------------------------------------------------
clear; clc; close all;

%% 1) Basic constants
fc      = 2.4e9;                 c = 3e8;        % carrier & light speed
lambda  = c/fc;                  d = lambda/2;   % element spacing
Nt      = 8;                     Nr = 8;         % antennas
B       = 5;                                    % beams
beamAz  = [-60 -30 0 30 60];                    % steering directions (deg)
Ptot    = 1;                                    % total TX power (W)

%% 2) Helper: array steering vector (ULA, half-wave spaced)
steer  = @(th) exp(1j*2*pi*d/lambda*(0:Nt-1).' * sind(th))/sqrt(Nt);

% Pre-compute beamforming gains g(b,k)=|a_t^H(beam)*a_t(obj)|^2
G = zeros(B, 1001);               % we will fill per frame for speed

%% 3) Object trajectories (angles in deg vs frame idx)
T       = 100;                                  % frames
objAng  = [ -30*ones(1,25) , linspace(-30,0,75) ;   % object 1 moves
             30*ones(1,T) ;                       % object 2 static
             60*ones(1,T) ];                      % object 3 static
K       = size(objAng,1);                        % # objects
sigma2n = 1e-3;                                  % noise power

RCS     = [1, 0.8, 1.2];                         % relative radar cross-sections

%% 4) Initialise power allocation & buffers
P   = ones(B,1)*Ptot/B;               % equal start
P_hist = zeros(B,T);                  % store for plotting
r_hist = zeros(B,T);
motionFlag = zeros(B,T);              % 1 where motion detected
alphaP = 0.8;                         % smoothing for power updates
smooth_r = zeros(B,1);  betaSmooth = 0.7;

%% 5) Simulation loop
for t = 1:T
    %% 5-a) Rayleigh channels for all objects this frame
    h = (randn(K,1)+1j*randn(K,1))/sqrt(2);      % K×1
    
    %% 5-b) Build per-beam gains to every object (constant over frame)
    for b = 1:B
        ab  = steer(beamAz(b));
        gSum = 0;
        for k = 1:K
            ak = steer(objAng(k,t));
            G(b,k) = abs(ab' * ak)^2;            % beam gain factor
        end
    end
    
    %% 5-c) Received power in every beam with current P(b)
    r = zeros(B,1);
    for b = 1:B
        signal = 0;
        for k = 1:K
            signal = signal + sqrt(P(b))*sqrt(G(b,k))*sqrt(RCS(k))*h(k);
        end
        r(b) = abs(signal)^2 + sigma2n/2*randn;  % power + AWGN
    end
    
    %% 5-d) Movement detection (simple energy change after smoothing)
    smooth_r = betaSmooth*smooth_r + (1-betaSmooth)*r;
    if t>1
        delta = abs(smooth_r - r_hist(:,t-1));
        thresh = 3*sqrt(sigma2n);                % coarse threshold
        motionFlag(:,t) = delta > thresh;
    end
    
    %% 5-e) Power adaptation: proportional-to-SNR rule
    snrEst = max(r - sigma2n, 0);                % quick-and-dirty SNR proxy
    if sum(snrEst)==0, snrEst = ones(B,1); end   % avoid div-by-zero
    Pnew   = Ptot * snrEst / sum(snrEst);        % redistribute
    P      = alphaP*P + (1-alphaP)*Pnew;         % smooth
       
    %% 5-f) Store histories
    P_hist(:,t) = P;
    r_hist(:,t) = smooth_r;
end

%% 6) Visualisation --------------------------------------------------
tAxis = 1:T;

figure('Position',[100 100 900 600]);

% (a) TX power allocation
subplot(3,1,1);
area(tAxis, P_hist','LineStyle','none'); grid on;
title('Adaptive transmit-power allocation per beam'); ylabel('Power (W)');
legend(arrayfun(@(x) sprintf('%d°',x),beamAz,'uni',0),"Location","eastoutside");

% (b) Smoothed received power envelopes
subplot(3,1,2);
plot(tAxis, 10*log10(r_hist'),'LineWidth',1.3); grid on;
ylabel('Received power (dB)'); title('Per-beam return envelopes');

% (c) Motion detection events
subplot(3,1,3);
stem(tAxis, motionFlag','filled'); grid on;
xlabel('Frame index'); ylabel('Motion flag');
yticks([0 1]); ylim([-0.2 1.2]);
title('Detected movement (1 = yes)');