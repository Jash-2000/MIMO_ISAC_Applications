% ------------------------------------------------------------
% gesture_beam_demo.m  --  Adaptive beam power + gesture tagging
% ------------------------------------------------------------
clear; clc; close all;

%% PARAMETERS -------------------------------------------------
fc   = 2.4e9;      c = 3e8;   lambda = c/fc;   d = lambda/2;
Nt   = 8;          Nr = 8;    B = 5;           beamAz = [-60 -30 0 30 60];
Ptot = 1;                                           % total TX power (W)
Fs   = 1000;                                        % "fast-time" sample rate
N    = 4096;                                       % total fast-time samples
adaptInt = 64;                                     % re-allocate every 64 samp

% Object trajectories (deg vs time index) & gesture labels
K = 3;
objAng = [ linspace(-30,0,N) ;                      % Obj-1 moves
           30*ones(1,N) ;                           % Obj-2 static
           60*ones(1,N) ];                          % Obj-3 static
gesture = ["waving","static","static"];             % assign gestures
RCS     = [1, 0.8, 1.2];                            % relative RCS
sigma2n = 1e-3;                                     % noise power

% Doppler model for "waving"
fd_max   = 80;      fm = 1;                         % 1 Hz hand motion

%% PRE-COMPUTE beam & steering vectors ------------------------
steer = @(th) exp(1j*2*pi*d/lambda*(0:Nt-1).'*sind(th))/sqrt(Nt);
aBeam = zeros(Nt,B);  
for b=1:B, aBeam(:,b)=steer(beamAz(b)); end

%% STATE ------------------------------------------------------
P      = ones(B,1)*Ptot/B;          % start equal power
yBeam  = zeros(B,N);                % raw beamformed base-band
rWin   = zeros(B,adaptInt);         % window buffer for SNR
ptr    = 1;                         % pointer inside window

%% MAIN FAST-TIME LOOP ----------------------------------------
for n = 1:N
    % Superpose object echoes into every beam ****************************
    ySnap = zeros(B,1);                                 % B×1 for this sample
    for k = 1:K
        ak = steer(objAng(k,n));
        
        % Gesture-dependent micro-Doppler
        switch gesture(k)
            case "static"
                sk = 1;
            case "waving"
                fd = fd_max * sin(2*pi*fm*n/Fs);        % time-varying Doppler
                sk = exp(1j*2*pi*fd*(n/Fs));
        end
        
        h = (randn + 1j*randn)/sqrt(2);                 % Rayleigh each sample
        for b = 1:B
            gain = abs(aBeam(:,b)'*ak)^2;               % |a_b^H a_k|^2
            ySnap(b) = ySnap(b) + sqrt(P(b))*sqrt(gain)* ...
                        sqrt(RCS(k))*h*sk;
        end
    end
    % Add AWGN
    ySnap = ySnap + sqrt(sigma2n/2)*(randn(B,1)+1j*randn(B,1));
    yBeam(:,n) = ySnap;                                 % store
    
    % ========= windowed power for adaptation ===========================
    rWin(:,ptr) = abs(ySnap).^2;
    ptr = ptr + 1;
    
    if ptr > adaptInt        % time to re-allocate power
        snrEst = mean(rWin,2) - sigma2n;
        snrEst = max(snrEst,0);
        if all(snrEst==0), snrEst(:)=1; end
        P  = 0.8*P + 0.2*(Ptot*snrEst/sum(snrEst));  % smoothed update
        
        % reset window
        rWin(:,:) = 0;
        ptr = 1;
    end
end

%% GESTURE CLASSIFICATION PER BEAM -------------------------------------
% STFT settings
wlen = 256;  nover = 192;   nfft = 512;
freqAxis = Fs*(0:(nfft/2))/nfft;

beamGesture = strings(1,B);

for b = 1:B
    [S,F,Tspec] = spectrogram(yBeam(b,:),hann(wlen),nover,nfft,Fs,'centered');
    Pxx = abs(S).^2;
    DopplerPow = mean(Pxx,2);                          % average over slow-time
    
    % crude peak-finder
    [pk,idx] = max(DopplerPow);
    fdPeak = F(idx);
    
    % classify
    if abs(fdPeak) < 10
        beamGesture(b) = "static";
    elseif abs(fdPeak) > 60 && abs(fdPeak) < 100
        beamGesture(b) = "waving";
    else
        beamGesture(b) = "unknown";
    end
    
    % ---------- Plot spectrogram for visual confirmation ------------
    figure(100+b); clf;
    imagesc(Tspec, F, 10*log10(Pxx/max(Pxx(:))));
    axis xy; colormap jet; colorbar;
    xlabel('slow-time (s)'); ylabel('Doppler (Hz)');
    title(sprintf('Beam @ %d°  →  classified as %s', beamAz(b), beamGesture(b)));
end

%% SUMMARY TABLE --------------------------------------------------------
fprintf('\n=== Gesture assignment (dominant object per beam) ===\n');
for b = 1:B
    fprintf('Beam %-2d @ %+3d°  :  %s\n', b, beamAz(b), beamGesture(b));
end
