% ============================================================
% Run CELP codecs: 16 kbps, 11 kbps, 9.6 kbps
% On your own real female voice WAV file
% ============================================================

clc; clear; close all;

%% ----------------------------------------------------------
% 1. Load YOUR speech file (replace filename here)
% ----------------------------------------------------------

input_file = '/MATLAB Drive/FileExchange/CELP codec-1.0.0.0/CELP_done/female_voice_8k.wav';   % <-- put your file here
fprintf("Loading: %s\n", input_file);

[x, Fs] = audioread(input_file);
x = x(:);   % mono

%% ----------------------------------------------------------
% 2. CELP REQUIRES 8000 Hz — resample if needed
% ----------------------------------------------------------
if Fs ~= 8000
    fprintf("Resampling %d → 8000 Hz for CELP compatibility…\n", Fs);
    x = resample(x, 8000, Fs);
    Fs = 8000;
end

% Save a clean 8 kHz copy for reference
audiowrite('female_voice_8k.wav', x, Fs);

%% ----------------------------------------------------------
% 3. CELP parameters (DO NOT CHANGE)
% ----------------------------------------------------------
N = 160;      % Frame length = 20 ms @ 8k
L = 40;       % Subframe length = 5 ms
M = 12;       % LPC order
c = 0.85;     % perceptual weight
Pidx = [16 160];

%% Gaussian codebook
rng(0);
cb = randn(L, 1024);

%% ----------------------------------------------------------
% 4. Run CELP codecs
% ----------------------------------------------------------
tic
fprintf("Running 11 kbps CELP...\n");
[xhat_11000, ~, ~, ~, ~, ~] = celp11k(x, N, L, M, c, cb, Pidx);

runtime_sec = toc;
%% ----------------------------------------------------------
% 5. Playback
% ----------------------------------------------------------
disp('Playing original audio...');
sound(x, Fs); pause(length(x)/Fs + 0.5);

disp('Playing CELP 11 kbps...');
sound(xhat_11000, Fs); pause(length(xhat_11000)/Fs + 0.5);

%% ----------------------------------------------------------
% 6. Plot waveforms
% ----------------------------------------------------------
figure;
subplot(4,1,1); plot(x); title('Original (8 kHz)');

subplot(4,1,3); plot(xhat_11000); title('CELP 11 kbps');


%% ----------------------------------------------------------
% 7. Comparison plots
% ----------------------------------------------------------


figure; plot([x; NaN; xhat_11000]);
title('Original vs CELP 11 kbps');
legend('Original','11 kbps');


disp('Done.');


fprintf("Runtime of encoder   : %.3f seconds\n", runtime_sec);