%% ========================================================================
%   Multi-File Voice-Excited LPC Vocoder (Residual + DCT)
%   Processes 4 input speech samples and produces 4 output WAV files
%   Computes SEGSNR for each pair
% ========================================================================

clear all; close all; clc;

%% ---------------------------- INPUT FILES --------------------------------
input_files = {
    '/MATLAB Drive/assignment_2/Speech_sample/male_voice_1_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/male_voice_2_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/female_voice_1_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/female_voice_2_8k.wav'
};

output_prefix = 've_lpc_';   % Default prefix

%% ---------------------- LPC–DCT VOCODER PARAMETERS -----------------------
sr      = 16000;       % Processing rate
L       = 10;          % LPC order
fr      = 20;          % frame rate (ms)
fs      = 30;          % frame size (ms)
preemp  = 0.9378;      % Pre-emphasis
frame_samples = round(sr * fs / 1000);

% DCT residual coding parameters
K        = 27;         % keep K DCT coefficients
dct_bits = 7;          
parcor_bits = 6;
gain_bits   = 8;

%% ---------------------- START TIMER -----------------------
tic;

%% ----------------- STORAGE FOR ALL FILES ------------------------
speech_all = cell(4,1);
aCoeff_all = cell(4,1);
resid_all = cell(4,1);
pitch_all = cell(4,1);
G_all = cell(4,1);
parcor_all = cell(4,1);
stream_all = cell(4,1);

%% ---------------------------- READ ALL FILES ---------------------------------
fprintf("\n========================================================\n");
fprintf("Reading all 4 input files...\n");
fprintf("========================================================\n");

for f = 1:4
    fprintf("Reading File %d: %s\n", f, input_files{f});
    [speech, fs_original] = audioread(input_files{f});
    speech = speech(:,1);

    if fs_original ~= sr
        fprintf("  Resampling %d → %d Hz\n", fs_original, sr);
        speech = resample(speech, sr, fs_original);
    end
    
    speech_all{f} = speech;
end

%% ---------------------------- LPC ANALYSIS FOR ALL FILES --------------------------------
fprintf("\n========================================================\n");
fprintf("Performing LPC analysis on all files...\n");
fprintf("========================================================\n");

for f = 1:4
    fprintf("Analyzing File %d...\n", f);
    [aCoeff, resid, pitch, G, parcor, stream] = ...
        proclpc(speech_all{f}, sr, L, fr, fs, preemp);
    
    aCoeff_all{f} = aCoeff;
    resid_all{f} = resid;
    pitch_all{f} = pitch;
    G_all{f} = G;
    parcor_all{f} = parcor;
    stream_all{f} = stream;
end

%% ---------------------------- DCT ENCODING FOR ALL FILES --------------------------
fprintf("\n========================================================\n");
fprintf("DCT encoding for all files...\n");
fprintf("========================================================\n");

dct_coeffs_all = cell(4,1);
for f = 1:4
    fprintf("DCT encoding File %d...\n", f);
    nframes = size(aCoeff_all{f},2);
    residual_frames = resid_all{f};
    
    dct_coeffs = zeros(K, nframes);
    for i = 1:nframes
        tmp = dct(residual_frames(:,i));
        dct_coeffs(:,i) = tmp(1:K);
    end
    dct_coeffs_all{f} = dct_coeffs;
end

%% ------------------------ QUANTIZATION FOR ALL FILES -------------------------------
fprintf("\n========================================================\n");
fprintf("Quantizing all files...\n");
fprintf("========================================================\n");

parcor_q_all = cell(4,1);
dct_q_all = cell(4,1);
G_q_all = cell(4,1);
dct_min_all = zeros(4,1);
dct_max_all = zeros(4,1);
G_min_all = zeros(4,1);
G_max_all = zeros(4,1);

parcor_min = -1; parcor_max = 1;
parcor_levels = 2^parcor_bits;
dct_levels = 2^dct_bits;
G_levels = 2^gain_bits;

for f = 1:4
    fprintf("Quantizing File %d...\n", f);
    
    % PARCOR quantization
    parcor_q = round((parcor_all{f} - parcor_min)/(parcor_max-parcor_min)*(parcor_levels-1));
    parcor_q = min(max(parcor_q,0),parcor_levels-1);
    parcor_q_all{f} = parcor_q;
    
    % DCT quantization
    dct_min = min(dct_coeffs_all{f}(:));  
    dct_max = max(dct_coeffs_all{f}(:))+eps;
    dct_q = round( (dct_coeffs_all{f}-dct_min)/(dct_max-dct_min)*(dct_levels-1) );
    dct_q = min(max(dct_q,0),dct_levels-1);
    dct_q_all{f} = dct_q;
    dct_min_all(f) = dct_min;
    dct_max_all(f) = dct_max;
    
    % Gain quantization
    G_log = log(G_all{f}+eps);
    G_min = min(G_log); 
    G_max = max(G_log)+eps;
    G_q = round((G_log-G_min)/(G_max-G_min)*(G_levels-1));
    G_q_all{f} = G_q;
    G_min_all(f) = G_min;
    G_max_all(f) = G_max;
end

%% --------------------- DECODER: DEQUANTIZATION FOR ALL FILES -----------------------
fprintf("\n========================================================\n");
fprintf("Dequantizing all files...\n");
fprintf("========================================================\n");

aCoeff_dec_all = cell(4,1);
resid_dec_all = cell(4,1);
G_dec_all = cell(4,1);

for f = 1:4
    fprintf("Dequantizing File %d...\n", f);
    
    parcor_dec = parcor_q_all{f}/(parcor_levels-1)*(parcor_max-parcor_min)+parcor_min;
    dct_dec    = dct_q_all{f}/(dct_levels-1)*(dct_max_all(f)-dct_min_all(f))+dct_min_all(f);
    G_log_dec  = G_q_all{f}/(G_levels-1)*(G_max_all(f)-G_min_all(f))+G_min_all(f);
    G_dec      = exp(G_log_dec);
    G_dec_all{f} = G_dec;
    
    % Convert PARCOR → LPC
    nframes = size(aCoeff_all{f},2);
    aCoeff_dec = zeros(L+1, nframes);
    for i = 1:nframes
        aCoeff_dec(:,i) = rc2poly(parcor_dec(:,i));
    end
    aCoeff_dec_all{f} = aCoeff_dec;
    
    % Reconstruct residual
    resid_dec = zeros(frame_samples,nframes);
    for i = 1:nframes
        tmp = zeros(frame_samples,1);
        tmp(1:K) = dct_dec(:,i);
        resid_dec(:,i) = idct(tmp);
    end
    resid_dec_all{f} = resid_dec;
end

%% ------------------------ LPC SYNTHESIS FOR ALL FILES ------------------------------
fprintf("\n========================================================\n");
fprintf("Synthesizing speech for all files...\n");
fprintf("========================================================\n");

synth_all = cell(4,1);
for f = 1:4
    fprintf("Synthesizing File %d...\n", f);
    synth = synlpc(aCoeff_dec_all{f}, resid_dec_all{f}, sr, G_dec_all{f}, fr, fs, preemp);
    synth = synth / max(abs(synth)+eps);
    synth_all{f} = synth;
end

%% ------------------------- SAVE ALL OUTPUTS -------------------------------
fprintf("\n========================================================\n");
fprintf("Saving all output files...\n");
fprintf("========================================================\n");

output_files = cell(4,1);
seg_values = zeros(4,1);

for f = 1:4
    outfile = sprintf('%soutput_%d.wav', output_prefix, f);
    audiowrite(outfile, synth_all{f}, sr);
    output_files{f} = outfile;
    
    % SEGSNR computation
    Lmin = min(length(speech_all{f}), length(synth_all{f}));
    seg_values(f) = segsnr(speech_all{f}(1:Lmin), synth_all{f}(1:Lmin), sr);
    
    fprintf("File %d: %s (SEGSNR = %.3f dB)\n", f, outfile, seg_values(f));
end

%% ---------------------- STOP TIMER AND CALCULATE RUNTIME -----------------------
total_runtime = toc;

%% ========================================================================
%                    BIT-RATE CALCULATION
% ========================================================================

fprintf("\n================= BIT-RATE ESTIMATION ==================\n");

% Calculate for each file
bitrates = zeros(4,1);
durations = zeros(4,1);

for f = 1:4
    nframes = size(aCoeff_all{f}, 2);
    
    % Bits per frame calculation
    parcor_bits_per_frame = L * parcor_bits;          % PARCOR coefficients
    dct_bits_per_frame = K * dct_bits;                % DCT coefficients
    gain_bits_per_frame = gain_bits;                  % Gain
    pitch_bits_per_frame = 8;                         % Pitch (assuming 8 bits)
    
    total_bits_per_frame = parcor_bits_per_frame + dct_bits_per_frame + ...
                           gain_bits_per_frame + pitch_bits_per_frame;
    
    % Total bits for the file
    total_bits = total_bits_per_frame * nframes;
    
    % Duration of speech
    duration = length(speech_all{f}) / sr;
    durations(f) = duration;
    
    % Bit-rate in kbps
    bitrate = total_bits / duration / 1000;
    bitrates(f) = bitrate;
    
    fprintf("File %d:\n", f);
    fprintf("  Frames: %d\n", nframes);
    fprintf("  Duration: %.2f seconds\n", duration);
    fprintf("  Bits per frame: %d (PARCOR:%d + DCT:%d + Gain:%d + Pitch:%d)\n", ...
            total_bits_per_frame, parcor_bits_per_frame, dct_bits_per_frame, ...
            gain_bits_per_frame, pitch_bits_per_frame);
    fprintf("  Total bits: %d\n", total_bits);
    fprintf("  Bit-rate: %.2f kbps\n\n", bitrate);
end

avg_bitrate = mean(bitrates);
fprintf("Average Bit-rate across all files: %.2f kbps\n", avg_bitrate);

%% ========================================================================
%                             RESULTS TABLE
% ========================================================================

fprintf("\n================= FINAL RESULTS ==================\n");

T = table(input_files, output_files, seg_values, bitrates, durations, ...
    'VariableNames', {'Original_File','Output_File','SEGSNR_dB','Bitrate_kbps','Duration_sec'});

disp(T);

fprintf("\n================= RUNTIME ANALYSIS ==================\n");
fprintf("Total Runtime: %.3f seconds\n", total_runtime);
fprintf("Average time per file: %.3f seconds\n", total_runtime/4);

fprintf("\nProcessing complete for all 4 files.\n");