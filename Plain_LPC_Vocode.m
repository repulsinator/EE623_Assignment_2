%% Plain LPC Vocoder (LPC-10) - Batch Processing for 4 Files

clear all; close all; clc;

%% ============== INPUT FILES ==============
input_files = {
    '/MATLAB Drive/assignment_2/Speech_sample/male_voice_1_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/male_voice_2_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/female_voice_1_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/female_voice_2_8k.wav'
};

output_prefix = 'plain_lpc_';
bitstream_prefix = 'plain_lpc_bitstream_';

%% ============== PARAMETERS ==============
% LPC Parameters
sr = 8000;              % Sampling rate (Hz) - 8kHz for narrowband
L = 15;                 % LPC order (10 coefficients for LPC-10) Modified to 15 
fr = 20;                % Frame rate (ms)
fs = 30;                % Frame size (ms)
preemp = 0.9378;        % Pre-emphasis coefficient

% Quantization Parameters (based on uniform/non-uniform quantization theory)
parcor_bits = 5;        % Bits per PARCOR coefficient (total: 10*5 = 50 bits)
pitch_bits = 6;         % Bits for pitch (0-63 range)
gain_bits = 5;          % Bits for log gain (non-uniform quantization)
vuv_bits = 1;           % Voiced/Unvoiced flag

total_bits_per_frame = L*parcor_bits + pitch_bits + gain_bits + vuv_bits;

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

%% ============== READ ALL FILES ==============
fprintf('\n========================================================\n');
fprintf('Reading all 4 input files...\n');
fprintf('========================================================\n');

for f = 1:4
    fprintf('Reading File %d: %s\n', f, input_files{f});
    [speech, fs_original] = audioread(input_files{f});
    speech = speech(:,1);  % Convert to mono if stereo

    % Resample if necessary
    if fs_original ~= sr
        speech = resample(speech, sr, fs_original);
        fprintf('  Resampled from %d Hz to %d Hz\n', fs_original, sr);
    end
    
    speech_all{f} = speech;
end

%% ============== ENCODER (ANALYSIS) FOR ALL FILES ==============
fprintf('\n========================================================\n');
fprintf('=== Plain LPC Vocoder Encoder ===\n');
fprintf('Using PARCOR coefficients (reflection coefficients)\n');
fprintf('Performing LPC analysis on all files...\n');
fprintf('========================================================\n');

for f = 1:4
    fprintf('\nAnalyzing File %d...\n', f);
    [aCoeff, resid, pitch, G, parcor, stream] = proclpc(speech_all{f}, sr, L, fr, fs, preemp);
    
    aCoeff_all{f} = aCoeff;
    resid_all{f} = resid;
    pitch_all{f} = pitch;
    G_all{f} = G;
    parcor_all{f} = parcor;
    stream_all{f} = stream;
    
    nframes = size(aCoeff, 2);
    fprintf('  Number of frames: %d\n', nframes);
end

%% ============== QUANTIZATION FOR ALL FILES ==============
fprintf('\n========================================================\n');
fprintf('--- Quantization (Scalar Quantization) ---\n');
fprintf('Quantizing all files...\n');
fprintf('========================================================\n');

parcor_quant_all = cell(4,1);
pitch_quant_all = cell(4,1);
G_quant_all = cell(4,1);
vuv_all = cell(4,1);
G_min_all = zeros(4,1);
G_max_all = zeros(4,1);

parcor_min = -1; parcor_max = 1;
parcor_levels = 2^parcor_bits;
pitch_max = 147;  % Maximum pitch period in samples at 8 kHz

for f = 1:4
    fprintf('\nQuantizing File %d...\n', f);
    
    % Quantize PARCOR coefficients (reflection coefficients)
    parcor_quant = round((parcor_all{f} - parcor_min) / (parcor_max - parcor_min) * (parcor_levels - 1));
    parcor_quant = min(max(parcor_quant, 0), parcor_levels - 1);
    parcor_quant_all{f} = parcor_quant;
    
    % Quantize pitch
    pitch_quant = round(pitch_all{f} / pitch_max * (2^pitch_bits - 1));
    pitch_quant = min(max(pitch_quant, 0), 2^pitch_bits - 1);
    pitch_quant_all{f} = pitch_quant;
    
    % Quantize gain using logarithmic quantization
    G_log = log(G_all{f} + eps);
    G_min = min(G_log); G_max = max(G_log);
    G_quant = round((G_log - G_min) / (G_max - G_min + eps) * (2^gain_bits - 1));
    G_quant_all{f} = G_quant;
    G_min_all(f) = G_min;
    G_max_all(f) = G_max;
    
    % Voiced/Unvoiced decision
    vuv = (pitch_all{f} > 0);
    vuv_all{f} = vuv;
    
    fprintf('  PARCOR: %d bits per coefficient × %d = %d bits\n', parcor_bits, L, L*parcor_bits);
    fprintf('  Pitch: %d bits\n', pitch_bits);
    fprintf('  Gain: %d bits (logarithmic)\n', gain_bits);
    fprintf('  V/UV: %d bit\n', vuv_bits);
    fprintf('  Total: %d bits per frame\n', total_bits_per_frame);
end

%% ============== PACK BITSTREAMS FOR ALL FILES ==============
fprintf('\n========================================================\n');
fprintf('--- Packing bitstreams ---\n');
fprintf('========================================================\n');

bitstream_files = cell(4,1);
for f = 1:4
    bitstream_file = sprintf('%s%d.bin', bitstream_prefix, f);
    bitstream_files{f} = bitstream_file;
    
    fprintf('\nPacking File %d to: %s\n', f, bitstream_file);
    fid = fopen(bitstream_file, 'w');
    
    nframes = size(aCoeff_all{f}, 2);
    
    % Write header
    fwrite(fid, nframes, 'uint32');
    fwrite(fid, sr, 'uint32');
    fwrite(fid, L, 'uint8');
    fwrite(fid, G_min_all(f), 'float32');
    fwrite(fid, G_max_all(f), 'float32');
    
    % Write quantized parameters frame by frame
    for i = 1:nframes
        fwrite(fid, parcor_quant_all{f}(:, i), 'uint8');
        fwrite(fid, pitch_quant_all{f}(i), 'uint8');
        fwrite(fid, G_quant_all{f}(i), 'uint8');
        fwrite(fid, vuv_all{f}(i), 'uint8');
    end
    fclose(fid);
end

%% ============== DECODER (SYNTHESIS) FOR ALL FILES ==============
fprintf('\n========================================================\n');
fprintf('=== Plain LPC Vocoder Decoder ===\n');
fprintf('Decoding all files...\n');
fprintf('========================================================\n');

parcor_dec_all = cell(4,1);
pitch_dec_all = cell(4,1);
G_dec_all = cell(4,1);
aCoeff_dec_all = cell(4,1);

for f = 1:4
    fprintf('\nDecoding File %d: %s\n', f, bitstream_files{f});
    
    % Read bitstream
    fid = fopen(bitstream_files{f}, 'r');
    nframes_dec = fread(fid, 1, 'uint32');
    sr_dec = fread(fid, 1, 'uint32');
    L_dec = fread(fid, 1, 'uint8');
    G_min_dec = fread(fid, 1, 'float32');
    G_max_dec = fread(fid, 1, 'float32');
    
    % Read quantized parameters
    parcor_quant_dec = zeros(L_dec, nframes_dec);
    pitch_quant_dec = zeros(1, nframes_dec);
    G_quant_dec = zeros(1, nframes_dec);
    vuv_dec = zeros(1, nframes_dec);
    
    for i = 1:nframes_dec
        parcor_quant_dec(:, i) = fread(fid, L_dec, 'uint8');
        pitch_quant_dec(i) = fread(fid, 1, 'uint8');
        G_quant_dec(i) = fread(fid, 1, 'uint8');
        vuv_dec(i) = fread(fid, 1, 'uint8');
    end
    fclose(fid);
    
    % Dequantize parameters
    parcor_dec = parcor_quant_dec / (parcor_levels - 1) * (parcor_max - parcor_min) + parcor_min;
    pitch_dec = pitch_quant_dec / (2^pitch_bits - 1) * pitch_max;
    G_log_dec = G_quant_dec / (2^gain_bits - 1) * (G_max_dec - G_min_dec) + G_min_dec;
    G_dec = exp(G_log_dec);
    
    parcor_dec_all{f} = parcor_dec;
    pitch_dec_all{f} = pitch_dec;
    G_dec_all{f} = G_dec;
    
    % Convert PARCOR back to LPC coefficients
    aCoeff_dec = zeros(L_dec + 1, nframes_dec);
    for i = 1:nframes_dec
        lpc_poly = rc2poly(parcor_dec(:, i));
        aCoeff_dec(:, i) = lpc_poly(:);
    end
    aCoeff_dec_all{f} = aCoeff_dec;
end

%% ============== SYNTHESIS FOR ALL FILES ==============
fprintf('\n========================================================\n');
fprintf('Synthesizing speech for all files...\n');
fprintf('========================================================\n');

synth_all = cell(4,1);
for f = 1:4
    fprintf('Synthesizing File %d...\n', f);
    synth_speech = synlpc(aCoeff_dec_all{f}, pitch_dec_all{f}, sr, G_dec_all{f}, fr, fs, preemp);
    synth_speech_norm = synth_speech / max(abs(synth_speech));
    synth_all{f} = synth_speech_norm;
end

%% ============== SAVE ALL OUTPUTS ==============
fprintf('\n========================================================\n');
fprintf('Saving all output files...\n');
fprintf('========================================================\n');

output_files = cell(4,1);
for f = 1:4
    output_file = sprintf('%soutput_%d.wav', output_prefix, f);
    audiowrite(output_file, synth_all{f}, sr);
    output_files{f} = output_file;
    fprintf('File %d saved: %s\n', f, output_file);
end

%% ---------------------- STOP TIMER -----------------------
total_runtime = toc;

%% ============== BIT-RATE CALCULATION ==============
fprintf('\n========================================================\n');
fprintf('================= BIT-RATE ESTIMATION ==================\n');
fprintf('========================================================\n');

bitrates_theoretical = zeros(4,1);
bitrates_actual = zeros(4,1);
durations = zeros(4,1);
seg_snr_values = zeros(4,1);

frame_rate = 1000 / fr;  % frames per second

for f = 1:4
    nframes = size(aCoeff_all{f}, 2);
    duration = length(speech_all{f}) / sr;
    durations(f) = duration;
    
    % Theoretical bit-rate
    theoretical_bitrate = total_bits_per_frame * frame_rate / 1000;  % kbps
    bitrates_theoretical(f) = theoretical_bitrate;
    
    % Actual bit-rate (including header)
    file_info = dir(bitstream_files{f});
    file_size_bits = file_info.bytes * 8;
    actual_bitrate = file_size_bits / duration / 1000;  % kbps
    bitrates_actual(f) = actual_bitrate;
    
    % Segmental SNR
    min_len = min(length(speech_all{f}), length(synth_all{f}));
    seg_snr_values(f) = segsnr(speech_all{f}(1:min_len), synth_all{f}(1:min_len), sr);
    
    fprintf('\nFile %d:\n', f);
    fprintf('  Frames: %d\n', nframes);
    fprintf('  Duration: %.2f seconds\n', duration);
    fprintf('  Bits per frame: %d (PARCOR:%d + Pitch:%d + Gain:%d + V/UV:%d)\n', ...
            total_bits_per_frame, L*parcor_bits, pitch_bits, gain_bits, vuv_bits);
    fprintf('  Theoretical bit-rate: %.2f kbps (%d bits/frame × %.0f frames/s)\n', ...
            theoretical_bitrate, total_bits_per_frame, frame_rate);
    fprintf('  Actual bit-rate: %.2f kbps (including header)\n', actual_bitrate);
    fprintf('  SEGSNR: %.2f dB\n', seg_snr_values(f));
end

avg_theoretical = mean(bitrates_theoretical);
avg_actual = mean(bitrates_actual);

fprintf('\n--- Average Bit-rates ---\n');
fprintf('  Average Theoretical Bit-rate: %.2f kbps\n', avg_theoretical);
fprintf('  Average Actual Bit-rate: %.2f kbps\n', avg_actual);

%% ============== RESULTS TABLE ==============
fprintf('\n========================================================\n');
fprintf('================= FINAL RESULTS ==================\n');
fprintf('========================================================\n');

T = table(input_files, output_files, seg_snr_values, bitrates_theoretical, bitrates_actual, durations, ...
    'VariableNames', {'Original_File','Output_File','SEGSNR_dB','Bitrate_Theoretical_kbps','Bitrate_Actual_kbps','Duration_sec'});

disp(T);

fprintf('\n================= RUNTIME ANALYSIS ==================\n');
fprintf('Total Runtime: %.3f seconds\n', total_runtime);
fprintf('Average time per file: %.3f seconds\n', total_runtime/4);

fprintf('\nProcessing complete for all 4 files.\n');