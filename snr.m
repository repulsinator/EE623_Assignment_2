% ======================================================================
% SEGSNR computation for 2 male + 2 female samples
% For ONE vocoder (e.g., CELP 11 kbps)
% ======================================================================

clc; clear; close all;

% -------------------------------------------------------------
% 1. Specify file names (EDIT THESE ONLY)
% -------------------------------------------------------------

% --- Original files (2 male, 2 female) ---

orig_files = {
    '/MATLAB Drive/assignment_2/Speech_sample/male_voice_1_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/male_voice_2_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/female_voice_1_8k.wav'
    '/MATLAB Drive/assignment_2/Speech_sample/female_voice_2_8k.wav'
};
% --- Vocoder output files (only ONE vocoder) ---
v_files = {
    '/MATLAB Drive/assignment_2/Output_files/plain_lpc_male_voice_1_output.wav'
    '/MATLAB Drive/assignment_2/Output_files/plain_lpc_male_voice_2_output.wav'
    '/MATLAB Drive/assignment_2/Output_files/plain_lpc_female_voice_1_output.wav'
    '/MATLAB Drive/assignment_2/Output_files/plain_lpc_female_voice_2_output.wav'
    
};

% -------------------------------------------------------------
% 2. Compute SEGSNR for all files
% -------------------------------------------------------------
num_files = 4;
seg = zeros(num_files,1);

fprintf("\n================  SEGSNR RESULTS  ==================\n");

for i = 1:num_files

    % Load original
    [x, Fs] = audioread(orig_files{i});
    x = x(:);

    % Load synthesized vocoder output
    [y, Fs2] = audioread(v_files{i});
    y = y(:);

    % Resample if needed
    if Fs2 ~= Fs
        y = resample(y, Fs, Fs2);
    end

    % Match length
    L = min(length(x), length(y));
    x = x(1:L);
    y = y(1:L);

    % Compute SEGSNR
    seg(i) = segsnr(x, y, Fs);

    fprintf("File %d (%s): SEGSNR = %.3f dB\n", i, orig_files{i}, seg(i));
end

% -------------------------------------------------------------
% 3. Display table
% -------------------------------------------------------------
T = table(orig_files, v_files, seg, ...
    'VariableNames', {'Original','Vocoder','SEGSNR_dB'});

fprintf("\n\n================  SEGSNR TABLE  ==================\n");
disp(T);




