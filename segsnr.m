function Sseg = segsnr(x, y, Fs)
% SEGSNR Segmental SNR between reference x and processed y
%   Sseg = segsnr(x,y,Fs)
%   x,y : column vectors (same length)
%   Fs  : sampling rate (Hz)

% Ensure column vectors
x = x(:);
y = y(:);
N = min(length(x), length(y));
x = x(1:N);
y = y(1:N);

% Frame parameters
frame_ms = 20;           % 20 ms frames
frame_len = round(Fs * frame_ms/1000);
if frame_len < 1
    error('Sampling rate too low for frame length.');
end
hop = floor(frame_len/2); % 50% overlap

eps_energy = 1e-12;      % threshold to treat frame as silent
min_snr = -10;           % lower clip (dB)
max_snr =  35;           % upper clip (dB)

idx = 1:hop:(N - frame_len + 1);
numFrames = numel(idx);
snr_frames = nan(numFrames,1);

for k = 1:numFrames
    i = idx(k);
    xf = x(i:i+frame_len-1);
    yf = y(i:i+frame_len-1);
    en = sum(xf.^2);
    if en <= eps_energy
        % ignore very low-energy frames
        continue;
    end
    noise_en = sum((xf - yf).^2);
    if noise_en <= eps_energy
        frame_snr = max_snr;
    else
        frame_snr = 10*log10(en / noise_en);
        % clip
        frame_snr = min(max(frame_snr, min_snr), max_snr);
    end
    snr_frames(k) = frame_snr;
end

% Remove NaNs (silent frames)
valid = ~isnan(snr_frames);
if ~any(valid)
    Sseg = -Inf;
else
    Sseg = mean(snr_frames(valid));
end
end
