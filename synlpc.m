function synWave = synlpc(aCoeff,pitch,sr,G,fr,fs,preemp)
% USAGE: synWave = synlpc(aCoeff,pitch,sr,G,fr,fs,preemp);
%
% (original header kept unchanged)

if (nargin < 5), fr = 20; end;
if (nargin < 6), fs = 30; end;
if (nargin < 7), preemp = .9378; end;

msfs = round(sr*fs/1000); % framesize in samples
msfr = round(sr*fr/1000); % framerate in samples
msoverlap = msfs - msfr;
ramp = [0:1/(msoverlap-1):1]';
[L1 nframe] = size(aCoeff); % L1 = 1+number of LPC coeffs

for frameIndex=1:nframe
    A = aCoeff(:,frameIndex);
    if ( pitch(frameIndex) ~= 0 )
        t = 0 : 1/sr : fs*10^(-3);
        d = 0 : 1/pitch(frameIndex) : 1;
        residFrame = (pulstran(t, d, 'tripuls', 0.001))';
        residFrame = residFrame + 0.01*randn(msfs+1,1);
    else
        residFrame = [];
        for m = 1:msfs
            residFrame = [residFrame; randn];
        end
    end;

    synFrame = filter(G(frameIndex), A', residFrame);

    if(frameIndex==1)
        synWave = synFrame(1:msfr);
    else
        synWave = [synWave;
                   overlap+synFrame(1:msoverlap).*ramp;
                   synFrame(msoverlap+1:msfr)];
    end

    if(frameIndex==nframe)
        synWave = [synWave; synFrame(msfr+1:msfs)];
    else
        overlap = synFrame(msfr+1:msfs).*flipud(ramp);
    end
end;

synWave = filter(1, [1 -preemp], synWave);
