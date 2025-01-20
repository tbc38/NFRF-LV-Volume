function [SNR] = signalSNR(ncs, fSig, fNcs, fNoise)
% signalSNR: return the snr from the FFT at a given frequency range

if ~exist('fNoise', 'var')
    fNoise = [10 25];
end


[P, f] = ncsFFT(ncs, fNcs, 0, [0.5 20]);

idxF = f > fSig(1) & f < fSig(2);
idxNoise = f > fNoise(1) & f < fNoise(2);
% not the max in the range, the highest peak in the range

pk = findpeaks(P(idxF));

if isempty(pk); pk = 0; end

% SNR = max(pk)/median(P(idxNoise));
SNR = 10*log10(max(pk)/median(P(idxNoise)));

end

