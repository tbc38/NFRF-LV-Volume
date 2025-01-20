function [R, T] = segmentECGV3(ecg, f)
% input arguments: : ecg, HP filtered ECG | tecg, time axis for ECG samples
% also extracting T wave

% Normalize
ecg = zscore(ecg);

%% R wave
% Wavelet transform with 'sym4' wavelet to level 5 (looks like an ECG)
wt = modwt(ecg, 'sym4', 5);
% Define reconstruction wiht only the higher frequencies (levels 4 and 5)
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
% Take inverse wavelet to get only higher frequencies
y = imodwt(wtrec,'sym4');
% Take squared magnitude to improve distinction between R-peaks and adjacent values
y = abs(y).^2;
% Find R waves
[Rpeaks,Rlocs] = findpeaks(y, 'MinPeakHeight', 2, 'MinPeakDistance', 0.20*f);

% Move peaks such that they align with the real R-peak, rather than just nearby (sometimes the
% magnitude of the q/s-wave is greater than R-wave after wavelet transform so this is where peak is marked)
R = zeros(size(Rlocs));
for i = 1:length(Rlocs)
    lowerIdx = max([1, Rlocs(i)-100]);
    upperIdx = min([length(ecg), Rlocs(i)+100]);
    [pks, pklocs] = findpeaks(ecg(lowerIdx:upperIdx));
    [highestpk, highestpkidx] = max(pks);
    R(i) = lowerIdx+pklocs(highestpkidx)-1;
end

R = R(2:end-1);

%% T wave
bnds = [0.2 0.7]; % percentage of RR interval in which to look for T wave
T = zeros(1, length(R)-1);
for i = 1:length(R)-1
    R_i = R(i); 
    R_next = R(i+1);
    R_R = R_next - R_i;
    lowerIdx = round(R_i + R_R*bnds(1));
    upperIdx = round(R_i + R_R*bnds(2));
    
    [~, l] = max(ecg(lowerIdx:upperIdx));
    T(i) = lowerIdx+l-1;
end

R = R(1:end-1);


end