function [avgWaveform, allWaveforms, times] = ...
    segmentNCSWaveformsV2(ncs, tNcs, tR, segLen, ncsArtifacts, mean_win)

% function takes in the ncs and R-waves and returns the average waveform and each individual cut
% waveform. It removes those in the window of ncsArtifact. It also will combine the all-waveforms
% into a mean-filtered window fi desired

% segLen = 500;
allWaveforms = zeros(segLen, length(tR)-1);
allWaveformsNorm = zeros(segLen, length(tR)-1);

% Cut each waveform at R-wave
for i = 1:length(tR)-1
    tR_0 = tR(i); tR_1 = tR(i+1);
    idxNcs = tNcs > tR_0 & tNcs < tR_1; [~,loc] = max(idxNcs); idxStart(i) = loc;
    ncsCut = ncs(idxNcs);
    ncsCutInterp = interp1(1:length(ncsCut), ncsCut, linspace(1, length(ncsCut), segLen));
    allWaveforms(:, i) = ncsCutInterp;
    try
        allWaveformsNorm(:, i) = rescale(ncsCutInterp, -1, 1);
    catch
        allWaveformsNorm(:, i) = ncsCutInterp;
    end
    times(i) = tR_0;
end

% remove waveform at marked points
idxRem = true(size(idxStart));
for i = 1:length(idxStart)-1
    if sum(idxStart(i) < ncsArtifacts & idxStart(i+1) > ncsArtifacts)
        idxRem(i) = false;
    end
end

allWaveforms = allWaveforms(:, idxRem);
allWaveformsNorm = allWaveformsNorm(:, idxRem);
times = times(idxRem);
avgWaveform = mean(allWaveformsNorm, 2);

if mean_win ~= 1
    allWaveformsWin = zeros(size(allWaveforms, 1), size(allWaveforms, 2)-mean_win+1);
    for i = 1:size(allWaveformsWin, 2)
        allWaveformsWin(:,i) = mean(allWaveforms(:,i:i+mean_win-1), 2);
    end
    allWaveforms = allWaveformsWin;
    times = times(1:end-mean_win+1);
end

end