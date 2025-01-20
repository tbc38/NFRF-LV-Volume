function idx = ECGArrhythmiaCheck(bio, varargin)

thr = 0.20;

% Check if a second argument is provided
if nargin > 1
    thr = varargin{1};
end

% two checks: odd timing that doesn't recur, and odd heights that don't match surrounding
% given as input is the ecg waveform array and their times. 
tEcg = bio.t_allECGWaveforms;
idx = true(size(tEcg)); idx(end) = false;

% 5% RR interval difference to median beat in 'win' around it - as long as not mostly ARRs, should
% work?
win = 10;
RR = diff(tEcg); tRR = tEcg(1:end-1);
% sdthr = 2;
for i = 1:length(RR)
    
    RRtest = RR(i);
    
    lowerbnd = max([1 i-win]); upperbnd = min([i+win length(RR)]);
% 5% RR interval difference to any beat in 'win' around it - doesn't work cause ARRs come in bunches
%     highestRR_win = max([RR(lowerbnd:i-1) RR(i+1:upperbnd)]);
%     lowestRR_win = min([RR(lowerbnd:i-1) RR(i+1:upperbnd)]);
%     if RRtest > (1+thr)*highestRR_win || RRtest < (1-thr)*lowestRR_win
%         idx(i) = false;
%     end
% 3*SD across a large window
%     sd = std([RR(lowerbnd:i-1) RR(i+1:upperbnd)]);
%     med = median([RR(lowerbnd:i-1) RR(i+1:upperbnd)]);
%     if RRtest > med+sdthr*sd || RRtest < med-sdthr*sd
%         idx(i) = false;
%     end
% 15% RR interval difference to median beat in 'win' around it - as long as not mostly ARRs
% stepwise decreasing thr so very extreme beats won't impact median much
    for thr_ = 3*thr:-thr:thr
        RR_win = [RR(lowerbnd:i-1) RR(i+1:upperbnd)];
        idx_win = [idx(lowerbnd:i-1) idx(i+1:upperbnd)];
        med = median(RR_win(idx_win));
        if RRtest > (1+thr_)*med || RRtest < (1-thr_)*med
            idx(i) = false;
        end
    end
end

% figure; 
% ax1 = subplot(211);
% stem(tRR, RR); hold on; stem(tRR(idx), RR(idx));
% 
% ax2 = subplot(212); 
% plot(bio.tBio, bio.ecgFilt, 'k'); hold on;
% plot(bio.tBio(bio.RWave), bio.ecgFilt(bio.RWave), 'ro');
% R = bio.RWave(1:end-1);
% plot(bio.tBio(R(~idx)), bio.ecgFilt(R(~idx)), 'b^');
% linkaxes([ax1 ax2], 'x');
end