function [idx_rem, idxFlag] = PP_removeOutliers_V2(ncsPP, tncsPP, ncs, tncs)
% this function takes in the ncs peak-to-peak of each waveform as well as their times, it also takes
% in the ncs signal and its time reference
% the function returns the PP and tPP removed, as well as the indeces of removal
% V2: returns a flag if more than 90% of beats are removed or if a certain amount
% are jumps

jumpThr = 0.1; remThr = 0.2;

% cutting orignal signal to the window of the ncsPP
idxT = tncs > tncsPP(1) & tncsPP(end) > tncs;
ncs = ncs(idxT); ncs = ncs-mean(ncs);
tncs = tncs(idxT);

% standard outlier removal
idx_outlier = isoutlier(ncsPP, 'movmedian', 10); % isoutlier(ncsPP);

% finding jumps on ncs signal and removing surrounding beats
% threshold: based on average peak-to-peak on input
% window: 50ms
meanPP = 5*mean(ncsPP(~idx_outlier));
jump_indexes = findSignalJumps(ncs, 50, meanPP);
idxJumps = false(size(tncsPP));
% turn idxJumps into a logical array for the peak-to-peaks
for i = 1:length(tncsPP)-1 % for every heartbeat
    idxT = tncs > tncsPP(i) & tncs < tncsPP(i+1); % define logical array on time axis where its true
    for j = 1:length(idxT)
        for k = 1:length(jump_indexes)
            if jump_indexes(k) == j && idxT(j) % if the jump index lies within this heatbeat
                idxJumps(i) = true; % set jump to true
            end
        end
        % if ismember(j, jump_indexes) && idxT(j)
        %     idxJumps(i) = true;
        % end
    end
end

% remove one beat each side of a jump
idxJumps_copy = idxJumps;
for i = 2:length(idxJumps)-1
    if idxJumps(i)
        idxJumps_copy(i-1) = true; idxJumps_copy(i+1) = true;
    end
end
idxJumps = idxJumps_copy;

% figure; plot(tncs, ncs-mean(ncs)); hold on; 
% plot(tncs(jump_indexes), ncs(jump_indexes), 'ro'); 
% plot(tncsPP(idx_outlier), mean(ncs)*ones(size(tncsPP(idx_outlier))), 'b*'); 
% plot(tncsPP(idxJumps), mean(ncs)*ones(size(tncsPP(idxJumps))), 'g*'); 
% plot(tncsPP, ncsPP);

idx_rem = idx_outlier | idxJumps;

if sum(idx_rem)/length(idx_rem) > remThr || sum(idxJumps)/length(idxJumps) > jumpThr
    % figure; plot(tncs, ncs-mean(ncs)); hold on; 
    % plot(tncs(jump_indexes), ncs(jump_indexes), 'ro'); 
    % plot(tncsPP(idx_outlier), mean(ncs)*ones(size(tncsPP(idx_outlier))), 'b*'); 
    % plot(tncsPP(idxJumps), mean(ncs)*ones(size(tncsPP(idxJumps))), 'g*'); 
    % plot(tncsPP, ncsPP);
    % title([num2str(sum(idxJumps)) ' ' num2str(sum(idxJumps)/length(idxJumps))])
    disp('')
end

% idx_outlier_test = isoutlier(ncsPP, 'movmedian', 10);
% if sum(idx_outlier) > sum(idx_outlier_test)+3
%     figure; plot(tncs, ncs-mean(ncs)); hold on; 
%     plot(tncs(jump_indexes), ncs(jump_indexes), 'ro'); 
%     plot(tncsPP(idx_outlier), mean(ncs)*ones(size(tncsPP(idx_outlier))), 'b*'); 
%     plot(tncsPP(idxJumps), mean(ncs)*ones(size(tncsPP(idxJumps))), 'g*'); 
%     plot(tncsPP, ncsPP);
%     title([num2str(sum(idx_outlier)) ' ' num2str(sum(idx_outlier_test))])
%     disp('')
% end

% Quality checks:
% If all beats removed, then flag and make only one of them false so one
% beat stays present for later processing
if sum(idx_rem) == length(idx_rem)
    idx_rem(1:2) = false;
    idxFlag = true;
    return
    % flag if over x% removed or over y% are jumps
elseif sum(idx_rem)/length(idx_rem) > remThr || sum(idxJumps)/length(idxJumps) > jumpThr
    idxFlag = true;
    return
else
    idxFlag = false;
    return
end

end
