function [powerlab_sync] = alignPowerlabWithBIOV2_6(powerlab, powerlab_ch, bio_sig_orig, bio_date, showplot)
% V2: Syncrhonize a given BIO trace with a given powerlab trace

% Algorithm: Trying cross correlation? Remove outliers, find cross correlation for rough alignment
% after normalization, then tune? Will probably need to test amongts different clock offsets too:
% look at the p for the previous V1 method. Potential pitfall: the xcorr will be influenced by the
% baseline more heavily than by the actual waveforms. Can highpass filter, or implement a
% normalization on every subsection of the xcorr. Add time scaling (according to prev matched
% points) and account for record skipping
% Tried to split the powerlab recordings and indivdually find the best lags: 
%    this didn't work because when the shifted copy isn't totally over the signal then it will always have lower xcorr.
% Might be better to find the overall alignment, then if there's a break in it, address the
% break afterwards
% Currenlty addressing the break by using the timestamp gap between the two records and inputing
% zeros during the gap

% Shortcomings: 
% - if there are 2+ gaps in the same file, not yet addressed (trying to address in V2.1)
% V2.1: if there are 2+ gaps, it assumes a clock offset and that the gap starts at the beginning (eg
% the flag is 'beforeBreak'.) Keep an eye on if this fails at some point.
% V2.2: using the times from the powerlab and bio to guess a rough window where the signals align:
% should improve when the xcorr method fails on its own. the method of choosing the first powerlab
% signal to mark the time is imperfect, however it works for now. Same shortcoming from V2.1 about
% pre-post clock offset breaks exists here too.
% V2.3 edit: split the signals into pre/post break, find optimal correlation on the bio signal
% independently, then fill the necessary gap between with zeros. Check if the signal filled gap
% and the expected gap are similar, error/warn if they aren't.
% V2.4 edit: update one-break sync to ensure the pre-break/post-break signal combine to the proper
% length
% V2.5: improve the 2+ gap method; this method is unreliable often, so will make a window-wise sync
% similar to V2.3
% V2.6: improve the 2 gap method: rather than sticking zeros in the gap,
% stick as much of the signal as is aligned. Right now only applies to the
% middle section into the last section, but can propagate further if
% necessary

% extract shared powerlab signal
powerlab_sig = powerlab.data{powerlab_ch};
bio_sig = bio_sig_orig;

% if a certain channel is shorter than the rest of the plab channels (likely due to error), append zeros to
% make it long enough, give a warning, and put the warning in the powerlab.channel_labels struct
% v2-3: put zeros in front to see if that addresses PA/RVP issues in 2-7
for i = 1:length(powerlab.data)
    if length(powerlab.data{i}) < length(powerlab_sig)
        powerlab.data{i} = [zeros(length(powerlab_sig) - length(powerlab.data{i}), 1); powerlab.data{i}];
        powerlab.channel_labels(i).name = ['CHANNEL WARNING: ' powerlab.channel_labels(i).name];
        disp(['Number of samples in powerlab channel ' num2str(i) ', ' powerlab.channel_labels(i).name ', less than expected']);
    end
end

% remove outliers
outliers = isoutlier(powerlab_sig, 'movmedian', 100000, 'ThresholdFactor', 8);
powerlab_sig(outliers) = median(powerlab_sig);

outliers = isoutlier(bio_sig_orig, 'movmedian', 100000, 'ThresholdFactor', 8);
bio_sig(outliers) = median(bio_sig_orig);

% added highpass filter to remove baseline trending affecting xcorr
powerlab_sig = zscore(highpass(powerlab_sig, 1, 1000));
bio_sig = zscore(highpass(bio_sig, 1, 1000));

% bio_sig = bio_sig(1000:end-1000);

% extract effective start time of bio and powerlab
bio_time = bio_date.filetime; 
% need to skip some early, short recordings that don't impact the start time
for i = 1:length(powerlab.record_info)-1
    % Attempt 1: mark time as first record that takes up x% of distance to next record
    % if this doesn't work, can add up the 'missing' time from the start and define the search
    % window as that. This would have downside of drastically increasing the size of 05-12 where plab
    % date starts the day before for some reason
    
    d_samps = 1000*seconds(powerlab.record_info(i+1).data_start_dt-powerlab.record_info(i).data_start_dt);
    perc_thr = 0.9;
%     disp(d_samps); disp(powerlab.record_info(i).n_ticks); disp(d_samps*perc_thr);
    if powerlab.record_info(i).n_ticks > d_samps*perc_thr
        plab_time = powerlab.record_info(i).data_start_dt;
        disp(['Using ' num2str(i) ' powerlab record']);
        break
    end
end
seconds_offset = seconds(bio_time - plab_time);

% do cross correlation across whole trace while tuning to different clock offsets
clockOffsetBnds = 1.0:0.00005:1.0001; k = 1; % typical: 1.0001
clkCorrs = zeros(size(clockOffsetBnds)); clkIncs = zeros(size(clockOffsetBnds));
for clk = clockOffsetBnds
    % define a new axis of indeces scaled by the slight offset
    idx_foo = (1:length(powerlab_sig))*clk;
    % interpolate a new powerlab signal onto the original indeces
    foo_interp = interp1(idx_foo, powerlab_sig, 1:length(powerlab_sig), 'linear', 0);
       
    % with new interpolated plab signal, find peak correlation & index. Store.
    [r, lags] = xcorr(foo_interp', bio_sig);
    
    % V2.2 edit: keep only correlations in a window of expected lags, based on datetimes
    % allow a 50 min window where alignment can be found
    samples_offset = seconds_offset * 1000 * clk; search_width = 1000 * clk * 60 * 50; 
    expectedLags = [samples_offset - search_width; samples_offset + search_width];
    
    idxLags = lags > expectedLags(1) & lags < expectedLags(2); 
    lags = lags(idxLags); 
    r = r(idxLags);
    
    [clkCorrs(k), l] = max(r); 
    clkIncs(k) = lags(l);
    k = k+1;
end
[~, l] = max(clkCorrs); clkOffset = clockOffsetBnds(l); inc = clkIncs(l);
idx_foo = (1:length(powerlab_sig))*clkOffset;
powerlab_sig_interp = interp1(idx_foo, powerlab_sig, 1:length(powerlab_sig), 'linear', 0)';
disp(['Syncing: clock offset determined to be ' num2str(clkOffset)]);

% if signal is too short (ie shorter than comparison signal), pad some zeros to the end
if length(powerlab_sig_interp) < inc+length(bio_sig)
    dif = inc+length(bio_sig) - length(powerlab_sig_interp);
    powerlab_sig_interp = [powerlab_sig_interp; zeros(2*dif, 1)];
    for i = 1:length(powerlab.data)
        powerlab.data{i} = [powerlab.data{i}; zeros(2*dif, 1)];
    end 
end

% Check if there's a break in the signal, and address
bnds = [inc inc+length(bio_sig)];
a = cumsum(powerlab.channel_labels(1).n_samples(1:end-1))*clkOffset;
if sum(a > bnds(1) & a < bnds(2)) == 2
    % If there are 2+ breaks over optimal xcorr, redo the process with a hardcoded clock offset then
    % tune the gaps accordingly
    warning('2 record gaps in powerlab signal contained within BIO trace. Sync might be less reliable.');
    
    % clock offset because no longer reliable. Hardcoded to 1.0001
    clkOffset = 1.0001; inc = clkIncs(clockOffsetBnds == 1.0001);
    idx_foo = (1:length(powerlab_sig))*clkOffset;
    powerlab_sig_interp = interp1(idx_foo, powerlab_sig, 1:length(powerlab_sig), 'linear', 0)';
    disp(['Using hardcoded clock offset of ' num2str(clkOffset)]);
    
    % if signal is too short (ie shorter than comparison signal), pad some zeros to the end
    if length(powerlab_sig_interp) < inc+length(bio_sig)
        dif = inc+length(bio_sig) - length(powerlab_sig_interp);
        powerlab_sig_interp = [powerlab_sig_interp; zeros(2*dif, 1)];
        for i = 1:length(powerlab.data)
            powerlab.data{i} = [powerlab.data{i}; zeros(2*dif, 1)];
        end 
    end
    
    % Check where the breaks are in the interpolated signal
    bnds = [inc inc+length(bio_sig)];
    a = cumsum(powerlab.channel_labels(1).n_samples(1:end-1))*clkOffset;
    breaks = (a > bnds(1) & a < bnds(2));
    
    % V2.5 edit: find and tune breaks in signal: based on V2.3: record-wise search and correlation
    % find index of breaks (should be of size 2)
    idxBreak = round(a(breaks));

    % split the signals into pre/mid/post break, find optimal correlation on the bio signal
    % independently, then fill the necessary gap between with zeros.
    segPreBreak = powerlab_sig_interp(inc:idxBreak(1)-1); bioPreBreak = bio_sig(1:idxBreak(1)-inc-1);
    [r, lags] = xcorr(segPreBreak, bioPreBreak);
    [~, k] = max(r); optLag_PreBreak = lags(k);
    
    segMidBreak = powerlab_sig_interp(idxBreak(1):idxBreak(2)); bioMidBreak = bio_sig(idxBreak(1)-inc:idxBreak(2)-inc-1);
    [r, lags] = xcorr(segMidBreak, bioMidBreak);
    [~, k] = max(r); optLag_MidBreak = lags(k);

    segPostBreak = powerlab_sig_interp(idxBreak(2):inc+length(bio_sig)); bioPostBreak = bio_sig(idxBreak(2)-inc:end);
    [r, lags] = xcorr(segPostBreak, bioPostBreak);
    [~, k] = max(r); optLag_PostBreak = lags(k);

    % Not yet checking with timestamps as well; maybe we should?

    % apply the break and offset to all powerlab traces
    for i = 1:length(powerlab.data)
        idx_foo = (1:length(powerlab.data{i}))*clkOffset;
        sig_interp = interp1(idx_foo, powerlab.data{i}, 1:length(powerlab.data{i}), 'linear', 0)';
        % segPreBreak
        if optLag_PreBreak > 0
            foo = [sig_interp(inc+optLag_PreBreak:idxBreak(1)-1); zeros(optLag_PreBreak, 1)];
            foo = foo(1:length(bioPreBreak));
        else
            foo = [zeros(abs(optLag_PreBreak), 1); sig_interp(inc:idxBreak(1)-1)];
            foo = foo(1:length(bioPreBreak));
        end
        % segMidBreak
        if optLag_MidBreak > 0
            bar = [sig_interp(idxBreak(1)+optLag_MidBreak:idxBreak(2)-1); zeros(optLag_MidBreak, 1)];
            bar = bar(1:length(bioMidBreak));
        else
            bar = [zeros(abs(optLag_MidBreak), 1); sig_interp(idxBreak(1):idxBreak(2)-1)];
            % v2.6 special case for 2/5 since I cut off so much VCO, will stick in the post break signal
            cutMidBreakSig = bar(length(bioMidBreak):end); 
            bar = bar(1:length(bioMidBreak)); 
        end
        % segPostBreak
        if optLag_PostBreak > 0
            baz = [sig_interp(idxBreak(2)+optLag_PostBreak:inc+length(bio_sig)); zeros(optLag_PostBreak, 1)];
            baz = baz(1:length(bioPostBreak));
        else
            % V2.6 only use this case if the mid break signal was moved forward adjacent to this signal
            if optLag_MidBreak < 0 && length(cutMidBreakSig) < abs(optLag_PostBreak)
                % in this case, put the cut off signal in front of the
                % zeros
                baz = [cutMidBreakSig; zeros(abs(optLag_PostBreak)-length(cutMidBreakSig), 1); sig_interp(idxBreak(2):inc+length(bio_sig))];
                baz = baz(1:length(bioPostBreak));
            else
                baz = [zeros(abs(optLag_PostBreak), 1); sig_interp(idxBreak(2):inc+length(bio_sig))];
                baz = baz(1:length(bioPostBreak));
            end
        end
        powerlab_sync.data{i} = [foo; bar; baz];
    end

elseif sum(a > bnds(1) & a < bnds(2)) == 1
    % If there's only one break, impute the zeros and apply to all powerlab traces
    disp('One break in signal');
    % find index of break
    [~, l] = max(a > bnds(1) & a < bnds(2));
    idxBreak = round(a(l));
    
    % V2.3 edit: split the signals into pre/post break, find optimal correlation on the bio signal
    % independently, then fill the necessary gap between with zeros. Check if the signal filled gap
    % and the expected gap are similar, error/warn if they aren't.
    segPreBreak = powerlab_sig_interp(inc:idxBreak-1); bioPreBreak = bio_sig(1:idxBreak-inc);
    [r, lags] = xcorr(segPreBreak, bioPreBreak);
    [~, k] = max(r); optLag_PreBreak = lags(k);
    segPreBreak_tuned = powerlab_sig_interp(inc+optLag_PreBreak:idxBreak-1);
    
    segPostBreak = powerlab_sig_interp(idxBreak:inc+length(bio_sig)); bioPostBreak = bio_sig(idxBreak-inc+1:end);
    [r, lags] = xcorr(segPostBreak, bioPostBreak);
    [~, k] = max(r); optLag_PostBreak = lags(k);
    segPostBreak_tuned = powerlab_sig_interp(idxBreak:inc+length(bio_sig)+optLag_PostBreak);
    
    expectedGap = length(bio_sig)-length(segPreBreak_tuned)-length(segPostBreak_tuned);
    
    % find expected time delta between them based on timestamps
    startLag = seconds(powerlab.record_info(l+1).data_start_dt - powerlab.record_info(l).data_start_dt)/powerlab.record_info(l).tick_dt;
    timestampGap = (startLag - powerlab.record_info(l).n_ticks)*clkOffset;
    
    % V2.4 edit: if sync fails and the resulting waveform is wrong (doesn't match expected gap),
    % just use the timestamp gap and xcorr naively
    if abs(expectedGap - timestampGap) > 1000
        warning('One-break synchronization exceeds tolerance to expected gap. Using timestamp gap naively.')
        % 10000 sample heuristic window applied to widen search space
        foo = [powerlab_sig_interp(inc-10000:idxBreak-1);...
            zeros(round(timestampGap), 1);...
            powerlab_sig_interp(idxBreak:inc+length(bio_sig)+10000)];
        [r, lags] = xcorr(foo, bio_sig);
        [~, k] = max(r); optLag = lags(k);
        for i = 1:length(powerlab.data)
            idx_foo = (1:length(powerlab.data{i}))*clkOffset;
            sig_interp = interp1(idx_foo, powerlab.data{i}, 1:length(powerlab.data{i}), 'linear', 0)';
            foo = [sig_interp(inc-10000:idxBreak-1);...
                zeros(round(timestampGap), 1);...
                sig_interp(idxBreak:inc+length(bio_sig)+10000)];
            foo = foo(optLag:optLag+length(bio_sig)-1);
            powerlab_sync.data{i} = foo;
        end
        
    else % otherwise: working fine just apply to all
        % apply the break and offset to all powerlab traces
        for i = 1:length(powerlab.data)
            idx_foo = (1:length(powerlab.data{i}))*clkOffset;
            sig_interp = interp1(idx_foo, powerlab.data{i}, 1:length(powerlab.data{i}), 'linear', 0)';
            foo = [sig_interp(inc+optLag_PreBreak:idxBreak-1);...
                zeros(round(expectedGap), 1);...
                sig_interp(idxBreak:inc+length(bio_sig)+optLag_PostBreak)];
            powerlab_sync.data{i} = foo;
        end
    end
elseif sum(a > bnds(1) & a < bnds(2)) > 2
    error('3+ breaks in signal: sync too complicated');
else
    % No breaks: just impute the offset on all powerlab traces
    disp('No breaks in signal');
    for i = 1:length(powerlab.data)
        idx_foo = (1:length(powerlab.data{i}))*clkOffset;
        sig_interp = interp1(idx_foo, powerlab.data{i}, 1:length(powerlab.data{i}), 'linear', 0)';
        powerlab_sync.data{i} = sig_interp(inc:inc+length(bio_sig)-1);
    end
end

% convert to double
for i = 1:length(powerlab.data)
    powerlab_sync.data{i} = double(powerlab_sync.data{i}');
end

powerlab_sync.channel_labels = powerlab.channel_labels;

if showplot
    figure('Position', [50 100 600 300]);
    plot(zscore(powerlab_sync.data{powerlab_ch})); hold on; plot(zscore(bio_sig_orig));
    legend('PLab', 'BioPac');
end

end