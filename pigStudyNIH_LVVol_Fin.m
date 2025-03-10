%% pigStudyNIH_EDV_Fin.m

filePath = ['']; % insert file path here
inc = 1;
fileDate = '01-26-24'; routines = {'dexmed'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 0; 
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1; 
    clearvars tBnds
end

fileDate = '02-05-24'; routines = {'dexmed', 'phenylephrine-2'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 1; 
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1; 
    clearvars tBnds
end

fileDate = '02-06-24'; routines = {'dexmed'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 2;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end
 
fileDate = '03-04-24'; routines = {'dexmed', 'phenylephrine'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 5;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '03-05-24'; routines = {'dexmed', 'phenylephrine'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 6;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '03-06-24'; routines = {'phenylephrine'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 7; 
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '03-07-24'; routines = {'dexmed', 'phenylephrine2'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 8;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

clearvars -except all*

for i = 1:length(all_dates)
    dFileDate = all_dates{i}; interventionType = all_ints{i};
    switch dFileDate
        case '01-26-24'
            switch interventionType
                case 'dexmed'
                    all_bnds{i} = {[40806 40821]};
                case 'venacavaocclusion2' % note: no response at these bnds
                    all_bnds{i} = {[41544 41604]};
                case 'venacavaocclusion3'
                    all_bnds{i} = {[46640 46741]};
            end
        case '02-05-24'
            switch interventionType
                case 'atropine'
                    all_bnds{i} = {[50834 50934]};
                case 'dexmed'
                    all_bnds{i} = {[51161 51179], [51226 51243]};
                case 'dobutamine'
                    all_bnds{i} = {[50338 50390]}; 
                case 'phenylephrine-2'
                    all_bnds{i} = {[49812 49925]};
                case 'venacavaocclusion'
                    all_bnds{i} = {[52400 52452]};
            end
        case '02-06-24'
            switch interventionType
                case 'atropine'
                    all_bnds{i} = {[52500 52536]};
                case 'dexmed'
                    all_bnds{i} = {[52853 52945]};
                case 'dobutamine'
                    all_bnds{i} = {[52007 52067]};
                case 'phenylephrine2'
                    all_bnds{i} = {[51324 51380]};
                case 'venacavaocclusion'
                    all_bnds{i} = {[53762 53809]};
            end
        case '02-07-24'
            switch interventionType
                case 'atropine'
                    all_bnds{i} = {[54920 54979]};
                case 'venacavaocclusion'
                    all_bnds{i} = {[55224 55272]};
            end
        case '02-08-24'
            switch interventionType
                case 'atropine'
                    all_bnds{i} = {[53426 53481], [53503 53537]};
                case 'dexmed'
                    all_bnds{i} = {[53880 53902] [53935 53965]};
                case 'dobutamine'
                    all_bnds{i} = {[52995 53073]};
                case 'venacavaocclusion'
                    all_bnds{i} = {[54624 54663]};
            end
        case '03-04-24'
            switch interventionType
                case 'dexmed'
                    all_bnds{i} = {[57035 57048], [57060 57071]};
                case 'phenylephrine'
                    all_bnds{i} = {[60157 60205]};
                case 'venacavaocclusion4'
                    all_bnds{i} = {[60340 60370]};
            end
        case '03-05-24'
            switch interventionType
                case 'dexmed'
                    all_bnds{i} = {[47300 47364]};
                case 'dobutamine'
                    all_bnds{i} = {[49026 49082]};
                case 'phenylephrine'
                    all_bnds{i} = {[48389 48454]};
                case 'venacavaocclusion'
                    all_bnds{i} = {[50210 50238]};
                case 'venacavaocclusion2'
                    all_bnds{i} = {[50630 50661]};
            end
        case '03-06-24'
            switch interventionType
                case 'dobutamine'
                    all_bnds{i} = {[47386 47457]};
                case 'dobutamine2'
                    all_bnds{i} = {[47886 47949]};
                case 'phenylephrine'
                    all_bnds{i} = {[46184 46252]};
                case 'venacavaocclusion2' % pre occlude
                    all_bnds{i} = {[49335 49377]};
            end
        case '03-07-24'
            switch interventionType
                case 'dexmed'
                    all_bnds{i} = {[44701 44773]};
                case 'dobutamine'
                    all_bnds{i} = {[44348 44418]};
                case 'phenylephrine2'
                    all_bnds{i} = {[43790 43860]};
                case 'cavalocclusion' % pre occlude
                    all_bnds{i} = {[45845 45890]};
            end
            
    end
end
clearvars -except all*

% Analyze each routine
all_corrs = zeros(length(all_ncs{1}.chSel), length(all_ncs));
all_weightedCorrs = zeros(1, length(all_ncs)); 
all_ncs_weightedavg = cell(1, length(all_ncs));
all_LVCath = cell(1, length(all_ncs));
all_ecg = cell(1, length(all_ncs));
all_tncs = cell(1, length(all_ncs));
all_ncs_proc = cell(1, length(all_ncs));
all_snr = cell(1, length(all_ncs));
all_lin = cell(1, length(all_ncs));
for i = 1:length(all_bnds)
    dFileDate = all_dates{i}; 
    interventionType = all_ints{i};
    ch_sel = all_ncs{i}.chSel;
    ncs = zscore(all_ncs{i}.ncsBPFilt); 
    ecg = zscore(all_bio{i}.ecgNorm);
    
    % S1 routines
    if (contains(dFileDate, '03-05') && contains(interventionType, 'ph')) ||...
        (contains(dFileDate, '03-07') && contains(interventionType, 'ph')) ||...
        (contains(dFileDate, '03-05') && contains(interventionType, 'dex')) ||...
        (contains(dFileDate, '03-07') && contains(interventionType, 'dex'))  
        LVCath = zscore(all_bio{i}.S1_vol_BPfilt);
    else
        LVCath = zscore(all_bio{i}.LV_vol_BPfilt);
    end
    
    tncs = all_ncs{i}.tNcsDS; 
    tbio = all_bio{i}.tBio;
    bnds = all_bnds{i};

    idxBnds = all_bio{i}.tR > bnds{1}(1) & all_bio{i}.tR < bnds{end}(2);
    tR = all_bio{i}.tR(idxBnds);
    
    % extract waveforms within bnds
    idxBndsWave = logical(all_ncs{i}.t_allWaveforms);
    inc = 1;
    for t = all_ncs{i}.t_allWaveforms
        flag = false;
        for j = 1:length(bnds)
            if t > bnds{j}(1) && t < bnds{j}(2)
                flag = true;
            end
        end
        idxBndsWave(inc) = flag;
        inc = inc+1;
    end
    t_waveforms = all_ncs{i}.t_allWaveforms(idxBndsWave);

    % Normalize polarity based on HB systole direction
    HB_polarity = zeros(1, length(ch_sel));
    for ch = all_ncs{i}.chSel
        waveforms = all_ncs{i}.allBPWaveforms{ch}(:, idxBndsWave);
        Twaves = all_bio{i}.allidxTWaves(idxBndsWave);
        polarity = zeros(size(waveforms, 2), 1);
        for j = 1:length(polarity)
            polarity(j) = waveformPolarity(waveforms(:,j), ...
                Twaves(j), 1, 'mode');
        end
        HB_polarity(ch) = mode(polarity);
    end
    ncs = ncs(:, ch_sel);
    HB_polarity = HB_polarity(:, ch_sel);
    ncs = -ncs.*HB_polarity; % Normalize polarity

    % Cut signals to bnds: removing parts that bnds skip and stapleing
    % ends together    
    idxBnds = logical(tncs);
    for j = 1:length(bnds)
        idxBnds(j, :) = tncs >= bnds{j}(1) & tncs <= bnds{j}(2);
    end
    idxBnds = any(idxBnds, 1);

    ncs = ncs(idxBnds, :);
    tncs = tncs(idxBnds);
    
    % Interpolate ground truth signals to tncs
    LVCath = interp1(tbio, LVCath, tncs, 'linear', 'extrap')';
    ecg = interp1(tbio, ecg, tncs, 'linear', 'extrap')';
    
    % calculate snr/linearity of NCS signals
    RR = diff(t_waveforms); RR = RR(RR < 5); % filter out big gaps
    idx_outlier = isoutlier(RR); RR = RR(~idx_outlier);
    RRbnds = [0.9*min(RR) 1.05*max(RR)];
    ncs_snr = []; ncs_lin = [];
    f = all_ncs{i}.fNcsDS; % Sample rate
    for ch = 1:size(ncs, 2)
        try
            ncs_snr(ch) = signalSNR(ncs(:, ch), sort(1./RRbnds), f, [3 10]);
            ncs_lin(ch) = signalLinearity(ncs(:, ch), sort(1./RRbnds), f);
        catch
            ncs_snr(ch) = 0; ncs_lin(ch) = 0;
        end
    end
    
    % Calculate correlation metrics
    corrVals = zeros(size(ncs, 2), 1);
    for ch = 1:size(ncs, 2)
        corrVals(ch) = corr(ncs(:, ch), LVCath);
    end

    % Store correlations & signals
    all_corrs(:, i) = corrVals; 
    all_ncs_proc{i} = ncs;
    all_LVCath{i} = LVCath;
    all_ecg{i} = ecg;
    all_tncs{i} = tncs;
    all_snr{i} = ncs_snr;
    all_lin{i} = ncs_lin;

end
clearvars -except all*   

% Use CV to find optimal parameters for quality-weighted channel combining (with channel selection)
interventions = 1:length(all_ncs); n = interventions(end)-1;
CV_combinations = nchoosek(interventions, n);
channel_selections = 1:8; % Example: channels 1 through 7

% define parameter search
parameter_range = 0:0.2:1;
[x1, x2, x3, x4] = ndgrid(parameter_range, parameter_range, parameter_range, parameter_range);
grid_search = [x1(:), x2(:), x3(:), x4(:)]';
idx_grid = logical(1:size(grid_search, 2));
for j = 1:size(grid_search, 2)
    if grid_search(1, j) >= grid_search(2, j) || grid_search(3, j) >= grid_search(4, j)
        idx_grid(j) = false;
    end
end
grid_search = grid_search(:, idx_grid);
optimal_parameter = [];
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    mean_corr_CV = [];
    for j = 1:size(grid_search, 2) % search over possible hyperparameters
        corr_CV = [];
        for k = 1:size(CV_combinations, 2) % loop over each routine in CV combinations with a set of hyperparameters 
            ncs_lin = all_lin{CV_combinations(i, k)}(channel_selections); % Apply channel selection
            ncs_snr = all_snr{CV_combinations(i, k)}(channel_selections); % Apply channel selection
            ncs_proc = all_ncs_proc{CV_combinations(i, k)}(:, channel_selections); % Apply channel selection
            LVCath = all_LVCath{CV_combinations(i, k)};
            
            ncs_lin_rescale = rescale(ncs_lin, grid_search(1, j), grid_search(2, j));
            ncs_snr_rescale = rescale(ncs_snr, grid_search(3, j), grid_search(4, j));
            ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;
            ncs_proc_weight = zeros(length(ncs_proc(:,1)), length(ncs_snr)); 
            for ch = 1:length(ncs_snr) 
                ncs_proc_weight(:, ch) = ncs_proc(:,ch) .* ncs_linsnr_weight(ch);
            end
            ncs_proc_weight_mean = sum(ncs_proc_weight, 2);
            ncs_proc_weight_mean = zscore(ncs_proc_weight_mean / sum(ncs_linsnr_weight));
            % optimizing for pearson's correlation coeff
            corr_CV(k) = corr(ncs_proc_weight_mean, LVCath);
        end
        mean_corr_CV(j) = mean(corr_CV); 
    end
    [m, l] = max(mean_corr_CV);
    optimal_parameter(:, i) = grid_search(:, l);
end
clearvars -except all* data* optimal_parameter CV_combinations interventions channel_selections

% evaluate
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    % Test LOO correlation with optimal parameters
    LOOSub = setdiff(interventions, CV_combinations(i,:));
    ncs_lin = all_lin{LOOSub}(channel_selections); % Apply channel selection
    ncs_snr = all_snr{LOOSub}(channel_selections); % Apply channel selection
    ncs_proc = all_ncs_proc{LOOSub}(:, channel_selections); % Apply channel selection
    LVCath = all_LVCath{LOOSub};
    ncs_lin_rescale = rescale(ncs_lin, optimal_parameter(1, i), optimal_parameter(2, i));
    ncs_snr_rescale = rescale(ncs_snr, optimal_parameter(3, i), optimal_parameter(4, i));
    ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;
    ncs_proc_weight = zeros(length(ncs_proc(:,1)), length(ncs_snr));
    for ch = 1:length(ncs_snr)
        ncs_proc_weight(:, ch) = ncs_proc(:,ch) .* ncs_linsnr_weight(ch);
    end
    ncs_proc_weight_mean = sum(ncs_proc_weight, 2);
    ncs_proc_weight_mean = zscore(ncs_proc_weight_mean / sum(ncs_linsnr_weight));
    all_weightedCorrs(i) = corr(ncs_proc_weight_mean, LVCath);
end
sprintf('%.2f +- %.2f', mean(all_weightedCorrs), std(all_weightedCorrs))


% Correlation histograms
be = -1:0.1:1;
mean_corrs = mean(all_corrs, 2);
max_corrs = max(all_corrs);

% Fit Gaussian distribution to the data
g_all = fitdist(all_corrs(:), 'Normal');
g_mean = fitdist(mean_corrs(:), 'Normal');
g_max = fitdist(max_corrs(:), 'Normal');
g_snr = fitdist(all_weightedCorrs(:), 'Normal');

% Define the range for the x-axis
x = linspace(min(be), max(be), 100);

% Calculate the Gaussian distribution for each data set
yall = pdf(g_all, x);
ymean = pdf(g_mean, x);
ymax = pdf(g_max, x);
ysnr = pdf(g_snr, x);

% Calculate medians for each data set
mdall = median(all_corrs(:));
mdmax = median(max_corrs);
mdmean = median(mean_corrs);
mdsnr = median(all_weightedCorrs);

figure('Position', [100 0 600 800])

% Plot histogram of PCA component correlations
subplot(2, 1, 1);
histogram(all_weightedCorrs, 'Normalization', 'pdf', 'FaceColor', 'r', 'BinEdges', be); hold on;
plot(x, ysnr, 'r-', 'LineWidth', 2);
xline(mdsnr, 'r--', 'LineWidth', 2); hold off;
title('Histogram of Weighted Average Correlations');
xlabel('Correlation'); xlim([be(1) be(end)]);
ylabel('Probability');
legend('Weighted Corr', sprintf('\\mu=%.2f, \\sigma=%.2f', g_snr.mu, g_snr.sigma), sprintf('Med.=%.2f', mdsnr),...
    'Location', 'northwest', 'AutoUpdate', 'off');

% Plot histogram of channel-wise correlations
subplot(2, 1, 2);
histogram(mean_corrs, 'Normalization', 'pdf', 'FaceColor', 'g', 'BinEdges', be); hold on;
plot(x, ymean, 'g-', 'LineWidth', 2);
xline(mdmean, 'g--', 'LineWidth', 2);
histogram(max(all_corrs), 'Normalization', 'pdf', 'FaceColor', 'm', 'BinEdges', be);
plot(x, ymax, 'm-', 'LineWidth', 2);
xline(mdmax, 'm--', 'LineWidth', 2);
histogram(all_corrs(:), 'Normalization', 'pdf', 'FaceColor', 'k', 'BinEdges', be);
plot(x, yall, 'k-', 'LineWidth', 2);
xline(mdall, 'k--', 'LineWidth', 2);

hold off;
title('Histogram of Channel-Wise Correlations');
xlabel('Mean Correlation'); xlim([be(1) be(end)]);
ylabel('Probability');
legend('Mean Ch Corr', sprintf('\\mu=%.2f, \\sigma=%.2f', g_mean.mu, g_mean.sigma), sprintf('Med.=%.2f', mdmean),...
    'Max Ch Corr',  sprintf('\\mu=%.2f, \\sigma=%.2f', g_max.mu, g_max.sigma), sprintf('Med.=%.2f', mdmax),...
    'All Ch Corr',  sprintf('\\mu=%.2f, \\sigma=%.2f', g_all.mu, g_all.sigma), sprintf('Med.=%.2f', mdall),...
    'Location', 'northwest', 'AutoUpdate', 'off');


