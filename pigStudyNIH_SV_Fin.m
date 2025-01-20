%% pigStudyNIH_SV_Fin.m

filePath = ['']; % insert file path here
inc = 1;
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

fileDate = '02-08-24';  routines = {'dexmed'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 4;
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

fileDate = '03-07-24';  routines = {'dexmed', 'phenylephrine2'};
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
        case '02-05-24'
            switch interventionType
                case 'dexmed'
                    all_bnds{i} = {[51161 51179], [51226 51243]};
                case 'dobutamine'
                    all_bnds{i} = {[50338 50390]};
                case 'phenylephrine-2'
                    all_bnds{i} = {[49812 49925]};
            end
        case '02-06-24'
            switch interventionType
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
                    all_bnds{i} = {[57035 57048], [57060 57071], [57079 57086]};
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
                case 'venacavaocclusion2'
                    all_bnds{i} = {[49398 49436]};
            end
        case '03-07-24'
            switch interventionType
                case 'dexmed'
                    all_bnds{i} = {[44701 44773]};
                case 'dobutamine'
                    all_bnds{i} = {[44348 44418]};
                case 'phenylephrine2'
                    all_bnds{i} = {[43790 43860]};
                case 'cavalocclusion'
                    all_bnds{i} = {[46318 46363]};
            end
    end
end
clearvars -except all*

% Analyze P-P per-routine - extract all P-Ps
all_data = {}; inc = 1; 
for i = 1:length(all_ncs)
    win = 21;
    ncs = all_ncs{i}; bio = all_bio{i}; rt = all_ints{i}; dt = all_dates{i};
    tBnds = all_bnds{i}; power_sync = all_power{i};
    if contains(rt, 'dex') || contains(rt, 'phenyl')
        tPP = []; tNcsSig = []; tBioSig = [];
        bioS1_PP = []; bioLV_PP = []; bioFem_PP = []; 
        bioS1_sig = []; bioLV_sig = []; bioFem_sig = [];
        ncs_PP = cell(1, length(ncs.chSel));
        ncs_sig = cell(1, length(ncs.chSel));
        PPvar = ncs.PP_BP; Fem_SV = bio.SV_Fem(1:end-1);
        % apply arrhythmia mask
        idxArr = ECGArrhythmiaCheck(bio);
        for j = 1:length(ncs.chSel)
            chsel = ncs.chSel(j);
            PPvar{chsel} = PPvar{chsel}(idxArr);
        end
        S1_PP_arr_rem = bio.S1LV_SV(idxArr);
        LV_PP_arr_rem = bio.LV_SV(idxArr);
        Fem_SV_arr_rem = Fem_SV(idxArr);
        tPP_arrRem = ncs.t_allWaveforms(idxArr);
        % process within each bound
        filtType = 'mean'; win_filt = [];
        for t = 1:length(tBnds)
            % cut signal in the selected bounds
            mn = tBnds{t}(1); mx = tBnds{t}(2);
            idx = tPP_arrRem > mn & tPP_arrRem < mx;
            idxNcsSig = ncs.tNcsDS > mn & ncs.tNcsDS < mx;
            idxBioSig = bio.tBio > mn & bio.tBio < mx;
            tPPinterp = mn:0.5:mx;
            % for each channel, cut, remove outliers, interpolate, and filter
            for j = 1:length(ncs.chSel)
                chsel = ncs.chSel(j);
                ncs_temp = PPvar{chsel}(idx);
                tPP_temp = tPP_arrRem(idx);

                [idx_outlier, idxFlag(j)] = PP_removeOutliers_V2(ncs_temp, tPP_temp, ncs.ncsDS(:, chsel), ncs.tNcsDS);

                win_filt(t) = min([sum(~idx_outlier), win]);
                ncs_temp = interp1_customextrap(tPP_temp(~idx_outlier), ncs_temp(~idx_outlier), tPPinterp, 'linear', [ncs_temp(1) ncs_temp(end)]);
                if strcmp(filtType, 'mean')
                    b = (1/win_filt(t))*ones(1,win_filt(t)); a = 1;
                    foo = [ones(1, win_filt(t))*mean(ncs_temp(1:win_filt(t))) ncs_temp];
                    foo = filter(b, a, foo);
                    ncs_PP{j} = [ncs_PP{j} foo(win_filt(t)+1:end)];
                elseif strcmp(filtType, 'med')
                    ncs_PP{j} = [ncs_PP{j} medfilt1(ncs_temp, win_filt(t), 'truncate')];
                end
            end
            tPP = [tPP tPPinterp];

            % for conductance cath PPs, remove outliers and interpoalte and filter
            S1_temp = S1_PP_arr_rem(idx);
            idx_outlier = isoutlier(S1_temp, 'movmedian', 10);
            S1_temp = interp1_customextrap(tPP_temp(~idx_outlier), S1_temp(~idx_outlier), tPPinterp, 'linear', [S1_temp(1) S1_temp(end)]);
            LVV_temp = LV_PP_arr_rem(idx);
            idx_outlier = isoutlier(LVV_temp, 'movmedian', 10);
            LVV_temp = interp1_customextrap(tPP_temp(~idx_outlier), LVV_temp(~idx_outlier), tPPinterp, 'linear', [LVV_temp(1) LVV_temp(end)]);
            FemSV_temp = Fem_SV_arr_rem(idx);
            idx_outlier = isoutlier(FemSV_temp, 'movmedian', 10);
            FemSV_temp = interp1_customextrap(tPP_temp(~idx_outlier), FemSV_temp(~idx_outlier), tPPinterp, 'linear', [FemSV_temp(1) FemSV_temp(end)]);
            if strcmp(filtType, 'mean')
                foo = [ones(1, win_filt(t))*mean(S1_temp(1:win_filt(t))) S1_temp]; foo = filter(b, a, foo);
                bioS1_PP = [bioS1_PP foo(win_filt(t)+1:end)];
                foo = [ones(1, win_filt(t))*mean(LVV_temp(1:win_filt(t))) LVV_temp]; foo = filter(b, a, foo);
                bioLV_PP = [bioLV_PP foo(win_filt(t)+1:end)];
                foo = [ones(1, win_filt(t))*mean(FemSV_temp(1:win_filt(t))) FemSV_temp]; foo = filter(b, a, foo);
                bioFem_PP = [bioFem_PP foo(win_filt(t)+1:end)];
            elseif strcmp(filtType, 'med')
                bioS1_PP = [bioS1_PP medfilt1(S1_temp, win_filt(t), 'truncate')];
                bioLV_PP = [bioLV_PP medfilt1(LVV_temp, win_filt(t), 'truncate')];
                bioFem_PP = [bioFem_PP medfilt1(FemSV_temp, win_filt(t), 'truncate')];
            end

            % store synced ncs and gt signals
            for j = 1:length(ncs.chSel)
                chsel = ncs.chSel(j);
                ncs_sig{j} = [ncs_sig{j}; ncs.ncsBPFilt(idxNcsSig, chsel)];
            end
            bioS1_sig = [bioS1_sig power_sync.data{4}(idxBioSig)];
            bioLV_sig = [bioLV_sig bio.LV_vol(idxBioSig)];
            bioFem_sig = [bioFem_sig power_sync.data{9}(idxBioSig)];
            tNcsSig = [tNcsSig ncs.tNcsDS(idxNcsSig)];
            tBioSig = [tBioSig bio.tBio(idxBioSig)];
        end

        % store
        all_data{inc}.ncs_PP = ncs_PP;
        all_data{inc}.tPP = tPP;
        all_data{inc}.tPP_arrRem = tPP_arrRem;
        all_data{inc}.bioLV_PP = bioLV_PP;
        all_data{inc}.bioS1_PP = bioS1_PP;
        all_data{inc}.bioFem_PP = bioFem_PP;
        all_data{inc}.rt = rt;
        all_data{inc}.dt = dt;
        all_data{inc}.ncs_sig = ncs_sig;
        all_data{inc}.bioS1_sig = bioS1_sig;
        all_data{inc}.bioLV_sig = bioLV_sig;
        all_data{inc}.bioFem_sig = bioFem_sig;
        all_data{inc}.tNcsSig = tNcsSig;
        all_data{inc}.tBioSig = tBioSig;

        all_data{inc}.idxChannelFlag = idxFlag;

        all_data{inc}.win_filt = win_filt;
        all_data{inc}.tBnds = tBnds;
        all_data{inc}.t_allECGWaveforms = bio.t_allECGWaveforms;
        all_data{inc}.fNcsDS = ncs.fNcsDS;
        inc = inc+1;
    end
    clearvars -except all* inc
end
clearvars -except all*

% Calibrate all peak-to-peaks and calculate SNR/Linearity
clearvars -except all*
inc = 1;
for i = 1:length(all_data)
    % calibrate all PPs and normalize signals: calibrated to the mean of the first window
    win = all_data{i}.win_filt(1);
    LV_PP_calib_level = mean(all_data{i}.bioLV_PP(1:win)); 
    bioLV_PP_calib = all_data{i}.bioLV_PP./LV_PP_calib_level;
    bioLV_sig_norm = zscore(all_data{i}.bioLV_sig); 

    Fem_PP_calib_level = mean(all_data{i}.bioFem_PP(1:win)); 
    bioFem_PP_calib = all_data{i}.bioFem_PP./Fem_PP_calib_level;
    bioFem_sig_norm = zscore(all_data{i}.bioFem_sig); 

    S1_PP_calib_level = mean(all_data{i}.bioS1_PP(1:win)); 
    bioS1_PP_calib = all_data{i}.bioS1_PP./S1_PP_calib_level;
    bioS1_sig_norm = zscore(all_data{i}.bioS1_sig);
    for j = 1:length(all_data{i}.ncs_PP)
        ncs_PP_calib_level = mean(all_data{i}.ncs_PP{j}(1:win)); 
        ncs_PP_calib{j} = all_data{i}.ncs_PP{j}./ncs_PP_calib_level; 
        ncs_sig_norm{j} = zscore(all_data{i}.ncs_sig{j});
    end

    % calculate snr/linearity of NCS signals
    mn = all_data{i}.tBnds{1}(1); mx = all_data{i}.tBnds{end}(2); 
    idx = all_data{i}.tPP_arrRem > mn & all_data{i}.tPP_arrRem < mx;
    RR = diff(all_data{i}.t_allECGWaveforms); RR = [RR inf]; RR = RR(idx);
    idx_outlier = isoutlier(RR); RR = RR(~idx_outlier);
    RRbnds = [0.95*min(RR) 1.05*max(RR)];
    ncs_snr = []; ncs_lin = [];
    for j = 1:length(all_data{i}.ncs_PP)
        try
            if all_data{inc}.idxChannelFlag(j)
                error('channel flag');
            end
            ncs_snr(j) = signalSNR(ncs_sig_norm{j}, sort(1./RRbnds), all_data{i}.fNcsDS, [3 10]);
            ncs_lin(j) = signalLinearity(ncs_sig_norm{j}, sort(1./RRbnds), all_data{i}.fNcsDS);
        catch
            ncs_snr(j) = 0; ncs_lin(j) = 0;
        end
    end

    % store relevant data
    all_data{inc}.ncs_PP_calib = ncs_PP_calib;
    all_data{inc}.bioLV_PP_calib = bioLV_PP_calib;
    all_data{inc}.bioS1_PP_calib = bioS1_PP_calib;
    all_data{inc}.bioFem_PP_calib = bioFem_PP_calib;

    all_data{inc}.LV_PP_calib_level = LV_PP_calib_level;
    all_data{inc}.Fem_PP_calib_level = Fem_PP_calib_level;
    all_data{inc}.S1_PP_calib_level = S1_PP_calib_level;

    all_data{inc}.ncs_snr = ncs_snr;
    all_data{inc}.ncs_lin = ncs_lin;

    all_data{inc}.bioLV_sig_norm = bioLV_sig_norm;
    all_data{inc}.bioS1_sig_norm = bioS1_sig_norm;
    all_data{inc}.bioFem_sig_norm = bioFem_sig_norm;
    all_data{inc}.ncs_sig_norm = ncs_sig_norm;

    inc = inc+1;
    clearvars -except all* inc
end
clearvars -except all*

% Gather data for plotting and selecting appropriate S1 volume segments
inc = 1;
for i = 1:length(all_data)
    % replace routines where we use S1 instead of composite volume here
    data_ncs_PP{inc} = all_data{i}.ncs_PP_calib;
    data_ncs_snr{inc} = all_data{i}.ncs_snr;
    data_ncs_lin{inc} = all_data{i}.ncs_lin;
    data_win{inc} = all_data{i}.win_filt;
    data_tPP{inc} = all_data{i}.tPP;
    data_rt{inc} = [all_data{i}.dt all_data{i}.rt];
    if (contains(all_data{i}.dt, '03-05') && contains(all_data{i}.rt, 'ph')) ||...
        (contains(all_data{i}.dt, '03-07') && contains(all_data{i}.rt, 'ph')) ||...
        (contains(all_data{i}.dt, '03-05') && contains(all_data{i}.rt, 'dex')) ||...
        (contains(all_data{i}.dt, '03-07') && contains(all_data{i}.rt, 'dex'))  
        data_bio_PP_calib{inc} = all_data{i}.bioS1_PP_calib;
        data_bio_PP_calib_level{inc} = all_data{i}.S1_PP_calib_level;
    else
        data_bio_PP_calib{inc} = all_data{i}.bioLV_PP_calib;
        data_bio_PP_calib_level{inc} = all_data{i}.LV_PP_calib_level;
    end
    inc = inc+1;
end
clearvars -except all* data*

% Use CV to find optimal parameters for scaling P-P: selected channels only
interventions = 1:11; n = 10;
CV_combinations = nchoosek(interventions, n);
chs = 1:8;

% define parameter search
parameter_range = 0:0.2:2;
[x1, x2, x3, x4] = ndgrid(parameter_range, parameter_range, parameter_range, parameter_range);
grid_search = [x1(:), x2(:), x3(:), x4(:)]';
idx_grid = logical(1:size(grid_search, 2));
for j = 1:size(grid_search, 2)
    if grid_search(1, j) > grid_search(2, j) || grid_search(3, j) > grid_search(4, j)
        idx_grid(j) = false;
    end
end
grid_search = grid_search(:, idx_grid);

optimal_parameter = []; train_err = []; test_err = [];
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    mean_mae_end_CV = [];
    for j = 1:size(grid_search, 2) % search over possible hyperparameters
        mae_end_CV = [];
        for k = 1:size(CV_combinations, 2) % loop over each index in CV combinations with a set of hyperparameters 
            ncs_lin = data_ncs_lin{CV_combinations(i, k)}(chs); 
            ncs_snr = data_ncs_snr{CV_combinations(i, k)}(chs); 
            ncs_PP_calib = data_ncs_PP{CV_combinations(i, k)}(chs); 
            bio_PP_calib = data_bio_PP_calib{CV_combinations(i, k)};
            win = data_win{CV_combinations(i, k)};
            
            ncs_lin_rescale = rescale(ncs_lin, grid_search(1, j), grid_search(2, j));
            ncs_snr_rescale = rescale(ncs_snr, grid_search(3, j), grid_search(4, j));
            ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;

            ncs_PP_weight = zeros(length(ncs_snr), length(ncs_PP_calib{1}));
            for ch = 1:length(ncs_PP_calib)
                ncs_PP_weight(ch, :) = ncs_PP_calib{ch} .* ncs_linsnr_weight(ch);
            end
            ncs_PP_weight_mean = sum(ncs_PP_weight, 1);
            ncs_PP_weight_mean = ncs_PP_weight_mean / sum(ncs_linsnr_weight);
            ncs_PP_weight_mean = ncs_PP_weight_mean/mean(ncs_PP_weight_mean(1:win(1))); % calibrate
            % optimizing for mean absolute error
            mae_end_CV(k) = similarity(ncs_PP_weight_mean(round(0.9*length(ncs_PP_weight_mean)):end), ...
                bio_PP_calib(round(0.9*length(bio_PP_calib)):end), 'mae');
        end
        mean_mae_end_CV(j) = mean(mae_end_CV); 
    end
    [m, l] = min(mean_mae_end_CV);
    optimal_parameter(:, i) = grid_search(:, l);
end

clearvars -except all* data* optimal_parameter CV_combinations interventions chs

% Evaluate
all_SV_error = [];
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    % Test LOO mean error with optimal parameters
    LOOSub = setdiff(interventions, CV_combinations(i,:));
    ncs_lin = data_ncs_lin{LOOSub};
    ncs_snr = data_ncs_snr{LOOSub}; 
    calib_level = data_bio_PP_calib_level{LOOSub};
    ncs_PP_calib = data_ncs_PP{LOOSub};
    bio_PP_calib = data_bio_PP_calib{LOOSub}*calib_level; 
    win = data_win{LOOSub};
    ncs_lin_rescale = rescale(ncs_lin, optimal_parameter(1, i), optimal_parameter(2, i));
    ncs_snr_rescale = rescale(ncs_snr, optimal_parameter(3, i), optimal_parameter(4, i));
    ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;
    ncs_PP_weight = zeros(length(ncs_snr), length(ncs_PP_calib{1}));
    for ch = 1:length(ncs_PP_calib)
        ncs_PP_weight(ch, :) = ncs_PP_calib{ch} .* ncs_linsnr_weight(ch);
    end
    ncs_PP_weight_mean = sum(ncs_PP_weight, 1);
    ncs_PP_weight_mean = ncs_PP_weight_mean / sum(ncs_linsnr_weight);
    ncs_PP_weight_mean = ncs_PP_weight_mean/mean(ncs_PP_weight_mean(1:win(1)));
    ncs_PP_weight_mean = ncs_PP_weight_mean*calib_level;
    all_SV_error(i) = similarity(ncs_PP_weight_mean(round(0.9*length(ncs_PP_weight_mean)):end), ...
        bio_PP_calib(round(0.9*length(bio_PP_calib)):end), 'merr');
end

sprintf('%.2f +- %.2f', mean(all_SV_error), std(all_SV_error))
[~, q] = iqr(all_SV_error); 
sprintf('Median %.2f (IQR %.2f - %.2f)', median(all_SV_error), q(1), q(2))