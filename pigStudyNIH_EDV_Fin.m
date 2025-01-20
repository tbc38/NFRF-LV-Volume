%% pigStudyNIH_EDV_Fin.m

filePath = ['']; % insert file path here
inc = 1;

fileDate = '01-26-24'; routines = {'venacavaocclusion2'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 0; 
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1; 
    clearvars tBnds
end

fileDate = '02-05-24'; routines = {'venacavaocclusion'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 1; 
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1; 
    clearvars tBnds
end

fileDate = '02-06-24'; routines = {'venacavaocclusion'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 2;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '02-07-24'; routines = {'venacavaocclusion'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 3;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '02-08-24';  routines = {'venacavaocclusion'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 4;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end
 
fileDate = '03-04-24'; routines = {'venacavaocclusion4'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 5;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '03-05-24'; routines = {'venacavaocclusion2'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 6;
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '03-06-24'; routines = {'venacavaocclusion2'};
for rt = routines
    load([filePath fileDate ' ' rt{1} version '.mat'], 'ncs', 'bio', 'power_sync');
    all_ncs{inc} = ncs; all_bio{inc} = bio; all_subs{inc} = 7; 
    all_dates{inc} = fileDate; all_ints{inc} = rt{1};
    all_power{inc} = power_sync;
    inc = inc+1;
    clearvars tBnds
end

fileDate = '03-07-24';  routines = {'cavalocclusion'};
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
                case 'venacavaocclusion2'
                    all_bnds{i} = {[41771 41815]};
                case 'venacavaocclusion3'
                    all_bnds{i} = {[46640 46741]};
            end
        case '02-05-24'
            switch interventionType
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

% Analyze EDV & polarity per-routine - extract all 
all_data = {}; inc = 1; 
for i = 1:length(all_ncs)
    win = 21;
    ncs = all_ncs{i}; bio = all_bio{i}; rt = all_ints{i}; dt = all_dates{i};
    tBnds = all_bnds{i}; power_sync = all_power{i};
    if contains(rt, 'occ')
        t_EDV = []; tNcsSig = []; tBioSig = [];
        bioS1_PP = []; bioLV_PP = [];
        bioS1_EDV = []; bioLV_EDV = [];
        bioS1_sig = []; bioLV_sig = [];
        
        ncs_waveforms = cell(1, length(ncs.chSel));
        T_Wave_waveformidx = cell(1, length(ncs.chSel));
        t_waveforms = cell(1, length(ncs.chSel));
        ncs_EDV = cell(1, length(ncs.chSel));
        ncs_PP = cell(1, length(ncs.chSel));
        ncs_sig = cell(1, length(ncs.chSel));
        ncs_BP_sig = cell(1, length(ncs.chSel));
        EDVvar = ncs.EDV_DS; PPvar = ncs.PP_BP;
        ncsWaveformVar = ncs.allBPWaveforms;
        % apply arrhythmia mask
        idxArr = ECGArrhythmiaCheck(bio, 0.25);
        for j = 1:length(ncs.chSel)
            chsel = ncs.chSel(j);
            EDVvar{chsel} = EDVvar{chsel}(idxArr);
            PPvar{chsel} = PPvar{chsel}(idxArr);
            ncsWaveformVar{chsel} = ncsWaveformVar{chsel}(:, idxArr);
        end
        S1_PP_arr_rem = bio.S1LV_SV(idxArr);
        LV_PP_arr_rem = bio.LV_SV(idxArr);
        S1_EDV_arr_rem = bio.S1LV_EDV(idxArr);
        LV_EDV_arr_rem = bio.LV_EDV(idxArr);
        twave_idx = bio.allidxTWaves(idxArr);
        t_arrRem = ncs.t_allWaveforms(idxArr);
        filtType = 'mean'; 
        % cut signal in the selected bounds
        mn = tBnds{1}(1); mx = tBnds{1}(2);
        idx = t_arrRem > mn & t_arrRem < mx;
        idxNcsSig = ncs.tNcsDS > mn & ncs.tNcsDS < mx;
        idxBioSig = bio.tBio > mn & bio.tBio < mx;
        t_interp = mn:0.5:mx;
        % for each channel, cut, remove outliers, interpolate, and filter
        for j = 1:length(ncs.chSel)
            chsel = ncs.chSel(j);
            ncsEDV_temp = EDVvar{chsel}(idx);
            ncsPP_temp = PPvar{chsel}(idx);
            ncsWaveform_temp = ncsWaveformVar{chsel}(:, idx);
            tWave_temp = twave_idx(idx);
            t_temp = t_arrRem(idx);

            [idx_outlier, idxFlag(j)] = PP_removeOutliers_V2(ncsPP_temp, t_temp, ncs.ncsDS(:, chsel), ncs.tNcsDS);

            win_filt = min([sum(~idx_outlier), win]);
            ncsPP_temp = interp1_customextrap(t_temp(~idx_outlier), ncsPP_temp(~idx_outlier), t_interp, 'linear', [ncsPP_temp(1) ncsPP_temp(end)]);
            ncsEDV_temp = interp1_customextrap(t_temp(~idx_outlier), ncsEDV_temp(~idx_outlier), t_interp, 'linear', [ncsEDV_temp(1) ncsEDV_temp(end)]);
            ncs_waveforms{j} = ncsWaveform_temp(:, ~idx_outlier);
            t_waveforms{j} = t_temp(~idx_outlier);
            T_Wave_waveformidx{j} = tWave_temp(~idx_outlier);
            if strcmp(filtType, 'mean')
                b = (1/win_filt)*ones(1,win_filt); a = 1;
                foo = [ones(1, win_filt)*mean(ncsPP_temp(1:win_filt)) ncsPP_temp];
                foo = filter(b, a, foo);
                ncs_PP{j} = foo(win_filt+1:end);
                foo = [ones(1, win_filt)*mean(ncsEDV_temp(1:win_filt)) ncsEDV_temp];
                foo = filter(b, a, foo);
                ncs_EDV{j} = foo(win_filt+1:end);
            elseif strcmp(filtType, 'med')
                ncs_EDV{j} = medfilt1(ncsEDV_temp, win_filt, 'truncate');
                ncs_PP{j} = medfilt1(ncsPP_temp, win_filt, 'truncate');
            end
        end
        % store interpolated time axis
        t_EDV = t_interp;

        % for conductance cath waveforms, remove outliers and
        % interpoalte and filter PPs
        S1_temp = S1_PP_arr_rem(idx);
        idx_outlier = isoutlier(S1_temp, 'movmedian', 10);
        S1_temp = interp1_customextrap(t_temp(~idx_outlier), S1_temp(~idx_outlier), t_interp, 'linear', [S1_temp(1) S1_temp(end)]);
        LVV_temp = LV_PP_arr_rem(idx);
        idx_outlier = isoutlier(LVV_temp, 'movmedian', 10);
        LVV_temp = interp1_customextrap(t_temp(~idx_outlier), LVV_temp(~idx_outlier), t_interp, 'linear', [LVV_temp(1) LVV_temp(end)]);
        if strcmp(filtType, 'mean')
            foo = [ones(1, win_filt)*mean(S1_temp(1:win_filt)) S1_temp]; foo = filter(b, a, foo);
            bioS1_PP = foo(win_filt+1:end);
            foo = [ones(1, win_filt)*mean(LVV_temp(1:win_filt)) LVV_temp]; foo = filter(b, a, foo);
            bioLV_PP = foo(win_filt+1:end);
        elseif strcmp(filtType, 'med')
            bioS1_PP = medfilt1(S1_temp, win_filt, 'truncate');
            bioLV_PP = medfilt1(LVV_temp, win_filt, 'truncate');
        end

        % for conductance cath waveforms, remove outliers and
        % interpoalte and filter EDVs
        S1_temp = S1_EDV_arr_rem(idx);
        idx_outlier = isoutlier(S1_temp, 'movmedian', 10);
        S1_temp = interp1_customextrap(t_temp(~idx_outlier), S1_temp(~idx_outlier), t_interp, 'linear', [S1_temp(1) S1_temp(end)]);
        LVV_temp = LV_EDV_arr_rem(idx);
        idx_outlier = isoutlier(LVV_temp, 'movmedian', 10);
        LVV_temp = interp1_customextrap(t_temp(~idx_outlier), LVV_temp(~idx_outlier), t_interp, 'linear', [LVV_temp(1) LVV_temp(end)]);
        if strcmp(filtType, 'mean')
            foo = [ones(1, win_filt)*mean(S1_temp(1:win_filt)) S1_temp]; foo = filter(b, a, foo);
            bioS1_EDV = foo(win_filt+1:end);
            foo = [ones(1, win_filt)*mean(LVV_temp(1:win_filt)) LVV_temp]; foo = filter(b, a, foo);
            bioLV_EDV = foo(win_filt+1:end);
        elseif strcmp(filtType, 'med')
            bioS1_EDV = medfilt1(S1_temp, win_filt, 'truncate');
            bioLV_EDV = medfilt1(LVV_temp, win_filt, 'truncate');
        end
        
        % store synced ncs and gt signals
        for j = 1:length(ncs.chSel)
            chsel = ncs.chSel(j);
            ncs_sig{j} = ncs.ncsDS(idxNcsSig, chsel);
            ncs_BP_sig{j} = ncs.ncsBPFilt(idxNcsSig, chsel);
        end
        bioS1_sig = power_sync.data{4}(idxBioSig);
        bioLV_sig = bio.LV_vol(idxBioSig);
        tNcsSig = ncs.tNcsDS(idxNcsSig);
        tBioSig = bio.tBio(idxBioSig);
        
        % store
        all_data{inc}.ncs_EDV = ncs_EDV;
        all_data{inc}.t_EDV = t_EDV;
        all_data{inc}.t_waveforms = t_waveforms;
        all_data{inc}.t_arrRem = t_arrRem;
        all_data{inc}.idxArr = idxArr;
        all_data{inc}.bioLV_EDV = bioLV_EDV;
        all_data{inc}.bioS1_EDV = bioS1_EDV;
        all_data{inc}.T_Wave_Idxs = T_Wave_waveformidx;

        all_data{inc}.idxChannelFlag = idxFlag;
        
        all_data{inc}.ncs_PP = ncs_PP;
        all_data{inc}.bioLV_PP = bioLV_PP;
        all_data{inc}.bioS1_PP = bioS1_PP;

        all_data{inc}.rt = rt;
        all_data{inc}.dt = dt;
        all_data{inc}.ncs_sig = ncs_sig;
        all_data{inc}.ncs_BP_sig = ncs_BP_sig;
        all_data{inc}.bioS1_sig = bioS1_sig;
        all_data{inc}.bioLV_sig = bioLV_sig;
        all_data{inc}.tNcsSig = tNcsSig;
        all_data{inc}.tBioSig = tBioSig;
        all_data{inc}.win_filt = win_filt;
        all_data{inc}.tBnds = tBnds;
        all_data{inc}.t_allECGWaveforms = bio.t_allECGWaveforms;
        all_data{inc}.fNcsDS = ncs.fNcsDS;
        all_data{inc}.ncs_waveforms = ncs_waveforms;
        inc = inc+1;
    end
    clearvars -except all* inc
end
clearvars -except all*

% Calibrate signal per-routine - extract all 
clearvars -except all*
inc = 1;
for i = 1:length(all_data)
    win = all_data{i}.win_filt(1);
    % LVCath
    LV_sig_calib_level = mean(all_data{i}.bioLV_sig(1:win)); 
    bioLV_sig_calib = all_data{i}.bioLV_sig./LV_sig_calib_level;
    bioLV_sig_norm = zscore(all_data{i}.bioLV_sig); 

    LV_PP_calib_level = mean(all_data{i}.bioLV_PP(1:win)); 
    bioLV_PP_calib = all_data{i}.bioLV_PP./LV_PP_calib_level;
    bioLV_PP_norm = zscore(all_data{i}.bioLV_PP);

    bioLV_EDV_PP_calib = ...
        (all_data{i}.bioLV_EDV - mean(all_data{i}.bioLV_EDV(1:win)))./LV_PP_calib_level;

    % S1cath
    S1_sig_calib_level = mean(all_data{i}.bioS1_sig(1:win)); 
    bioS1_sig_calib = all_data{i}.bioS1_sig./S1_sig_calib_level;
    bioS1_sig_norm = zscore(all_data{i}.bioS1_sig); 

    S1_PP_calib_level = mean(all_data{i}.bioS1_PP(1:win)); 
    bioS1_PP_calib = all_data{i}.bioS1_PP./S1_PP_calib_level;
    bioS1_PP_norm = zscore(all_data{i}.bioS1_PP);
    
    bioS1_EDV_PP_calib = ...
        (all_data{i}.bioS1_EDV - mean(all_data{i}.bioS1_EDV(1:win)))./S1_PP_calib_level;

    %RF
    for j = 1:length(all_data{i}.ncs_sig)
        ncs_sig_calib_level{j} = mean(all_data{i}.ncs_sig{j}(1:win)); 
        ncs_sig_calib{j} = all_data{i}.ncs_sig{j}./ncs_sig_calib_level{j}; 
        ncs_sig_norm{j} = zscore(all_data{i}.ncs_sig{j});
        ncs_BP_sig_norm{j} = zscore(all_data{i}.ncs_BP_sig{j});

        ncs_PP_calib_level{j} = mean(all_data{i}.ncs_PP{j}(1:win)); 
        ncs_PP_calib{j} = all_data{i}.ncs_PP{j}./ncs_PP_calib_level{j}; 
        ncs_PP_norm{j} = zscore(all_data{i}.ncs_PP{j});

        ncs_EDV_PP_calib{j} = ...
            (all_data{i}.ncs_EDV{j} - mean(all_data{i}.ncs_EDV{j}(1:win)))./ncs_PP_calib_level{j};
    end

    % store relevant data
    all_data{inc}.ncs_sig_calib = ncs_sig_calib;
    all_data{inc}.ncs_BP_sig_norm = ncs_BP_sig_norm;
    all_data{inc}.bioLV_sig_calib = bioLV_sig_calib;
    all_data{inc}.bioS1_sig_calib = bioS1_sig_calib;
    all_data{inc}.ncs_sig_norm = ncs_sig_norm;
    all_data{inc}.bioLV_sig_norm = bioLV_sig_norm;
    all_data{inc}.bioS1_sig_norm = bioS1_sig_norm;
    all_data{inc}.LV_sig_calib_level = LV_sig_calib_level;
    all_data{inc}.S1_sig_calib_level = S1_sig_calib_level;
    all_data{inc}.ncs_sig_calib_level = ncs_sig_calib_level;
    
    all_data{inc}.ncs_EDV_PP_calib = ncs_EDV_PP_calib;
    all_data{inc}.bioLV_EDV_PP_calib = bioLV_EDV_PP_calib;
    all_data{inc}.bioS1_EDV_PP_calib = bioS1_EDV_PP_calib;

    all_data{inc}.ncs_PP_calib = ncs_PP_calib;
    all_data{inc}.bioLV_PP_calib = bioLV_PP_calib;
    all_data{inc}.bioS1_PP_calib = bioS1_PP_calib;
    all_data{inc}.ncs_PP_norm = ncs_PP_norm;
    all_data{inc}.bioLV_PP_norm = bioLV_PP_norm;
    all_data{inc}.bioS1_PP_norm = bioS1_PP_norm;
    all_data{inc}.LV_PP_calib_level = LV_PP_calib_level;
    all_data{inc}.S1_PP_calib_level = S1_PP_calib_level;
    all_data{inc}.ncs_PP_calib_level = ncs_PP_calib_level;

    inc = inc+1;
    clearvars -except all* inc
end
clearvars -except all*

% Analyze calibrated signal per-routine - extract all 
clearvars -except all*
inc = 1;
for i = 1:length(all_data)

    t = all_data{i}.t_allECGWaveforms(all_data{i}.idxArr);
    % calculate snr/linearity of NCS signals
    mn = all_data{i}.tBnds{1}(1); mx = all_data{i}.tBnds{end}(2); 
    idx = t > mn & t < mx;
    RR = diff(t); RR = [RR inf]; RR = RR(idx);
    idx_outlier = isoutlier(RR); RR = RR(~idx_outlier);
    RRbnds = [0.9*min(RR) 1.05*max(RR)];
    ncs_snr = []; ncs_lin = [];
    for j = 1:length(all_data{i}.ncs_EDV)
        try
            if all_data{inc}.idxChannelFlag(j)
                error('channel flag');
            end
            ncs_snr(j) = signalSNR(all_data{inc}.ncs_sig_norm{j}, sort(1./RRbnds), all_data{i}.fNcsDS, [3 10]);
            ncs_lin(j) = signalLinearity(all_data{inc}.ncs_sig_norm{j}, sort(1./RRbnds), all_data{i}.fNcsDS);
        catch
            ncs_snr(j) = 0; ncs_lin(j) = 0;
        end
    end

    % mL change
    bioLV_EDV_mldiff = range(all_data{inc}.bioLV_EDV);
    bioS1_EDV_mldiff = range(all_data{inc}.bioS1_EDV);

    % change denominated in number of P-Ps
    bioLV_EDV_PPdiff = range(all_data{inc}.bioLV_EDV_PP_calib);
    bioS1_EDV_PPdiff = range(all_data{inc}.bioS1_EDV_PP_calib);
    ncs_EDV_PPdiff = [];
    for j = 1:length(all_data{i}.ncs_EDV)
        ncs_EDV_PPdiff(j) = range(all_data{i}.ncs_EDV_PP_calib{j});
    end

    % estimate polarity of EDV decrease
    ncs_EDV_polarity = [];
    [~, l] = min(all_data{i}.bioLV_EDV);
    for j = 1:length(all_data{i}.ncs_EDV)
        % method 2 of EDV Polarity: is NCS at min LV cath location higher
        % or lower than start point
        ncs_EDV_polarity(j) = sign(all_data{i}.ncs_EDV{j}(l) - all_data{i}.ncs_EDV{j}(1));
        % ncs_EDV_polarity(j) = mode(sign(diff(all_data{i}.ncs_EDV_PP_calib{j})));
    end

    % store relevant data
    all_data{inc}.ncs_snr = ncs_snr;
    all_data{inc}.ncs_lin = ncs_lin;

    all_data{inc}.bioLV_EDV_mldiff = bioLV_EDV_mldiff;
    all_data{inc}.bioS1_EDV_mldiff = bioS1_EDV_mldiff;

    all_data{inc}.bioLV_EDV_PPdiff = bioLV_EDV_PPdiff;
    all_data{inc}.bioS1_EDV_PPdiff = bioS1_EDV_PPdiff;
    all_data{inc}.ncs_EDV_PPdiff = ncs_EDV_PPdiff;

    all_data{inc}.ncs_EDV_polarity = ncs_EDV_polarity;    
    inc = inc+1;
    clearvars -except all* inc
end
clearvars -except all*

% Analyze beat-wise polarity
clearvars -except all*
inc = 1;
RTwidth = 1;
early_polarity_width = 0.2;
HB_polarity = []; HB_polarity_early = [];
for i = 1:length(all_data)
    % Polarity of each ch's heartbeats
    for ch = 1:length(all_data{i}.ncs_EDV)
        polarity = zeros(size(all_data{i}.t_waveforms{ch}));
        for j = 1:length(all_data{i}.t_waveforms{ch})
            polarity(j) = waveformPolarity(all_data{i}.ncs_waveforms{ch}(:,j), ...
                all_data{i}.T_Wave_Idxs{ch}(j), RTwidth, 'mode');
        end
        HB_polarity(ch) = mode(polarity);    
        all_data{i}.polarity{ch} = polarity;

        % extract HB polarity of first early_polarity_width% of routine
        t_min = min(all_data{i}.t_waveforms{ch});
        t_max = max(all_data{i}.t_waveforms{ch});
        d = (t_max - t_min)*early_polarity_width;
        idxT = all_data{i}.t_waveforms{ch} >= t_min & all_data{i}.t_waveforms{ch} < t_min+d;
        HB_polarity_early(ch) = mode(polarity(idxT));
    end

    all_data{i}.ncs_HB_polarity = HB_polarity;
    all_data{i}.ncs_HB_polarity_early = HB_polarity_early;

end
clearvars -except all*

% Gather data for polarity/EDV analysis
clearvars -except all*
inc = 1;
for i = 1:length(all_data)
    data_ncs_EDV_calib{inc} = all_data{i}.ncs_EDV_PP_calib;
    data_ncs_EDV_PPdiff(inc, :) = all_data{i}.ncs_EDV_PPdiff;
    data_ncs_snr(inc, :) = all_data{i}.ncs_snr;
    data_ncs_lin(inc, :) = all_data{i}.ncs_lin;
    data_HB_polarity(inc, :) = all_data{i}.ncs_HB_polarity_early;
    data_EDV_polarity(inc, :) = all_data{i}.ncs_EDV_polarity;
    data_win{inc} = all_data{i}.win_filt;
    data_t{inc} = all_data{i}.t_waveforms;
    data_tEDV{inc} = all_data{i}.t_EDV;
    data_rt{inc} = [all_data{i}.dt all_data{i}.rt];
    if contains(all_data{i}.dt, '03-05') && contains(all_data{i}.rt, 'sion2') % put S1 segments here
        data_bio_EDV_PP_calib{inc} = all_data{i}.bioS1_EDV_PP_calib;
        data_bio_EDV_PPdiff(inc)  = all_data{i}.bioS1_EDV_PPdiff;
    else
        data_bio_EDV_PP_calib{inc} = all_data{i}.bioLV_EDV_PP_calib;
        data_bio_EDV_PPdiff(inc)  = all_data{i}.bioLV_EDV_PPdiff;
    end
    inc = inc+1;
end
clearvars -except all* data*

% LOO CV SNR weighted avg labelling whether systole is in same direction as
% EDV change
interventions = 1:length(all_data); n = max(interventions)-1;
CV_combinations = nchoosek(interventions, n);

% define parameter search
parameter_range = 0:0.2:2;
[x1, x2] = ndgrid(parameter_range, parameter_range);
grid_search = [x1(:), x2(:)]';
idx_grid = logical(1:size(grid_search, 2));
for j = 1:size(grid_search, 2)
    if grid_search(1, j) > grid_search(2, j) ...
            || sum(grid_search(:, j)) == 0
        idx_grid(j) = false;
    end
end
grid_search = grid_search(:, idx_grid);

% find best scaling parameters for polarity
optimal_parameter = [];
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    mean_acc_CV = [];
    for j = 1:size(grid_search, 2) % search over possible hyperparameters
        acc_CV = [];
        for k = 1:size(CV_combinations, 2) % loop over each index in CV combinations with a set of hyperparameters 
            ncs_snr = data_ncs_snr(CV_combinations(i, k), :); 
            ncs_HB_pol = data_HB_polarity(CV_combinations(i, k), :);
            ncs_EDV_pol = data_EDV_polarity(CV_combinations(i, k), :);
            
            ncs_polarity_agreement = ncs_HB_pol == ncs_EDV_pol;

            ncs_snr_weight = rescale(ncs_snr, grid_search(1, j), grid_search(2, j));

            ncs_pol_weight = zeros(length(ncs_polarity_agreement), 1);
            for ch = 1:length(ncs_HB_pol)
                ncs_pol_weight(ch) = ncs_polarity_agreement(ch) .* ncs_snr_weight(ch);
            end
            ncs_weight_mean = sum(ncs_pol_weight, 1);
            ncs_weight_mean = ncs_weight_mean / sum(ncs_snr_weight);

            if ncs_weight_mean > 0.5
                acc_CV(k) = 1;
            else
                acc_CV(k) = 0;
            end
        end
        mean_acc_CV(j) = mean(acc_CV); 
    end
    [m, l] = max(mean_acc_CV);
    optimal_parameter(:, i) = grid_search(:, l);
    optimal_acc(i) = m;
end

% apply scaling to determine if its a volume increase or decrease
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    % Find train set mean error with optimal parameters 
    acc_CV = []; all_train_weighted_scores = [];
    for k = 1:size(CV_combinations, 2) % loop over each index in CV combinations with optimal hyperparameters
        ncs_snr = data_ncs_snr(CV_combinations(i, k), :); 
        ncs_snr(isinf(ncs_snr)) = 0;
        ncs_HB_pol = data_HB_polarity(CV_combinations(i, k), :);
        ncs_EDV_pol = data_EDV_polarity(CV_combinations(i, k), :);
        
        ncs_polarity_agreement = ncs_HB_pol == ncs_EDV_pol;

        ncs_snr_weight = rescale(ncs_snr, optimal_parameter(1, i), optimal_parameter(2, i));

        ncs_pol_weight = zeros(length(ncs_polarity_agreement), 1);
        for ch = 1:length(ncs_HB_pol)
            ncs_pol_weight(ch) = ncs_polarity_agreement(ch) .* ncs_snr_weight(ch);
        end
        ncs_weight_mean = sum(ncs_pol_weight, 1);
        ncs_weight_mean = ncs_weight_mean / sum(ncs_snr_weight);

        if ncs_weight_mean > 0.5
            acc_CV(k) = 1;
        else
            acc_CV(k) = 0;
        end
        all_train_weighted_scores(k) = ncs_weight_mean;
    end
    train_err(i) = mean(acc_CV); 
    
    % Test LOO mean error with optimal parameters
    LOOSub = setdiff(interventions, CV_combinations(i,:));
    ncs_snr = data_ncs_snr(LOOSub, :); 
    ncs_snr(isinf(ncs_snr)) = 0;
    ncs_HB_pol = data_HB_polarity(LOOSub, :);
    ncs_EDV_pol = data_EDV_polarity(LOOSub, :);
    
    ncs_polarity_agreement = ncs_HB_pol == ncs_EDV_pol;

    ncs_snr_weight = rescale(ncs_snr, optimal_parameter(1, i), optimal_parameter(2, i));

    ncs_pol_weight = zeros(length(ncs_polarity_agreement), 1);
    for ch = 1:length(ncs_HB_pol)
        ncs_pol_weight(ch) = ncs_polarity_agreement(ch) .* ncs_snr_weight(ch);
    end
    ncs_weight_mean = sum(ncs_pol_weight, 1);
    ncs_weight_mean = ncs_weight_mean / sum(ncs_snr_weight);
    % deterimine if, after weighting, the vol seems to be going up or down
    % optimizing for total correct in CV fold
    if ncs_weight_mean > 0.5
        acc_CV = 1;
    else
        acc_CV = 0;
    end
    vol_dec(i) = acc_CV;
end

clearvars -except all* data* optimal_parameter CV_combinations interventions vol_dec

% Normalizing NCS volume decrease direction based on systole weighted avg direction w optimal parameters
data_ncs_EDV_normdir = cell(size(data_ncs_EDV_calib));
for i = 1:length(vol_dec)
    if vol_dec(i)
        % this means volume is decreasing: make sure EDV direction is neg
        for j = 1:length(data_ncs_EDV_calib{i})
            if data_EDV_polarity(i, j) == 1
                data_ncs_EDV_normdir{i}(j, :) = -data_ncs_EDV_calib{i}{j};
            else
                data_ncs_EDV_normdir{i}(j, :) = data_ncs_EDV_calib{i}{j};
            end
        end
    else
        % this means volume is increasing: make sure EDV direction is pos
        for j = 1:length(data_ncs_EDV_calib{i})
            if data_EDV_polarity(i, j) == 1
                data_ncs_EDV_normdir{i}(j, :) = data_ncs_EDV_calib{i}{j};
            else
                data_ncs_EDV_normdir{i}(j, :) = -data_ncs_EDV_calib{i}{j};
            end
        end

    end
end
clearvars -except all* data* vol_dec

% Use CV to find optimal parameters for scaling EDV (with channel selection)
interventions = 1:length(all_data); n = max(interventions)-1;
CV_combinations = nchoosek(interventions, n);

% Define channel selection
channel_selections = 1:8;

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
optimal_parameter = [];

for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    mean_mae_end_CV = [];
    for j = 1:size(grid_search, 2) % search over possible hyperparameters
        mae_end_CV = [];
        for k = 1:size(CV_combinations, 2) % loop over each index in CV combinations with a set of hyperparameters 
            ncs_lin = data_ncs_lin(CV_combinations(i, k), channel_selections); % Apply channel selection
            ncs_snr = data_ncs_snr(CV_combinations(i, k), channel_selections); % Apply channel selection
            ncs_snr(ncs_snr == 0) = min(ncs_snr(ncs_snr ~= 0));
            ncs_EDV_calib = data_ncs_EDV_normdir{CV_combinations(i, k)}(channel_selections,:); % Apply channel selection

            bio_EDV_calib = data_bio_EDV_PP_calib{CV_combinations(i, k)};
            
            ncs_lin_rescale = rescale(ncs_lin, grid_search(1, j), grid_search(2, j));
            ncs_snr_rescale = rescale(ncs_snr, grid_search(3, j), grid_search(4, j));
            ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;
            ncs_EDV_weight = zeros(length(ncs_snr), size(ncs_EDV_calib, 2)); 
            for ch = 1:size(ncs_EDV_calib, 1) 
                ncs_EDV_weight(ch, :) = ncs_EDV_calib(ch, :) .* ncs_linsnr_weight(ch);
            end
            ncs_EDV_weight_mean = sum(ncs_EDV_weight, 1);
            ncs_EDV_weight_mean = ncs_EDV_weight_mean / sum(ncs_linsnr_weight);

            mae_end_CV(k) = similarity(ncs_EDV_weight_mean(round(0.9*length(ncs_EDV_weight_mean)):end), ...
                bio_EDV_calib(round(0.9*length(bio_EDV_calib)):end), 'mae');
        end
        mean_mae_end_CV(j) = mean(mae_end_CV); 
    end
    [m, l] = min(mean_mae_end_CV);
    optimal_parameter(:, i) = grid_search(:, l);
end
clearvars -except all* data* optimal_parameter CV_combinations interventions channel_selections


% Find EDV weighted mean trace with optimal parameters (with channel selection)
train_err = []; test_err = [];
data_ncs_EDV_weight_mean = cell(size(data_ncs_EDV_calib));
for i = 1:size(CV_combinations, 1) % loop over all CV combinations
    % Find train set mean error with optimal parameters 
    train_err_all = [];
    for k = 1:size(CV_combinations, 2) % loop over each index in CV combinations with optimal hyperparameters
        ncs_lin = data_ncs_lin(CV_combinations(i, k), channel_selections); % Apply channel selection
        ncs_snr = data_ncs_snr(CV_combinations(i, k), channel_selections); % Apply channel selection
        % for those set to 0, make equal to the min non-zero SNR
        ncs_snr(ncs_snr == 0) = min(ncs_snr(ncs_snr ~= 0));
        ncs_EDV_calib = data_ncs_EDV_normdir{CV_combinations(i, k)}(channel_selections,:); % Apply channel selection
        bio_EDV_calib = data_bio_EDV_PP_calib{CV_combinations(i, k)};

        ncs_lin_rescale = rescale(ncs_lin, optimal_parameter(1, i), optimal_parameter(2, i));
        ncs_snr_rescale = rescale(ncs_snr, optimal_parameter(3, i), optimal_parameter(4, i));
        ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;
        ncs_EDV_weight = zeros(length(ncs_snr), size(ncs_EDV_calib, 2)); %length of ncs_snr is now the selected channels
        for ch = 1:size(ncs_EDV_calib, 1) %size of ncs_EDV_calib is now the selected channels
            ncs_EDV_weight(ch, :) = ncs_EDV_calib(ch, :) .* ncs_linsnr_weight(ch);
        end
        ncs_EDV_weight_mean = sum(ncs_EDV_weight, 1);
        ncs_EDV_weight_mean = ncs_EDV_weight_mean / sum(ncs_linsnr_weight);
        train_err_all(k) = similarity(ncs_EDV_weight_mean(round(0.9*length(ncs_EDV_weight_mean)):end), ...
            bio_EDV_calib(round(0.9*length(bio_EDV_calib)):end), 'merr');
    end
    train_err(i, :) = train_err_all; 
    % Test LOO mean error with optimal parameters
    LOOSub = setdiff(interventions, CV_combinations(i,:));
    ncs_lin = data_ncs_lin(LOOSub, channel_selections); % Apply channel selection
    ncs_snr = data_ncs_snr(LOOSub, channel_selections); % Apply channel selection
    ncs_snr(ncs_snr == 0) = min(ncs_snr(ncs_snr ~= 0));
    ncs_EDV_calib = data_ncs_EDV_normdir{LOOSub}(channel_selections,:); % Apply channel selection
    bio_EDV_calib = data_bio_EDV_PP_calib{LOOSub};
    ncs_lin_rescale = rescale(ncs_lin, optimal_parameter(1, i), optimal_parameter(2, i));
    ncs_snr_rescale = rescale(ncs_snr, optimal_parameter(3, i), optimal_parameter(4, i));
    ncs_linsnr_weight = ncs_lin_rescale + ncs_snr_rescale;
    ncs_EDV_weight = zeros(length(ncs_snr), size(ncs_EDV_calib, 2)); %length of ncs_snr is now the selected channels
    for ch = 1:size(ncs_EDV_calib, 1) %size of ncs_EDV_calib is now the selected channels
        ncs_EDV_weight(ch, :) = ncs_EDV_calib(ch, :) .* ncs_linsnr_weight(ch);
    end
    ncs_EDV_weight_mean = sum(ncs_EDV_weight, 1);
    ncs_EDV_weight_mean = ncs_EDV_weight_mean / sum(ncs_linsnr_weight);
    test_err(i) = similarity(ncs_EDV_weight_mean(round(0.9*length(ncs_EDV_weight_mean)):end), ...
        bio_EDV_calib(round(0.9*length(bio_EDV_calib)):end), 'merr');
    data_ncs_EDV_weight_mean{LOOSub} = ncs_EDV_weight_mean;
end

% Histogram of EDV Correlation: Pearson's Correlation Coefficient
% Initialize correlation array
corr_ = [];
rfNorm = cell(size(data_ncs_EDV_weight_mean));
bioNorm = cell(size(data_bio_EDV_PP_calib));

% Calculate normalized correlations
for i = 1:length(data_ncs_EDV_weight_mean)
    rfNorm{i} = rescale(data_ncs_EDV_weight_mean{i}, -1, 0);
    bioNorm{i} = rescale(data_bio_EDV_PP_calib{i}, -1, 0);

    % Compute Pearson's correlation coefficient
    corr_(i) = corr(rfNorm{i}', bioNorm{i}');
end

% Create a histogram of the correlation coefficients
be = 0.9:0.01:1;
fn = 'Arial'; fs = 20;
figure('Position', [100, 100, 600, 350]);
histogram(corr_, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k',...
    'BinEdges', be, 'Normalization', 'probability');
xlabel("\rho_E_D_V_N_o_r_m", 'FontSize', fs+4);
ylabel('Probability');
yticks(0:0.1:0.5);
xlim([0.9 1]);
grid on;
set(gca, 'FontName', fn, 'FontSize', fs);

meanCorr = mean(corr_); stdCorr = std(corr_);