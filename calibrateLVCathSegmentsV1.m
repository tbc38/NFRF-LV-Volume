function power_sync_calib = calibrateLVCathSegmentsV1(power_sync)
% function takes in the powerlab, and calibrates the secondary segments
% such that they match the calbration of the composite volume

power_sync_calib = power_sync;

calib_vol = power_sync.data{3};
S1_uncalib = power_sync.data{4};

meds = ones(size(power_sync.data));
for i = 1:length(power_sync.data)
    meds(i) = median(power_sync.data{i});
end

% looking for and removing zero-ed data due to powerlab gaps
idx_gaps = false(size(calib_vol));
for i = 1:length(calib_vol)
    if calib_vol(i) == 0 && S1_uncalib(i) == 0
        idx_gaps(i) = true;
        calib_vol(i) = meds(3);
        for seg = 4:8
            power_sync.data{seg}(i) = meds(seg);
        end
    end
end

for seg = 4:8
    foo = rescale(power_sync.data{seg}', 0, 100);

    % solve least squares for: Ax=b
    A = [ones(length(foo), 1) foo];
    b = calib_vol';
    x = A\b;
    power_sync_calib.data{seg} = double(foo*x(2)+x(1))';
end

% re-inserting zeros in gaps
calib_vol(idx_gaps) = 0;
for seg = 4:8
    power_sync_calib.data{seg}(idx_gaps) = 0;
end

% figure;
% plot(power_sync_calib.data{3}, 'k', 'LineWidth', 1); hold on;
% plot(power_sync_calib.data{4}, 'g', 'LineWidth', 1); 
% plot(power_sync_calib.data{5}, 'm', 'LineWidth', 1); 
% plot(power_sync_calib.data{6}, 'c', 'LineWidth', 1);
% plot(power_sync_calib.data{7}, 'b', 'LineWidth', 1); 
% plot(power_sync_calib.data{8}, 'r', 'LineWidth', 1); 
% legend('composite', 'S1', 'S2', 'S3', 'S4', 'S5');

end