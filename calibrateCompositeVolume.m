function powerlab = calibrateCompositeVolume(powerlab, fileDate)
% Calibrate, the LV composite volume signal based on the scaling factors found in power lab CO Calib B2 
% only processing the appropriate recordings

% Determined that 3-6 and 3-7 need scaling
% 3-7: scaling by 4*, subtract 80 from volume
% 3-6 scaling by 2*, subtract 50 from volume


if strcmp(fileDate, '03-06-24')
    powerlab.data{3} = powerlab.data{3}*2-50;
elseif strcmp(fileDate, '03-07-24')
    powerlab.data{3} = powerlab.data{3}*4-80;
end



end