function [EDV, ESV] = interpretVolumesV2(waveforms)
% Return the peak and min following the peak of the volume trace: V2: run on waveforms

len = size(waveforms, 1);

search_width = round(len*0.75);

for i = 1:size(waveforms,2)
        [EDV(i), ~] = max(waveforms(1:search_width, i));

        [ESV(i), ~] = min(waveforms(1:search_width, i));
end

end