function polarity = waveformPolarity(waveform, idxT, width, method)
% returns the sign of the polarity of ventricular contraction to determine
% if EDV change is more or less volume

% method 1: mode of slope change in window from R-T wave
switch method
    case 'mode'
        polarity = mode(sign(diff(waveform(1:round(idxT*width)))));
    otherwise
        error('unknown pol extraction method')
end



end