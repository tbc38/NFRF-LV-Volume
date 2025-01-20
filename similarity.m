function simScore = similarity(sig1, sig2, simType)

switch simType
    case 'dot'
        simScore = dot(sig1, sig2)/length(sig1);
    case 'normxcorr'
        simScore = xcorr(sig1, sig2, 0, 'normalized');
    case 'mse'
        simScore = mean((sig1-sig2).^2);
    case 'mae'
        simScore = mean(abs(sig1-sig2));
    case 'merr'
        simScore = mean((sig1-sig2));
    otherwise
        error('invalid simType')
end




end