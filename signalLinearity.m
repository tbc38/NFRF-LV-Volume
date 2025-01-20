function [SNR_sig_1st_2nd_ratio] = signalLinearity(ncs_sig,fRange,Fs)
% Objective function 2 fundamental tone to 2nd harmonic ratio
    L_sig = length(ncs_sig);             % Length of signal
    Y_sig = fft(detrend(ncs_sig));
    P2 = abs(Y_sig/L_sig); 
    P1 = P2(1:round(L_sig/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L_sig/2))/L_sig;
    % F
    mag_1_peak = max(P1(fRange(2)>f&f>fRange(1)));
    mag_2_peak = max(P1(fRange(2)*2>f&f>fRange(1)*2));
    SNR_sig_1st_2nd_ratio = 10*log10(mag_1_peak /mag_2_peak );

end