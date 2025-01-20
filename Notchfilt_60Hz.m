function sig_filt = Notchfilt_60Hz(sig, Q, fs, showPlot)
% Thomas Conroy | June '23 | Notch filter (simualted comb filter) for breathing removal
% this works to remove breathing pretty well
% V2 should adaptively find peaks and tune the bandwidth/Q factor for their width. Or can stop combs
% once we reach HR frequency. Research into adaptive filtering would help here

wo = 60/(fs/2); bw = wo/(Q); [n, d] = iirnotch(wo, bw); sig_filt = filtfilt(n, d, sig);


if showPlot
    figure; 
    subplot(211); plot(sig_filt); hold on; plot(sig); 
    subplot(212); 
    [p,f] = ncsFFT(sig_filt, fs, 0); plot(f, p); hold on;
    [p,f] = ncsFFT(sig, fs, 0); plot(f, p); xlim([0.1 65]);
end


end