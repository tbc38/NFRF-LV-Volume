function [PPs, EDVs, ESVs, Polaritys] = P2P_V6_1(waveforms, T, segLen)
% will return the peak-to-peaks of on a per-waveform basis, using T-wave as marker
% generalizable to complex waveforms
% will find the max difference betwen the beginning search width and the search width around the
% T-wave
% V6: will also return the EDV/ESV, based on same locations as the PP
% V6.1: change the polarity estimate to be the mode of the R-T segment
% slope
PPs = zeros(1, size(waveforms, 2));
EDVs = zeros(1, size(waveforms, 2));
ESVs = zeros(1, size(waveforms, 2));
Polaritys = zeros(1, size(waveforms, 2));
searchWin = round(0.05*segLen);
for i = 1:size(waveforms, 2)
    sig = waveforms(:,i);
    r = 0;
    for s = 1:searchWin
        for e = -round(searchWin/2):round(searchWin/2)
            r_test = abs(sig(s) - sig(T(i)+e));
            if r_test > r
                r = r_test;
                s_best = s; e_best = e;
            end
        end
    end
    PPs(i) = r;

    % V6 change: if positive polarity, extract EDV as peak near
    % R-wave and ESV and min near T-wave. VV is negative pol
    % positive polarity
    if sig(s_best) > sig(T(i)+e_best)
        EDVs(i) = sig(s_best);
        ESVs(i) = sig(T(i)+e_best);
    else
        ESVs(i) = sig(s_best);
        EDVs(i) = sig(T(i)+e_best);
    end
    
    % V6.1 edit: polarity based on slope during systole
    pol = mode(sign(diff(sig(1:T(i)))));
    if pol > 0
        Polaritys(i) = 0;
    else
        Polaritys(i) = 1;
    end

end

end