function idxJumps = findSignalJumps(ncs, win, thr)
    if isempty(win); win = 20; end
    if isempty(thr); thr = 0.005; end
    % Detect signal jumps
    jumpUp = []; jumpDown = []; upInc = 1; downInc = 1;
    for i = 1:2:length(ncs)-win-1
        d = range(ncs(i:i+win, 1));
        if d > thr
            jumpUp(upInc) = i; upInc = upInc+1;
        elseif d < -thr
            jumpDown(downInc) = i; downInc = downInc+1;
        end
    end
    idxJumps = sort([jumpUp jumpDown]);
end