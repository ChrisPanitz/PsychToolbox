% Written to wait for and detect MRI triggers sent via serial port; uses
% keyboard logging commands. Logs system time when trigger is detected.
% keys = vector of acceptable keys as chars (e.g. ['a', 'b']) or ints ([4, 5])
function startTime = waitForScanTriggerKb(keys)
    % create keyCode object with 256 zeros (length of KbCheck keyCode output)
    keyCode = zeros(1,256);
    
    if ischar(keys)
        % translate key strings into DEC format
        keyInd = NaN(length(keys),1);
        for key = 1:length(keys)
            keyInd(key) = KbName(keys(key));
        end
    else
        % leave integers as is
        keyInd = keys;
    end
    
    % check for "key strokes" and log them in keyCode
    while sum(keyCode(keyInd)) == 0
        [~, ~, keyCode] = KbCheck;
        WaitSecs(0.001);
    end

    % When key stroke (aka trigger) has been detected, log system time
    startTime = GetSecs;

end % function