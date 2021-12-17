% Written to wait for and detect MRI pulses sent via serial port; uses
% keyboard logging commands. Logs system time when pulse is detected.
% keys = vector of acceptable keys as string (e.g. ['a', 'b'])

function startTime = waitForScanTriggerKb(keys)
    % create keyCode object with 256 zeros (length of KbCheck keyCode output)
    keyCode = zeros(1,256);
    
    % translate key strings into DEC format
    keyInd = NaN(length(keys),1);
    for key = 1:length(keys)
        keyInd(key) = KbName(keys(key));
    end
    
    % check for "key strokes" and log them in keyCode
    while sum(keyCode(keyInd)) == 0
        [~, ~, keyCode] = KbCheck;
        WaitSecs(0.001);
    end

    % When key stroke (aka pulse) has been detected, log system time
    startTime = GetSecs;

end % function