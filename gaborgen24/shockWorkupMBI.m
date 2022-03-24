function shockWorkupMBI()
    %% header
    nrPulses = 8;
    secBtwPulses = 1/120;
    
    %% initialize Biopac device & output port for electrical stimulation
    % in case it was already initialized 
    daqreset; 
    % create object for data acquisition toolbox device, add output channel, 
    % set scan rate, and set DC to baseline level
    d = daq('mcc');
    addoutput(d,'Board0','Ao0','Voltage');
    d.Rate = 960000;
    write(d, -10);
       
    % create square function for shock delivery
    biopacTrigger = repmat([10.*ones(1,round(secBtwPulses*d.Rate/2)), ... 
                           -10.*ones(1,round(secBtwPulses*d.Rate/2))]', ...
                           nrPulses, 1);
    
    %% apply shocks or end routine
    quitRoutine = false;
    disp('press [p] to send pulse or [q] to quit routine');
    while quitRoutine == false
        [~, ~, keyCode, ~] = KbCheck(-1);
        if keyCode(KbName('p'))
            disp('Sending Pulse');
            write(d, biopacTrigger);
            WaitSecs(.200);
            disp('press [p] to send pulse or [q] to quit routine');
        elseif keyCode(KbName('q'))
            disp('Quit Routine');
            quitRoutine = true;
        end
        WaitSecs(.010);
    end
end

