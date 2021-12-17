% presents a flickering image texture, sends a shock trigger via serial
% port output and a stimulus onset marker via parallel port. 
% flickVec = a one-dimensional vector with opaqueness values from 0 to
% 1 for each frame
% imageSize: can be empty (no scaling of original image) or have two values
% for size in pixels ([x, y])
% texture = texture object
% shockOnsetFrame & shockOnsetDur = onset and duration of trigger in frames
% lptObject = parallel port object created with io64
% lptAdresses = adresse(s) of LPT port(s) as DEC
% lptCodes = DEC codes to send via parallel port
% serialObject = serial port object created with serial, opened with fopen
% serialCode = string sent via serial port
% funciton returns actual flicker stimulus duration and stimulus onset time
% in seconds (system time)
function [actFlickDur, timeOnset] = presFlickShockMRI(window, flickVec, imageSize, texture, shockOnsetFrame, shockDurFrame, lptObject, lptAdresses, lptCodes, serialObject, serialCode)
    if isempty(imageSize)
        imageCoordinates = [];
    else
        [winSize(1), winSize(2)] = Screen('WindowSize', window);
        winCenter = winSize ./ 2;
        imageCoordinates = [winCenter(1)-imageSize(1)/2 winCenter(2)-imageSize(2)/2 winCenter(1)+imageSize(1)/2 winCenter(2)+imageSize(2)/2];
    end
    
    % send marker via LPT port(s) and log time for measuring flicker duration
    for portI = 1:size(lptAdresses,1)
        %io64(lptObject, lptAdresses(portI,:), lptCodes(portI));
    end
    timeOnset = GetSecs();
    
    % Present Flicker Stimulus
    for frame = 1:length(flickVec)
        Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], flickVec(frame));
        %Screen('Flip', window); % PUT IN AGAIN
        %%% FIND OUT HOW TO TRIGGER BIOPAC SHOCKS
        if (frame >= shockOnsetFrame) && (frame < shockOnsetFrame + shockDurFrame)
            %IOPort('Write', serialObject, serialCode);
            Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], [], [255 0 0]);
        end
        Screen('Flip', window); % TAKE OUT AGAIN
    end

    % check time passed since before flicker stim
    timeEnd = GetSecs();
    actFlickDur = timeEnd - timeOnset;
end