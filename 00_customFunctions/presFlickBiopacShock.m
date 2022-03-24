% presents a flickering image texture, sends a shock trigger via USB
% port output and a stimulus onset marker via parallel port. 
% flickVec = a one-dimensional vector with opaqueness values from 0 to
% 1 for each frame
% imageSize: can be empty (no scaling of original image) or have two values
% for size in pixels ([x, y])
% texture = texture object
% outputSignal = analog signal that goes into BIOPAC Stimulator
% daqObject = DataAcquistion (toolbox) device
% lptObject = parallel port object created with io64
% lptAdresses = adresse(s) of LPT port(s) as DEC
% lptCodes = DEC codes to send via parallel port
% function returns actual flicker stimulus duration and stimulus onset time
% in seconds (system time)
function [actFlickDur, startTime] = presFlickBiopacShock(window, flickVec, imageSize, texture, ...
                                                         outputSignal, daqObject, ...
                                                         lptObject, lptAdresses, lptCodes)
    if isempty(imageSize)
        imageCoordinates = [];
    else
        [winSize(1), winSize(2)] = Screen('WindowSize', window);
        winCenter = winSize ./ 2;
        imageCoordinates = [winCenter(1)-imageSize(1)/2 winCenter(2)-imageSize(2)/2, ...
                            winCenter(1)+imageSize(1)/2 winCenter(2)+imageSize(2)/2];
    end

    % send marker via LPT port(s) and log time for measuring flicker duration
    for portI = 1:size(lptAdresses,1)
        io64(lptObject, lptAdresses(portI,:), lptCodes(portI));
    end

    % send analog output in the background
    start(daqObject);
    write(daqObject, outputSignal);

    startTime = GetSecs();

    for frame = 1:length(flickVec)
        Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], flickVec(frame));
        Screen('Flip', window);
    end

    % check time passed since before flicker stim
    actFlickDur = GetSecs() - startTime;

    % stop output channel
    stop(daqObject);
end