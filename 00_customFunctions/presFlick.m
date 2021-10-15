function actFlickDur = presFlick(window, condI, flickVec, imageSize, textureVec)
    if isempty(imageSize)
        imageCoordinates = [];
    else
        [winSize(1), winSize(2)] = Screen('WindowSize', window);
        winCenter = winSize ./ 2;
        imageCoordinates = [winCenter(1)-imageSize(1)/2 winCenter(2)-imageSize(2)/2 winCenter(1)+imageSize(1)/2 winCenter(2)+imageSize(2)/2];
    end
    
    % send marker via IO Port and log time for measuring flicker duration
    %IOPort('Write', s3, 'fufufufu99fufufu');
    %[nwritten, when, errmsg, prewritetime, postwritetime, lastchecktime] = IOPort('Write', s3, 'fufufufu99fufufu');
    tic;
    
    % Present Flicker Stimulus
    for frame = 1:size(flickVec, 3)
        Screen('DrawTexture', window, textureVec(condI(3)), [], imageCoordinates, 0, [], flickVec(condI(1), condI(2), frame));
        Screen('Flip', window);
    end
    
    % check time passed since before flicker stim
    actFlickDur = toc;
end