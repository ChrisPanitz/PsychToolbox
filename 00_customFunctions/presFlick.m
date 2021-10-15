% presents a flickering image texture; flickVec is a one-dimensional vector
% with opaqueness values from 0 to 1 for each frame; imageSize can be empty
% (no scaling of original image) or have two values for size in pixels ([x
% y]); texture is texture object; returns actual flicker stimulus duration
% in seconds obtained with tic & toc commands
function actFlickDur = presFlick(window, flickVec, imageSize, texture)
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
        Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], flickVec(frame));
        Screen('Flip', window);
    end
    
    % check time passed since before flicker stim
    actFlickDur = toc;
end