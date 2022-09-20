Screen('Preference', 'SkipSyncTests', 0);
rand('state',sum(100*clock));

%%
backgroundCol = 0;

% US parameters
soundSamp = 44100; % sampling rate of sound presentation
soundDur = 1; % duration of noise stimulus in s
nrNoises = 5;
noiseJitter = [4 6];


%% Loud noise
soundIntensityMax = 1; % 92 dB
soundIntensityHigh = 10^(log10(soundIntensityMax) - (1/20)); % 91 dB
soundIntensityLow = 10^(log10(soundIntensityMax) - (4/20)); % 88 dB

whiteNoiseHigh = rand(1, soundSamp*soundDur) .* soundIntensityHigh;
whiteNoiseLow = rand(1, soundSamp*soundDur) .* soundIntensityLow;



%% get all info for screen & window

screens = Screen('Screens');
screenNumber = max(screens);
w = Screen('OpenWindow', screenNumber, backgroundCol);
[width,height] = Screen('WindowSize', w);
Screen('TextSize', w, 32); 


%% coordinates
yesCoords = [.2*width, .8*height, .3*width, .9*height];
noCoords = [.7*width, .8*height, .8*width, .9*height];
buttonCoords = [yesCoords', noCoords'];

%%
InitializePsychSound;


%%
quitRoutine = false;
routineIndex = 1;
validClick = false;

while quitRoutine == false
    if routineIndex == 1
        message = ['In a moment, you will hear the loud noise that we would like\n' ...
                   'to present to you multiple times in this experiment. Ideally,\n' ...
                   'you find this noise unpleasant but can tolerate it.\n\n' ...
                   'Please click the mouse button to play the sound.'];
        DrawFormattedText(w, message, 'center', 'center', 255);
        Screen('Flip', w);

        waitForClick();
        DrawFormattedText(w, '+', 'center', 'center', 255);
        Screen('Flip', w);
        
        WaitSecs(2);
        Snd('Play', whiteNoiseHigh, soundSamp, 16);
        WaitSecs(soundDur + 1);
        
        message = ['Do you think that you are okay with hearing this noise\n' ...
                   'multiple times today?\n\n' ...
                   'Please click [YES] or [NO].'];
        DrawFormattedText(w, message, 'center', 'center', 255);
        
        Screen('FillRect', w, [127.5, 127.5, 127.5], buttonCoords);
        DrawFormattedText(w, 'YES', 'center', 'center', 0, [], [], [], [], [], yesCoords);
        DrawFormattedText(w, 'NO', 'center', 'center', 0, [], [], [], [], [], noCoords);
        
        Screen('Flip', w);
        
        while validClick == false
            [x,y] = waitForClick();
            if x > yesCoords(1) && x < yesCoords(3) && y > yesCoords(2) && y < yesCoords(4)
                validClick = true;
                routineIndex = 2;
            elseif x > noCoords(1) && x < noCoords(3) && y > noCoords(2) && y < noCoords(4)
                validClick = true;
                routineIndex = 3;
            end
        end
        validClick = false;
        Screen('Flip', w);
        WaitSecs(1);
    end
    
    
    
    if routineIndex == 2
        message = ['Great! Just to be sure, we will play this sound a couple of\n' ...
                   'times in a row. We want to make sure that you can tolerate\n' ...
                   'the noise even if you listen to it multiple times.\n\n' ...
                   'Please click the mouse button to play the sounds.'];
        DrawFormattedText(w, message, 'center', 'center', 255);
        Screen('Flip', w);

        waitForClick();
        DrawFormattedText(w, '+', 'center', 'center');
        Screen('Flip', w);
        
        for noiseInd = 1:nrNoises
            WaitSecs(randi(noiseJitter));
            Snd('Play', whiteNoiseHigh, soundSamp, 16);
            WaitSecs(soundDur);
        end
        
        message = ['Are you okay with hearing this noise multiple\n' ...
                   'times today?\n\n' ...
                   'Please click [YES] or [NO].'];
        DrawFormattedText(w, message, 'center', 'center', 255);

        Screen('FillRect', w, [127.5, 127.5, 127.5], buttonCoords);
        DrawFormattedText(w, 'YES', 'center', 'center', 0, [], [], [], [], [], yesCoords);
        DrawFormattedText(w, 'NO', 'center', 'center', 0, [], [], [], [], [], noCoords);
        
        Screen('Flip', w);
        
        while validClick == false
            [x,y] = waitForClick();
            if x > yesCoords(1) && x < yesCoords(3) && y > yesCoords(2) && y < yesCoords(4)
                validClick = true;
                quitRoutine = true;
                useDB = 1;
            elseif x > noCoords(1) && x < noCoords(3) && y > noCoords(2) && y < noCoords(4)
                validClick = true;
                routineIndex = 3;
            end
        end
        validClick = false;
        Screen('Flip', w);
        WaitSecs(1);
    end
    
    
    
    if routineIndex == 3
        message = ['No problem! We can use a noise with lower volume. \n' ...
                   'It still should be unpleasant to you but you should be\n' ...
                   'able to tolerate it.\n\n' ...
                   'Are you okay with listening to it once?\n\n' ...
                   'Please click [YES] or [NO].'];
        DrawFormattedText(w, message, 'center', 'center', 255);

        Screen('FillRect', w, [127.5, 127.5, 127.5], buttonCoords);
        DrawFormattedText(w, 'YES', 'center', 'center', 0, [], [], [], [], [], yesCoords);
        DrawFormattedText(w, 'NO', 'center', 'center', 0, [], [], [], [], [], noCoords);
        
        Screen('Flip', w);
        
        while validClick == false
            [x,y] = waitForClick();
            if x > yesCoords(1) && x < yesCoords(3) && y > yesCoords(2) && y < yesCoords(4)
                message = 'yes';
                validClick = true;
                routineIndex = 4;
            elseif x > noCoords(1) && x < noCoords(3) && y > noCoords(2) && y < noCoords(4)
                message = 'no';
                validClick = true;
                quitRoutine = true;
                usedDB = 3;
            end
        end
        validClick = false;
        Screen('Flip', w);
        WaitSecs(1);
    end
    
    
    
    if routineIndex == 4
        message = ['Okay. You will hear the same noise as before once but with lower\n' ...
                   'volume. We want to present it to you multiple times in this experiment. \n' ...
                   'Ideally, you find this noise unpleasant but can tolerate it.\n\n' ...
                   'Please click the mouse button to play the sound.'];
        DrawFormattedText(w, message, 'center', 'center', 255);
        Screen('Flip', w);

        waitForClick();
        DrawFormattedText(w, '+', 'center', 'center', 255);
        Screen('Flip', w);
        
        WaitSecs(2);
        Snd('Play', whiteNoiseLow, soundSamp);
        WaitSecs(soundDur + 1);
        
        message = ['Do you think that you are okay with hearing this noise\n' ...
                   'multiple times today?\n\n' ...
                   'Please click [YES] or [NO].'];
        DrawFormattedText(w, message, 'center', 'center', 255);
        
        Screen('FillRect', w, [127.5, 127.5, 127.5], buttonCoords);
        DrawFormattedText(w, 'YES', 'center', 'center', 0, [], [], [], [], [], yesCoords);
        DrawFormattedText(w, 'NO', 'center', 'center', 0, [], [], [], [], [], noCoords);
        
        Screen('Flip', w);
        
        while validClick == false
            [x,y] = waitForClick();
            if x > yesCoords(1) && x < yesCoords(3) && y > yesCoords(2) && y < yesCoords(4)
                validClick = true;
                routineIndex = 5;
            elseif x > noCoords(1) && x < noCoords(3) && y > noCoords(2) && y < noCoords(4)
                validClick = true;
                quitRoutine = true;
                useDB = 3;
            end
        end
        validClick = false;
        Screen('Flip', w);
        WaitSecs(1);
    end
    
    
    
    if routineIndex == 5
        message = ['Great! Again, just to be sure, we will play this sound a couple of\n' ...
                   'times in a row. We want to make sure that you can tolerate\n' ...
                   'the noise even if you listen to it multiple times.\n\n' ...
                   'Please click the mouse button to play the sounds.'];
        DrawFormattedText(w, message, 'center', 'center', 255);
        Screen('Flip', w);

        waitForClick();
        DrawFormattedText(w, '+', 'center', 'center');
        Screen('Flip', w);
        
        for noiseInd = 1:nrNoises
            WaitSecs(randi(noiseJitter));
            Snd('Play', whiteNoiseLow, soundSamp);
            WaitSecs(soundDur);
        end
        
        message = ['Are you okay with hearing this noise multiple\n' ...
                   'times today?\n\n' ...
                   'Please click [YES] or [NO].'];
        DrawFormattedText(w, message, 'center', 'center', 255);

        Screen('FillRect', w, [127.5, 127.5, 127.5], buttonCoords);
        DrawFormattedText(w, 'YES', 'center', 'center', 0, [], [], [], [], [], yesCoords);
        DrawFormattedText(w, 'NO', 'center', 'center', 0, [], [], [], [], [], noCoords);
        
        Screen('Flip', w);
        
        while validClick == false
            [x,y] = waitForClick();
            if x > yesCoords(1) && x < yesCoords(3) && y > yesCoords(2) && y < yesCoords(4)
                validClick = true;
                quitRoutine = true;
                useDB = 2;
            elseif x > noCoords(1) && x < noCoords(3) && y > noCoords(2) && y < noCoords(4)
                validClick = true;
                quitRoutine = true;
                useDB = 3;
            end
        end
        validClick = false;
        Screen('Flip', w);
        WaitSecs(1);
    end    
end


if useDB == 1
    message = ['Thank you! You stated that you are okay with listening to\n' ...
               'the noise in the upcoming experiment.\n\n' ...
               'The experimenter will be with you in a moment.'];
elseif useDB == 2
    message = ['Thank you! You stated that you are okay with listening to\n' ...
               'the noise at a lower volume in the upcoming experiment.\n\n' ...
               'The experimenter will be with you in a moment.'];
elseif useDB == 3
    message = ['Thank you! You stated that you are not okay with listening to\n' ...
               'the noise in the upcoming experiment.\n\n' ...
               'The experimenter will be with you in a moment.'];
end

DrawFormattedText(w, message, 'center', 'center', 255);
Screen('Flip', w);
KbStrokeWait();

%%
Screen('CloseAll');




function [x, y, buttons] = waitForClick
% wait for any mouse click, detected via the GetMouse function in
% Psychtoolbox
    buttons = 0;
    while ~any(buttons) % wait for press
        [x, y, buttons] = GetMouse;
    
        % Wait 10 ms before checking the mouse again to prevent
        % overload of the machine at elevated Priority()
        WaitSecs(0.01);
    end
end