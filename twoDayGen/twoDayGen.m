function twoDayGen(subNo, csPerm)
%% SET TO ZERO FOR ACTUAL EXPERIMENT
Screen('Preference', 'SkipSyncTests', 1);

%% Header
% trial parameters 
nrBlocksHab = 1; % # of blocks in acquisition; ratings after each block
nrTrialsHab = 10; % # of trials per stimulus and per block

nrBlocksAcq = 2; % # of blocks in acquisition; ratings after each block
nrTrialsAcq = 15; % # of trials per stimulus and per block
nrUSAcq = 15; % # of US per block in acquisition

nrBlocksExt = 1; % # of blocks in acquisition; ratings after each block
nrTrialsExt = 15; % # of trials per stimulus and per block

minMaxItiSec = [5.5 15.4]; % min & max ITI in sec, taken from exponential distribution
meanItiSec = 7; % mean ITI duration in sec

% CS parameters
rotStim = [15 35 55 75]; % vector of orientation angles (clockwise from 0 = vertical)
cycPerDeg = 3.5; % spatial frequency of Gabors in cycles/degree visual angle
imSizeAng = [5 5]; % size of Gabor patches in degrees visual angle

% distance participant <-> screen
seatDistInch = 125/2.54;

% flicker parameters
% framPerCyc for 120Hz screen: 8 = 15 Hz
framPerCyc = 8; % duration of one cycle in frames
flickMod = 'box'; % how flicker stimulus is modulated. 'box' or 'sin'
flickDurFrames = 248; % duration of whole stimulus in frames; 248/120 = 2.067 s

% US parameters
shockOnsetFrame = 241; % onset frame of US; 241 = 2.0083s
shockDurFrame = 8; % US duration in frames; 8 = .0667s

% ports & codes
%%% NOT TESTED YET %%%
addressesLPTHex = ['E020'; 'E030']; % LPT port adress in HEX
lptCodes = [2 2]; % codes (DEC) sent via LPT when Gabors are presented
TRdur = 2; % Time of Repetition in sec
TRtriggerCodes = ['t','5']; % serial port codes from fMRI, coded as keyboard strokes

% rating & instruction parameters
pauseBtwRat = 1*TRdur; % pause duration between two ratings, in sec
ratingDur = 2*TRdur; % time out for ratings, in sec
instructDur = 2*TRdur; % fixed duration of instruction screens, in sec
buttonsLeft = KbName('a'); % DEC key codes for moving rating cursor left
buttonsRight = KbName('d'); % DEC key codes for moving rating cursor right
buttonsOK = KbName('w'); % DEC key codes for confirm rating selection

% aethetics
backgroundCol = 127; % bakcground color; 127 = mid gray
gratingDark = 0; % dark grating stripes = 0 = black
gratingBright = 255; % bright grating stripes = 0 = white
fontSize = 32; % well... font size

% for logfiles and rating files
logFileFolder = '/Users/christianpanitz/Desktop/GitHub/PsychToolbox/twoDayGen/logfiles/';
logFilePrefix = 'twoDayGen'; % TODO: come up with something more specific
ratFileFolder = '/Users/christianpanitz/Desktop/GitHub/PsychToolbox/twoDayGen/ratings/';
ratFilePrefix = 'twoDayGen'; % TODO: come up with something more specific

% strings of instructions to participants          
welcomeMsg_expStart = ['Welcome and thank you for participating in our\n' ...
                       'experiment. The experimenter will start the experiment soon.'];
goodbyeMsg = ['You have completed the task. Thank you for your participation!\n\n' ...
              'The experimenter will be with you in a moment.'];

          

%% clearing and initializing stuff
IOPort('CloseAll');
clc;

% Setting rng seed; "old" rand function (instead of rng) was used for
% reasons of compatibility with lab computer
rand('state',sum(100*clock));

AssertOpenGL; 

%%% TODO: INCLUDE PORT STATEMENTS WHEN PORT AVAILABLE %%%
% Serial port
%serPort = serial('COM3'); % NEED TO CHANGE FOR SCANNER COMPUTER?
%fopen(serPort); % connects serial port object to the serial port

% LPT port
%lptPortObj = io64;
%status = io64(lptPortObj);
adressesLPT = hex2dec(addressesLPTHex);
   
%% Prepare logging

% files to write trial & rating data + variable names for header line
% logfiles
logFileName = strcat(logFileFolder, logFilePrefix, '_', num2str(subNo), '_logfile.dat');
fIDLog = fopen(logFileName, 'w');
fprintf(fIDLog,'partNo,csPerm,trial,condCode,timeSinceStartPTB,timeSinceFirstTR,actFlickDur,itiDur\n');
fclose(fIDLog);
% rating files
ratFileName = strcat(ratFileFolder, ratFilePrefix, '_', num2str(subNo), '_ratings.dat');
fIDRat = fopen(ratFileName, 'w');
fprintf(fIDRat,['val_csp,val_g1,val_g2,val_g3,ar_csp,ar_g1,ar_g,ar_g3,' ...
                'fear_csp,fear_g1,fear_g2,fear_g3,exp_csp,exp_g1,exp_g2,exp_g3\n']);
fclose(fIDRat);

%% create vectors & matrices for stimuli & trials
% permutation: CS+ is always one of the outermost orientations; orientations are shuffled
% thereby, CS+ is always stim #1, GS1 = stim #2, GS2 = stim #3, GS3 = stim #4
if csPerm == 1
    % nothing - leave rotStim as is
elseif csPerm == 2
    rotStim = rotStim(length(rotStim):-1:1); % reverse order of rotStim
else
    error('csPerm must be either 1 or 2'); 
end

% compute some "intermediate" variables to make following code easier to follow
nrCS = length(rotStim); % # of CS
conds = 1:nrCS; % vector counting up for # of CS
trialsPerBlockHab = nrCS*nrTrialsHab; % how many trials per HAB block across all CS
trialsPerBlockAcq = nrCS*nrTrialsAcq; % how many trials per ACQ block across all CS
trialsPerBlockExt = nrCS*nrTrialsExt; % how many trials per EXT block across all CS

% initialize (empty) trialMat to store CS type, US [yes/no] & ITI duration
trialMat = NaN(nrBlocksHab*trialsPerBlockHab + ...
               nrBlocksAcq*trialsPerBlockAcq + ...
               nrBlocksExt*trialsPerBlockExt, 3);
% initialize (empty) ratingAfterTrial that contains indices for every
% block's last trial
ratingAfterTrial = NaN(nrBlocksHab+nrBlocksAcq+nrBlocksExt,1);

% some helping count variables
trialCount = 0;
ratingCount = 0;

% for habituation trials
for bl = 1:nrBlocksHab
    tempMat = zeros(trialsPerBlockHab,3);
    tempMat(:,1) = sort(repmat(conds', nrTrialsHab, 1)); % CS type column
    tempMat(:,3) = createExpIntDurs(trialsPerBlockHab, minMaxItiSec, meanItiSec*trialsPerBlockHab); % ITI columns
    tempMat = tempMat(randperm(trialsPerBlockHab),:); % shuffle trial order
    
    % write values into trialMat and ratingAfterTrial
    trialCount = trialCount + trialsPerBlockHab;
    ratingCount = ratingCount + 1;
    trialMat((trialCount-trialsPerBlockHab+1):trialCount,:) = tempMat;
    ratingAfterTrial(ratingCount) = trialCount;
end

% for acquisition trials
for bl = 1:nrBlocksAcq
    tempMat = zeros(trialsPerBlockAcq,3);
    tempMat(:,1) = sort(repmat(conds', nrTrialsAcq, 1)); % CS type column
    tempMat(1:nrUSAcq,2) = 1; % US [yes/no] column
    tempMat(:,3) = createExpIntDurs(trialsPerBlockAcq, minMaxItiSec, meanItiSec*trialsPerBlockAcq); % ITI column
    tempMat = tempMat(randperm(trialsPerBlockAcq),:); % shuffle trial order
    
    % write values into trialMat and ratingAfterTrial
    trialCount = trialCount + trialsPerBlockAcq;
    ratingCount = ratingCount + 1;
    trialMat((trialCount-trialsPerBlockAcq+1):trialCount,:) = tempMat;
    ratingAfterTrial(ratingCount) = trialCount;
end

% for extinction trials
for bl = 1:nrBlocksExt
    tempMat = zeros(trialsPerBlockExt,3);
    tempMat(:,1) = sort(repmat(conds', nrTrialsExt, 1)); % CS type column
    tempMat(:,3) = createExpIntDurs(trialsPerBlockExt, minMaxItiSec, meanItiSec*trialsPerBlockExt); % ITI column
    tempMat = tempMat(randperm(trialsPerBlockExt),:); % shuffle trial order
    
    % write values into trialMat and ratingAfterTrial
    trialCount = trialCount + trialsPerBlockExt;
    ratingCount = ratingCount + 1;
    trialMat((trialCount-trialsPerBlockExt+1):trialCount,:) = tempMat;
    ratingAfterTrial(ratingCount) = trialCount;
end

% flicker sequences (on/off frames)   
flickerVec = genFlickVec(framPerCyc, flickDurFrames, flickMod);
flickerVec = squeeze(flickerVec);


%% PTB code is in try loop
try
    %% get all info for screen & window
    screens = Screen('Screens');
    screenNumber = max(screens);
    w = Screen('OpenWindow', screenNumber, backgroundCol);
        
    % set alpha blending mode for manipulating picture transparency
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % Compute screen's PPI via window size (pix) and display size (mm)
    ppi = Screen('WindowSize', w) / Screen('DisplaySize', screenNumber) * 25.4;
    % translate visual angle into pixels for picture scaling
    imSizePix = pixFromAngle(imSizeAng, seatDistInch, ppi); % size of Gabors
    spatFreqPix = cycPerDeg/pixFromAngle(1, seatDistInch, ppi); % spatial frequency of Gabors

    % set font size
    Screen('TextSize', w, fontSize);    

    %% prepare Gabor patches

    textureVec = NaN(nrCS,1);
    for textInd = 1:nrCS
        textureVec(textInd) = createPseudoGaborTexture(w,imSizePix,rotStim(textInd),spatFreqPix,[gratingDark backgroundCol gratingBright]);
    end % textInd loop
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   HERE COMES THE ACTUAL EXPERIMENT   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HideCursor;
    
    % wait for the first MRI pulse and log system time
    timeFirstTR = waitForScanTriggerKb(TRtriggerCodes);

    % welcome screen, initiated by experimenter terminated by participant
    DrawFormattedText(w, welcomeMsg_expStart, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);
        
    %% trials start
    

    for trial = 1 : size(trialMat,1)
        % present trial while logging actual flicker duration
        presFix(w, trialMat(trial,3));
        
        %%% WORK AROUND UNTIL WE CAN INITIALIZE ACTUAL PORT OBJECTS
        serPort = 'dummy';
        lptPortObj = 'dummy';
        %%%
        
        % present CS/GS, log actual stimulus duration + sys time at onset
        if trialMat(trial, 2) == 1 % if there is a US to presented
            [actFlickDur, onsetTime] = presFlickShockMRI(w, ...
                flickerVec, [], textureVec(trialMat(trial,1)), ...
                shockOnsetFrame, shockDurFrame, ...
                lptPortObj, adressesLPT, lptCodes, serPort, 'fufufufu99fufufu');
        else % no US
            [actFlickDur, onsetTime] = presFlickShockMRI(w, ...
            flickerVec, [], textureVec(trialMat(trial,1)), ...
            0, 0, ...
            lptPortObj,  adressesLPT, lptCodes, serPort, 'fufufufu99fufufu');
        end % US vs NoUS conditional

        % write trial parameters to logfile
        % part ID, CS+ permutation, trial #, trial type, stimulus onset (sys time),
        % stimulus onset (since first TR), simulus duration, iti duration;
        % condCode will be 10 for CS+ w/o US, 11, for CS+ w/ US, 20/30/40/...
        % for GS1/GS2/GS3/...
        timeSinceFirstTR = onsetTime - timeFirstTR;
        condCode = 10*trialMat(trial,1) + trialMat(trial,2); 
        trialOutVec = [subNo, csPerm trial, condCode, onsetTime, ...
                       timeSinceFirstTR, actFlickDur, trialMat(trial,3)];
        dlmwrite(logFileName, trialOutVec, '-append');

        % check if this is the last trial of a block (==> ratings)
        if ismember(trial, ratingAfterTrial)
            % fix cross for the BOLD response to play out before ratings
            presFix(w, 3*TRdur);

            % wait for TR pulse to compensate for timing errors across participants
            waitForScanTriggerKb(TRtriggerCodes);

            % run rating routine and write output into file (append mode)
            currRatings = rateCSMRI(w, buttonsLeft, buttonsRight, buttonsOK, ...
                                    ratingDur, pauseBtwRat, textureVec);
            dlmwrite(ratFileName, currRatings, '-append');

            % another wait-for-TR to compensate for timing errors
            waitForScanTriggerKb(TRtriggerCodes);

        end % rating conditional
        

    end % trial loop
        
    % final screen, disappears automatically
    DrawFormattedText(w, goodbyeMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(intructDur);

    % tidy up
    Screen('CloseAll');
    ShowCursor;
    fclose('all');    
    IOPort('CloseAll');
    
    %%% that's it (if everything worked) %%%

    % in case of errors:
catch   
    % ready to end the exp   
    Screen('CloseAll');
    ShowCursor;
    fclose('all');

    psychrethrow(psychlasterror);
end

end % paradigm function

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   HERE COME THE CUSTOM FUCTIONS   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generates vectors of luminance/opaqueness values for flicker
% presentation (one value per frame); takes frames per cylce, total number
% of frames and the shape of the modulation ('box' for on/off, 'sin' for
% sinusoidal)

function [flickVec] = genFlickVec(framesPerCyc, framesTotal, shapeFunc)
    if strcmp(shapeFunc, 'box')
        flickVec = repmat([ones(framesPerCyc/2,1); zeros(framesPerCyc/2,1)], ...
                          framesTotal/framesPerCyc, 1)';
    elseif strcmp(shapeFunc, 'sin')
        flickVec = repmat(sin(1.5*pi:(2*pi/framesPerCyc):(3.5*pi-2*pi/framesPerCyc)) / 2 + .5, ...
                          1, framesTotal/framesPerCyc);
    end
end


% computes # of pixels required for picture based on the required visual
% angle, seating distance in inches, and the PPI of the screen; can take
% multiple and compute vector of values in one call
function nrPix = pixFromAngle(visAngle, distInch, ppi)
    nrPix = round(tan(visAngle/360*pi) * ppi * 2 * distInch, 0);
end


% creates a grating as texture object for Psychtoolbox on basis of a Gabor
% patch. Instead of smooth light-dark transitions, this grating's "bars" 
% have sharp edges. The grating has a circular shape.
% sizePix = [x, y] size in pixels
% angle = orientation of Grating in degrees with 0 == vertical orientation
% and clockwise rotation
% spatFreqPix = spatialFrequency of grating in units of cycles per pixel
% brightness = set of 3 values to determine brightness of stimulus components:
% [dark bars, background, bright bars]; e.g. [0 127 255] for black & white
% bars in front of gray background
function pseudoGaborTexture = createPseudoGaborTexture(window, sizePix, angle, spatFreqPix, brightness)

% to be safe, round size parameter
sizePix = round(sizePix,0);

% compute values for Gabor pattern; formula provided by Maeve Boylan (thanks!)
[x,y] = meshgrid(-sizePix(1)/2 : sizePix(1)/2, -sizePix(2)/2 : sizePix(2)/2); 
gaborFunction = (exp(-((x/sizePix(1)*4).^2)-((y/sizePix(2)*4).^2)) .* sin(cos(angle*pi/180)*(2*pi*spatFreqPix)*x + sin(angle*pi/180)*(2*pi*spatFreqPix)*y)); 

% "categorize" Gabor values into the dark, background and bright
gaborPixels = NaN(size(gaborFunction));
gaborPixels(gaborFunction < -.02) = brightness(1);
gaborPixels(gaborFunction >= -.02 & gaborFunction <= .02) = brightness(2);
gaborPixels(gaborFunction > .02) = brightness(3);

% create texture object
pseudoGaborTexture = Screen('MakeTexture', window, gaborPixels, 0, 4);
end


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


% present fixation cross; cross size is hardcoded and relative to window
% size (pixels), duration is given in seconds
function presFix(window, fixDurSec)
    [winSize(1), winSize(2)] = Screen('WindowSize', window);
    winCenter = winSize ./ 2;
    
    % hard-coded length of fixation cross arme relative to screen height
    fixL = round(winSize(2)*.02, 0);
    
    % Present fixation cross
    Screen('DrawLines', window, [-fixL fixL 0 0; 0 0 -fixL fixL], 7, 0, winCenter);
    Screen('Flip', window);
    WaitSecs(fixDurSec);
end


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


% PTB rating routine for a threat conditioning paradigm with habituation, 
% acquistion, and extinction; items have partricipants rate stimulus
% valence, arousal, their fear, and US expectancy; Likert scales with
% slider response mode;
% The routine will first acquire valence ratings for all CS, then arousal,
% fear, and expectancy; order of CS is randomized at beginning of rating
% routine (but constant within routine);
% item texts, scale ticks & labels can be easily adapted (see "header");
% scales everything relative to screen (except for font size)
% Ratings have fixed duration ==> rating = selection at time-out
% buttonsLeft, buttonsRight, buttonsOK are given in DEC; currently the OK
% buttons just make the marker change color
% maxDurInSec: fixed duration of each rating in seconds
% pauseDur: fixed duration of fixation cross after each rating in seconds
% csTextureVec: vector with textures for all CS to be rated
% output: ratingMat = 1 x nrRatings vector with all val ratings (stim1, stim2, ...),
% then all arousal, fear, expectancy ratings

function ratingMat = rateCSMRI(window, buttonsLeft, buttonsRight, buttonsOK, maxDurInSec, pauseDur, csTextureVec)
    % get system time at start of routine
    startTime = GetSecs;
    
    %% header with parameters
    textValence = 'How unpleasant was this pattern to you?';
    anchorsValence = {'very pleasant', 'very unpleasant'};
    textArousal = 'How arousing was this pattern to you?';
    anchorsArousal = {'not arousing at all', 'very arousing'};
    textFear = 'How fearful did you feel when you saw this pattern?';
    anchorsFear = {'not fearful at all', 'very fearful'};
    textExpectancy = 'How likely was it that a shocked followed this patern?';
    anchorsExpectancy = {'never followed', 'always followed'};
    
    labelsUnipolar = 0:10;
    labelsBipolar = -5:5;
    labelsExpectancy = 0:10:100;

    textCol = [0 0 0]; % color for text & scale
    noSelCol = [255 0 0]; % color for cursor before selection is confirmed
    selMadeCol = [0 255 0]; % color for cursor after selection is confirmed

    %% matrices that we need
    orderVec = randperm(length(csTextureVec));
    valenceMat = NaN(1, length(csTextureVec));
    arousalMat = NaN(1, length(csTextureVec));
    fearMat = NaN(1, length(csTextureVec));
    expectancyMat = NaN(1, length(csTextureVec));

    % automatic computation of sizes relative to screen
    % global values / coordinates
    [xTotal, yTotal] = Screen('WindowSize', window);
    xCenter = xTotal / 2; %yCenter = yTotal ./ 2;

    % x coordinates
    widthLeft = round(xTotal * .2, 0);
    widthRight = round(xTotal * .8, 0);
    widthTicksUnipolar = round(widthLeft : (widthRight-widthLeft)/(length(labelsUnipolar)-1) : widthRight, 0);
    widthTicksBipolar = round(widthLeft : (widthRight-widthLeft)/(length(labelsBipolar)-1) : widthRight, 0);
    widthTicksExpectancy = round(widthLeft : (widthRight-widthLeft)/(length(labelsExpectancy)-1) : widthRight, 0);

    % y coordinates
    heightStim = round(yTotal * .3, 0);
    heightQuestion = round(yTotal * .6, 0);
    heightScale = round(yTotal * .8, 0);
    heightLabels = round(yTotal * .9, 0);
    heightAnchors = round(yTotal * .95, 0);

    % stimulus size
    sizeStim = [yTotal*.3, yTotal*.3];
    destRectCs = [xCenter-sizeStim(1)/2, heightStim-sizeStim(2)/2, ...
                  xCenter+sizeStim(1)/2, heightStim+sizeStim(2)/2];
    
    % x and y coordinates for scale ticks to use with Screen('DrawLines')
    coordsTicksUnipolar = [sort(repmat(widthTicksUnipolar, 1, 2)); ...
                           repmat([heightScale-.05*yTotal,  heightScale+.05*yTotal], ...
                           1, length(labelsUnipolar))];
    coordsTicksBipolar = [sort(repmat(widthTicksBipolar, 1, 2)); ...
                           repmat([heightScale-.05*yTotal,  heightScale+.05*yTotal], ...
                           1, length(labelsBipolar))];
    coordsTicksExpectancy = [sort(repmat(widthTicksExpectancy, 1, 2)); ...
                             repmat([heightScale-.05*yTotal,  heightScale+.05*yTotal], ...
                             1, length(labelsExpectancy))];

    %% Valence Rating
    % selection index variables
    selMin = 1;
    selMax = length(labelsBipolar);
    selCent = ceil((selMax+selMin)/2);
    currTime = 0; % initialize currTime with some value that is smaller than any GetSecs

    for StimI = 1:length(csTextureVec)
        
        % for the very first rating we use the startTime of the function to
        % account for time taken up by initializing all variables
        if StimI > 1
            startTime = GetSecs;
        end
        
        % initiate cursor properties
        currSelect = selCent;
        pointCol = noSelCol;
        
        while (currTime-startTime) < maxDurInSec
            
            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textValence, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksBipolar, [], textCol);
            Screen('DrawLine', window, pointCol, ...
                               widthTicksBipolar(currSelect), heightScale-.05*yTotal, ...
                               widthTicksBipolar(currSelect), heightScale+.05*yTotal, 3);
            
            % Draw Labels
            for lab = 1:length(labelsBipolar)
                DrawFormattedText(window, int2str(labelsBipolar(lab)), ...
                                  'center', heightLabels, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicksBipolar(lab), heightLabels, widthTicksBipolar(lab), heightLabels]);
            end
            
            % Draw Anchors
            DrawFormattedText(window, anchorsValence{1}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksBipolar(1), heightAnchors, widthTicksBipolar(1), heightAnchors]);
            DrawFormattedText(window, anchorsValence{2}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksBipolar(end), heightAnchors, widthTicksBipolar(end), heightAnchors]);
            
            % Flip it!
            Screen('Flip', window)

            % wait for button press, change selected rating value & cursor color
            [~, buttonsPressed, ~] = KbStrokeWait([], startTime + maxDurInSec);
            if sum(buttonsPressed(buttonsLeft)) > 0 && currSelect > selMin
                currSelect = currSelect - 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsRight)) > 0 && currSelect < selMax
                currSelect = currSelect + 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsOK)) > 0
                pointCol = selMadeCol;
            end
            
            % log system time
            currTime = GetSecs;
        end
        
        % log rating
        valenceMat(orderVec(StimI)) = labelsBipolar(currSelect);
        
        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end


    %% Arousal Rating
    % selection index variables
    selMin = 1;
    selMax = length(labelsUnipolar);
    selCent = ceil((selMax+selMin)/2);
    
    for StimI = 1:length(csTextureVec)
        % log system time at rating onset and initiate cursor properties
        startTime = GetSecs;
        currSelect = selCent;
        pointCol = noSelCol;
        
        while (currTime-startTime) < maxDurInSec
            
            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textArousal, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksUnipolar, [], textCol); 
            Screen('DrawLine', window, pointCol, ...
                               widthTicksUnipolar(currSelect), heightScale-.05*yTotal, ...
                               widthTicksUnipolar(currSelect), heightScale+.05*yTotal, 3);
            
            % Draw Labels
            for lab = 1:length(labelsUnipolar)
                DrawFormattedText(window, int2str(labelsUnipolar(lab)), ...
                                  'center', heightLabels, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicksUnipolar(lab), heightLabels, widthTicksUnipolar(lab), heightLabels]);
            end
            
            % Draw Anchors
            DrawFormattedText(window, anchorsArousal{1}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksUnipolar(1), heightAnchors, widthTicksUnipolar(1), heightAnchors]);
            DrawFormattedText(window, anchorsArousal{2}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksUnipolar(end), heightAnchors, widthTicksUnipolar(end), heightAnchors]);
            
            % Flip it!
            Screen('Flip', window)

            % wait for button press, change selected rating value & cursor color
            [~, buttonsPressed, ~] = KbStrokeWait([], startTime + maxDurInSec);
            if sum(buttonsPressed(buttonsLeft)) > 0 && currSelect > selMin
                currSelect = currSelect - 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsRight)) > 0 && currSelect < selMax
                currSelect = currSelect + 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsOK)) > 0
                pointCol = selMadeCol;
            end

            % log system time
            currTime = GetSecs;
        end
        
        % log rating
        arousalMat(orderVec(StimI)) = labelsUnipolar(currSelect);
        
        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end


    %% Fear Rating
    % selection index variables
    selMin = 1;
    selMax = length(labelsUnipolar);
    selCent = ceil((selMax+selMin)/2);
    
    for StimI = 1:length(csTextureVec)
        % log system time at rating onset and initiate cursor properties
        startTime = GetSecs;
        currSelect = selCent;
        pointCol = noSelCol;
        
        while (currTime-startTime) < maxDurInSec
            
            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textFear, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksUnipolar, [], textCol); 
            Screen('DrawLine', window, pointCol, ...
                               widthTicksUnipolar(currSelect), heightScale-.05*yTotal, ...
                               widthTicksUnipolar(currSelect), heightScale+.05*yTotal, 3);
            
            % Draw Labels
            for lab = 1:length(labelsUnipolar)
                DrawFormattedText(window, int2str(labelsUnipolar(lab)), ...
                                  'center', heightLabels, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicksUnipolar(lab), heightLabels, widthTicksUnipolar(lab), heightLabels]);
            end
            
            % Draw Anchors
            DrawFormattedText(window, anchorsFear{1}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksUnipolar(1), heightAnchors, widthTicksUnipolar(1), heightAnchors]);
            DrawFormattedText(window, anchorsFear{2}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksUnipolar(end), heightAnchors, widthTicksUnipolar(end), heightAnchors]);
            
            % Flip it!
            Screen('Flip', window)

            % wait for button press, change selected rating value & cursor color            
            [~, buttonsPressed, ~] = KbStrokeWait([], startTime + maxDurInSec);
            if sum(buttonsPressed(buttonsLeft)) > 0 && currSelect > selMin
                currSelect = currSelect - 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsRight)) > 0 && currSelect < selMax
                currSelect = currSelect + 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsOK)) > 0
                pointCol = selMadeCol;
            end

            % log system time
            currTime = GetSecs;
        end
        
        % log rating
        fearMat(orderVec(StimI)) = labelsUnipolar(currSelect);
        
        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end

    %% Expectancy Rating
    % selection index variables
    selMin = 1;
    selMax = length(labelsUnipolar);
    selCent = ceil((selMax+selMin)/2);
    
    for StimI = 1:length(csTextureVec)
        % log system time at rating onset and initiate cursor properties
        startTime = GetSecs;
        currSelect = selCent;
        pointCol = noSelCol;
        
        while (currTime-startTime) < maxDurInSec
            
            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textExpectancy, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksExpectancy, [], textCol); 
            Screen('DrawLine', window, pointCol, ...
                               widthTicksExpectancy(currSelect), heightScale-.05*yTotal, ...
                               widthTicksExpectancy(currSelect), heightScale+.05*yTotal, 3);
            
            % Draw Labels
            for lab = 1:length(labelsExpectancy)
                DrawFormattedText(window, [int2str(labelsExpectancy(lab)) '%'], ...
                                  'center', heightLabels, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicksExpectancy(lab), heightLabels, widthTicksExpectancy(lab), heightLabels]);
            end
            
            % Draw Anchors
            DrawFormattedText(window, anchorsExpectancy{1}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksExpectancy(1), heightAnchors, widthTicksExpectancy(1), heightAnchors]);
            DrawFormattedText(window, anchorsExpectancy{2}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksExpectancy(end), heightAnchors, widthTicksExpectancy(end), heightAnchors]);
            
            % Flip it!
            Screen('Flip', window)

            % wait for button press, change selected rating value & cursor color
            [~, buttonsPressed, ~] = KbStrokeWait([], startTime + maxDurInSec);
            if sum(buttonsPressed(buttonsLeft)) > 0 && currSelect > selMin
                currSelect = currSelect - 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsRight)) > 0 && currSelect < selMax
                currSelect = currSelect + 1;
                pointCol = noSelCol;
            elseif sum(buttonsPressed(buttonsOK)) > 0
                pointCol = selMadeCol;
            end

            % log system time
            currTime = GetSecs;
        end
        
        % log rating
        expectancyMat(orderVec(StimI)) = labelsExpectancy(currSelect);
        
        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end    

    % concatenate ratings for the different items; result is a 1 x nrRatings vector
    ratingMat = cat(2, valenceMat, arousalMat, fearMat, expectancyMat);
end