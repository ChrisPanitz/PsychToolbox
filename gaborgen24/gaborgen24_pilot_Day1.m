function gaborgen24_pilot_Day1(subNo, csPerm)
%% SET TO ZERO FOR ACTUAL EXPERIMENT
Screen('Preference', 'SkipSyncTests', 0);


%% Header
% trial parameters 
nrBlocksHab = 1; % # of blocks in acquisition; ratings after each block
nrTrialsHab = 10; % # of trials per stimulus and per block

nrBlocksAcq = 2; % # of blocks in acquisition; ratings after each block
nrTrialsAcq = 16; % # of trials per stimulus and per block
percPairing = .5;
nrBoostersMassed = 3;
nrBoostersRandom = 3;
maxGSbtwBooster = 2;

nrBlocksExt = 1; % # of blocks in acquisition; ratings after each block
nrTrialsExt = 15; % # of trials per stimulus and per block

minMaxItiSec = [5.5 15.4]; % min & max ITI in sec, taken from exponential distribution
meanItiSec = 7; % mean ITI duration in sec

% How many sequential presentations of the same CS are allowed
% if no restrictions wanted ==> set to ridiculously high number
maxSeqCS = 2; 

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
flickDurFrames = 360; % duration of whole stimulus in frames; 360/120 = 3 s

% US parameters
soundSamp = 44100; % sampling rate of sound presentation
soundDur = 1; % duration of noise stimulus in s
soundIntensity = .95; % choose a value between 0 (muted) and 1 (max intensity)
%soundOnset = 241; % sound Onset in frames; 241 = after 2 seconds
soundOnset = 2; % sound Onset in seconds

% rating & instruction parameters
pauseBtwRat = .5; % pause duration between two ratings, in sec
buttonsOK = 1; % indeces of mouse buttons to select rating

% aethetics
backgroundCol = 127; % bakcground color; 127 = mid gray
gratingDark = 0; % dark grating stripes = 0 = black
gratingBright = 255; % bright grating stripes = 0 = white
fontSize = 32; % well... font size

% for logfiles and rating files
logFileFolder = '/home/andreaskeil/Desktop/As_Exps/gaborgen24/logfiles/';
logFilePrefix = 'gaborgen24_pilot_Day1'; % TODO: come up with something more specific
ratFileFolder = '/home/andreaskeil/Desktop/As_Exps/gaborgen24/ratings/';
ratFilePrefix = 'gaborgen24_pilot_Day1'; % TODO: come up with something more specific
ratAllFileFolder = '/home/andreaskeil/Desktop/As_Exps/gaborgen24/ratings/';
ratAllFilePrefix = 'gaborgen24_pilot_Day1'; % TODO: come up with something more specific

% strings of instructions to participants          
welcomeMsg_noiseTest = ['Welcome and thank you for participating in our\n' ...
                        'experiment. In a moment, you will hear a loud sound\n' ...
                        'that will be presented to you multiple times during\n' ...
                        'this experiment.\n\n' ...
                        'Please click any mouse button whenever you are ready.'];
welcomeMsg_expStart = ['Thank you!\n\n' ...
                       'The experimenter will start the experiment soon.'];
welcomeMsg_partStart = ['Everything has been set up now. In the upcoming\n' ...
                        'task you will be presented flickering pictures\n' ...
                        'in the center of the screen. Please pay attention\n' ...
                        'to them. The loud sound also may occur from time\n' ...
                        'to time. Every now and then you will be asked to\n' ...
                        'answer questions about how you experience the pictures.\n'...
                        'After each round of questions, you can take a \n' ...
                        'short break before continuing. The whole task will\n' ...
                        'last about 60 minutes.\n\n\n' ...
                        'If you have questions, please let the\n' ...
                        'experimenter know now. When you are ready\n' ...
                        'you can start the task by pressing any mouse button.'];
breakMsg = ['You can take a short break now.\n\n' ...
            'When you are ready, please continue the experiment\n' ...
            'by pressing any mouse button.'];
goodbyeMsg = ['You have completed the task. Thank you for your participation!\n\n' ...
              'The experimenter will be with you in a moment.'];

          

%% clearing and initializing stuff
IOPort('CloseAll');
clc;

% Setting rng seed; "old" rand function (instead of rng) was used for
% reasons of compatibility with lab computer
rand('state',sum(100*clock));

AssertOpenGL; 

InitializePsychSound(1)

%%% TODO: INCLUDE PORT STATEMENTS WHEN PORT AVAILABLE %%%
% Serial port
[serPort, ~] = IOPort('OpenSerialPort', '/dev/ttyS0', 'Baudrate = 9600');
   
%% Prepare logging

% files to write trial & rating data + variable names for header line
% logfiles
logFileName = strcat(logFileFolder, logFilePrefix, '_', num2str(subNo), '_logfile.dat');
fIDLog = fopen(logFileName, 'w');
fprintf(fIDLog,'partNo,csPerm,trial,condCode,actFlickDur,itiDur\n');
fclose(fIDLog);
% ratings
ratingsHeader = ['partInd,ratInd,' ...
                 'val_csp,val_g1,val_g2,val_g3,ar_csp,ar_g1,ar_g,ar_g3,' ...
                 'fear_csp,fear_g1,fear_g2,fear_g3,exp_csp,exp_g1,exp_g2,exp_g3,' ...
                 'valX_csp,valX_g1,valX_g2,valX_g3,arX_csp,arX_g1,arX_g,arX_g3,' ...
                 'fearX_csp,fearX_g1,fearX_g2,fearX_g3,expX_csp,expX_g1,expX_g2,expX_g3,' ...
                 'valDur_csp,valDur_g1,valDur_g2,valDur_g3,arDur_csp,arDur_g1,arDur_g,arDur_g3,' ... 
                 'fearDur_csp,fearDur_g1,fearDur_g2,fearDur_g3,expDur_csp,expDur_g1,expDur_g2,expDur_g3' ...
                 '\n'];
% individual rating files
ratFileName = strcat(ratFileFolder, ratFilePrefix, '_', num2str(subNo), '_ratings.dat');
fIDRat = fopen(ratFileName, 'w');
fprintf(fIDRat, ratingsHeader);
fclose(fIDRat);
% rating file for whole sample
ratAllFileName = strcat(ratAllFileFolder, ratAllFilePrefix, '_ratingsAll.dat');
if ~exist(ratAllFileName, 'file')
    fIDRatAll = fopen(ratAllFileName, 'w');
    fprintf(fIDRatAll, ratingsHeader);
    fclose(fIDRatAll);
end

currentRatingInd = 1; % for logging which round of ratings we're in

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
trialMat = NaN(nrBlocksHab*trialsPerBlockHab + ...partStart
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
    tempMat = NaN(trialsPerBlockHab,3);
    tempMat(:,1:2) = pseudoShuffleWithBooster(conds, nrTrialsHab, 0, ...
                     0, [], maxSeqCS, 0, true);
    % ITI columns
    tempMat(:,3) = createExpIntDurs(trialsPerBlockHab, minMaxItiSec, meanItiSec*trialsPerBlockHab); 
    
    % write values into trialMat and ratingAfterTrial
    trialCount = trialCount + trialsPerBlockHab;
    ratingCount = ratingCount + 1;
    trialMat((trialCount-trialsPerBlockHab+1):trialCount,:) = tempMat;
    ratingAfterTrial(ratingCount) = trialCount;
end

% for acquisition trials
for bl = 1:nrBlocksAcq
   tempMat = NaN(length(conds)*nrTrialsAcq,3);
   % CS & US column - boosters only in first block
   if bl == 1
       tempMat(:,1:2) = pseudoShuffleWithBooster(conds, nrTrialsAcq, nrBoostersMassed, ...
           nrBoostersRandom, maxGSbtwBooster, maxSeqCS, percPairing, true);
   elseif bl > 1
       tempMat(:,1:2) = pseudoShuffleWithBooster(conds, nrTrialsAcq, 0, ...
           0, maxGSbtwBooster, maxSeqCS, percPairing, true);
   end
   % ITI column
   tempMat(:,3) = createExpIntDurs(trialsPerBlockAcq, minMaxItiSec, meanItiSec*trialsPerBlockAcq);

   % write values into trialMat and ratingAfterTrial
   trialCount = trialCount + trialsPerBlockAcq;
   ratingCount = ratingCount + 1;
   trialMat((trialCount-trialsPerBlockAcq+1):trialCount,:) = tempMat;
   ratingAfterTrial(ratingCount) = trialCount;
end

% for extinction trials
for bl = 1:nrBlocksExt
    tempMat = NaN(trialsPerBlockExt,3);
    tempMat(:,1:2) = pseudoShuffleWithBooster(conds, nrTrialsExt, 0, ...
                     0, [], maxSeqCS, 0, true);
    % ITI columns
    tempMat(:,3) = createExpIntDurs(trialsPerBlockExt, minMaxItiSec, meanItiSec*trialsPerBlockExt);

    % write values into trialMat and ratingAfterTrial
    trialCount = trialCount + trialsPerBlockExt;
    ratingCount = ratingCount + 1;
    trialMat((trialCount-trialsPerBlockExt+1):trialCount,:) = tempMat;
    ratingAfterTrial(ratingCount) = trialCount;
end

% flicker sequences (on/off frames)   
flickerVec = genFlickVec(framPerCyc, flickDurFrames, flickMod);
flickerVec = squeeze(flickerVec)';

% white noise US
% if sound onset is given in frames
%whiteNoise = rand(1, soundSamp*soundDur) * soundIntensity; 
% if soundOnset is given in seconds
whiteNoise = [zeros(1,soundSamp*soundOnset), rand(1, soundSamp*soundDur)] * soundIntensity;


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
    
    HideCursor();
        
    % playing the US noise once upfront
    DrawFormattedText(w, welcomeMsg_noiseTest, 'center', 'center');
    Screen('Flip', w);
    waitForClick();
    DrawFormattedText(w, 'sound incoming...', 'center', 'center');
    Screen('Flip', w);
    Snd('Play', whiteNoise, soundSamp);
    WaitSecs(soundOnset + soundDur + 1); % works only correctly when sound latencies given in seconds
    
    % welcome screen, terminated by experimenter key press
    DrawFormattedText(w, welcomeMsg_expStart, 'center', 'center');
    Screen('Flip', w);
    KbStrokeWait();
    
    % instruction screen for task, terminated by participant mouse click
    DrawFormattedText(w, welcomeMsg_partStart, 'center', 'center');
    Screen('Flip', w);
    waitForClick();
    presFix(w, 5);
        
    
    %% trials start

    for trial = 1 : size(trialMat,1)
        % present CS/GS, log actual stimulus duration
        % if sound onset given in frames
%         if trialMat(trial, 2) == 1 % if there is a US to presented
%             actFlickDur = presFlickSound(w, ...
%                 flickerVec, [], textureVec(trialMat(trial,1)), ...
%                 whiteNoise, soundOnset, soundSamp, serPort);
        %if sound onset given in seconds
        if trialMat(trial, 2) == 1 % if there is a US to presented
            actFlickDur = presFlickSound(w, ...
                flickerVec, [], textureVec(trialMat(trial,1)), ...
                whiteNoise, 1, soundSamp, serPort);
        else % no US
            actFlickDur = presFlickSound(w, ...
                flickerVec, [], textureVec(trialMat(trial,1)), ...
                [], soundOnset, soundSamp, serPort);
        end % US vs NoUS conditional

        % present trial while logging actual flicker duration
        presFix(w, trialMat(trial,3));
        
        % write trial parameters to logfile
        % part ID, CS+ permutation, trial #, trial type, stimulus onset (sys time),
        % stimulus onset (since first TR), simulus duration, iti duration;
        % condCode will be 10 for CS+ w/o US, 11, for CS+ w/ US, 20/30/40/...
        % for GS1/GS2/GS3/...
        condCode = 10*trialMat(trial,1) + trialMat(trial,2); 
        trialOutVec = [subNo, csPerm, trial, condCode, actFlickDur, trialMat(trial,3)];
        dlmwrite(logFileName, trialOutVec, '-append');

        % check if this is the last trial of a block (==> ratings)
        if ismember(trial, ratingAfterTrial)
            % fix cross before ratings
            presFix(w, 5);
            ShowCursor();
            
            % run rating routine and write output into file (append mode)
            currRatings = rateCSmouse(w,buttonsOK,pauseBtwRat,false,textureVec);
            dlmwrite(ratFileName, [subNo, currentRatingInd, currRatings], '-append');
            dlmwrite(ratAllFileName, [subNo, currentRatingInd, currRatings], '-append');
            currentRatingInd = currentRatingInd+1;
            HideCursor();

            % break (but not after last trial)
            if trial < size(trialMat, 1)
                presFix(w, 1);
                DrawFormattedText(w, breakMsg, 'center', 'center');
                Screen('Flip', w);
                waitForClick();
            end
            
            % fix cross before continuing
            presFix(w, 5);
            
        end % rating conditional
        
    end % trial loop
        
    % final screen, disappears automatically
    DrawFormattedText(w, goodbyeMsg, 'center', 'center');
    Screen('Flip', w);
    KbStrokeWait();

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

function shuffledMat = pseudoShuffleWithBooster(condVec, trPerCS, nrBoostersMassed, nrBoostersRandom, maxTrBtwBooster, maxSeqCS, percPairing, pairingProbLocal) 
    
    clusteredCS = 1; % any value > 0
    while clusteredCS > 0
        % trials for the booster "phase"
        boostMassInd = NaN(nrBoostersMassed,1);
        if sum(nrBoostersMassed + nrBoostersRandom) > 0
            boostMassInd(1) = randi(maxTrBtwBooster+1);
        end
        for i = 2:length(boostMassInd)
            boostMassInd(i) = boostMassInd(i-1) + randi(maxTrBtwBooster+1);
        end
        boostMat = zeros(length(boostMassInd),2);
        boostMat(boostMassInd,1) = 1;
        boostMat(boostMat(:,1) == 0,1) = randi([2,max(condVec)], sum(boostMat(:,1) == 0), 1);

        remainTrials = repmat(trPerCS, 1, length(condVec)) - histcounts(boostMat(:,1), 1:condVec(end)+1);

        afterBoostMat = zeros(sum(remainTrials),2);
        for condI = condVec
            afterBoostMat(sum(remainTrials(1:condI-1))+1 : sum(remainTrials(1:condI)), 1) = condI;
        end
    
        % shuffle after-massed-booster trial order and concatenate all trials
        afterBoostMat = afterBoostMat(randperm(sum(remainTrials)),:); % shuffle trial order
        shuffledMat = cat(1, boostMat, afterBoostMat);

         % 0 if the same CS occurs twice in sequence
        diffTrials = diff(shuffledMat(:,1));
        % adds neighboring values in groups of [maxSeqCS]
        % => [maxSeqCS] times no change is too much
        % "abs" is to avoid sums of zero from other values than only zeros
        seqVec = filter(ones(maxSeqCS,1), 1, abs(diffTrials)); 
        % first [maxSeqCS-1] values are not based on the full set of [maxSeqCS] elements
        seqVec = seqVec(maxSeqCS:end); 
        clusteredCS = sum(seqVec == 0); % if there are zeros ==> keep shuffling
    end

    % US [yes/no] column
    if pairingProbLocal == true
        nrUSafterBooster = floor((trPerCS-nrBoostersMassed-nrBoostersRandom) * percPairing);
    elseif pairingProbLocal == false
        nrUSafterBooster = floor(trPerCS*percPairing - nrBoostersMassed - nrBoostersRandom);
    end
    usVecAftBoost = zeros(remainTrials(1)-nrBoostersRandom,1);
    usVecAftBoost(1:nrUSafterBooster) = 1;
    usVecAftBoost = usVecAftBoost(randperm(length(usVecAftBoost)));
    shuffledMat(shuffledMat(:,1) == 1, 2) = cat(1, ones(nrBoostersMassed+nrBoostersRandom,1), usVecAftBoost);
end



% wait for any mouse click, detected via the GetMouse function in
% Psychtoolbox
function [x, y, buttons] = waitForClick
    buttons = 0;
    while ~any(buttons) % wait for press
        [x, y, buttons] = GetMouse;
    
        % Wait 10 ms before checking the mouse again to prevent
        % overload of the machine at elevated Priority()
        WaitSecs(0.01);
    end
end



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


% presents a flickering image texture; flickVec is a one-dimensional vector
% with opaqueness values from 0 to 1 for each frame; imageSize can be empty
% (no scaling of original image) or have two values for size in pixels ([x
% y]); texture is texture object; returns actual flicker stimulus duration
% in seconds obtained with tic & toc commands
function actFlickDur = presFlickSound(window, flickVec, imageSize, texture, soundStim, soundOn, soundSamp, portOut)
    if isempty(imageSize)
        imageCoordinates = [];
    else
        [winSize(1), winSize(2)] = Screen('WindowSize', window);
        winCenter = winSize ./ 2;
        imageCoordinates = [winCenter(1)-imageSize(1)/2 winCenter(2)-imageSize(2)/2 winCenter(1)+imageSize(1)/2 winCenter(2)+imageSize(2)/2];
    end
    
    % send marker via IO Port and log time for measuring flicker duration
    %[nwritten, when, errmsg, prewritetime, postwritetime, lastchecktime] = IOPort('Write', s3, 'fufufufu99fufufu');
    IOPort('Write', portOut, 'fufufufu99fufufu');
    startTime = GetSecs();
    
    % Present Flicker Stimulus and US
    for frame = 1:soundOn-1
        Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], flickVec(frame));
        Screen('Flip', window);
    end
    
    if ~isempty(soundStim)
        Snd('Play', soundStim, soundSamp);
    end
    
    for frame = soundOn:size(flickVec, 1)
        Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], flickVec(frame));
        Screen('Flip', window);
    end
    
    % check time passed since before flicker stim
    actFlickDur = GetSecs() - startTime;
end




% returns random interval duration in seconds, following an exponential
% function with mu as mean of the base function (exprnd function in MATLAB) 
% and a vector with minimal and maximum duration of interval; the value of
% 5 for computing lambda was obtained by trying out
function intInSec = randExpoInt(minMaxSec)
    intInSec = NaN;
    while ~(intInSec >= minMaxSec(1) && intInSec <= minMaxSec(2))
        lambda = 5/diff(minMaxSec);
        intInSec = minMaxSec(1) - log(rand(1,1))/lambda;
    end
end



% creates a vector of interval durations that follow an exponential function
% (e.g. for ITI or trial duration)
% requires randExpoInt() function for the generation of single intervals
% nrIntervals = number of values to create
% minMaxDur = vector with 2 values: minimum and maximum duration
% totalDur = optional; if a value is set, the function forces all intervals
% to add up this value
% maxTries = if totalDur is set, maxTries can be set to determine how many
% attempts will be undertaken; default is 10,000. If the function fails to
% find a solution for given minMaxDur and totalDur, an error message is
% given with estimates of how much the expectancy value of the interval or
% the total duration has to be adapted.

function intDurVec = createExpIntDurs(nrIntervals, minMaxDur, totalDur, maxTries)
    if nargin < 4
        maxTries = 10000; % default for maximum tries is 10,000
    end
    if nargin < 3
        totalDur = []; % if no totalDur is given, just sample all intervals independently
    end

    % create required vectors with NaN
    intDurVec = NaN(nrIntervals,1);
    if ~isempty(totalDur)
        lastInts = NaN(maxTries,1);
        meanIntDur = NaN(maxTries,1);
    end
    
    currTry = 1; % initiate loop counter variable

    % create interval list until all values are within specified boundaries
    % (minMaxDur) and only while the maximum number of tries has not been
    % reached
    while (sum(intDurVec > minMaxDur(1)) ~= nrIntervals || ...
          sum(intDurVec < minMaxDur(2)) ~= nrIntervals) && ...
          currTry <= maxTries
        % loop through list of intervals until second to last
        for i = 1:nrIntervals-1
            intDurVec(i) = randExpoInt(minMaxDur);
        end
        if isempty(totalDur)
            % if no total duration is given, just compute the last interval
            % independently
            intDurVec(nrIntervals) = randExpoInt(minMaxDur);
        else
            % if total duration is given, the last interval is computed to
            % exactly achieve totalDur
            intDurVec(nrIntervals) = totalDur - sum(intDurVec(1:nrIntervals-1));
            meanIntDur(currTry) = mean(intDurVec(1:end-1)); % log mean interval duration
            lastInts(currTry) = intDurVec(nrIntervals); % log lengrth of last (forced) interval
        end
        currTry = currTry + 1; % loop counter +1
    end

    if currTry > maxTries % if function has not found a solution in specified number of tries
        totalDurOff = mean(lastInts) - mean(meanIntDur); % compute the average absolute discrepancy between last interval and the mean of the others
        meanDurOff = totalDurOff / nrIntervals; % compute the average mean discrepancy
        if totalDurOff < 0 % if last trials were too short
            error(['Interval not compatible with total duration. Either use ' ...
                   'interval with lower mean (ca. ' num2str(abs(meanDurOff)) ' less) or ' ...
                   'increase total duration (by ca. ' num2str(abs(totalDurOff)) ').'])
        elseif totalDurOff > 0 % if last trials were too long
            error(['Interval not compatible with total duration. Either use ' ...
                   'interval with higher mean (ca. ' num2str(meanDurOff) ' more) or ' ...
                   'decrease total duration (by ca. ' num2str(totalDurOff) ').'])
        end
    end
end



function ratingMat = rateCSmouse(window, buttonsOK, pauseDur, continuousRatings, csTextureVec)
% PTB rating routine for Likert scales with mouse in put & slider response mode;
% The routine will first acquire rating # 1 for all stimuli, then rating 2, 
% and so on; order of stimuli is randomized at beginning of rating
% turn (but constant across ratings within turn);
% item texts, scale ticks & labels can be easily adapted (see "header") within this function;
% scales everything relative to screen (except for font size)
% window: Psychtoolbox window pointer
% buttonsOK: indices of mouse buttons to select current rating
% pauseDur: fixed duration of fixation cross after each rating in seconds
% continuousRatings: false = Likert scale, only ticks can be selected; true
% = continuous ratings, i.e., participants can select values between ticks
% csTextureVec: vector with textures for all CS to be rated
% output: ratingMat = 1 x [nrRatings*3] vector with all ratings (rating1_stim1, 
% rating1_stim2 ... rating2_stim1 ... lastRating_2ndlastStim, lastRating_lastStim);
% after ratings, x coordinates of mouse cursor & rating durations (in sec) are
% logged (same order as ratings)


    %% header with parameters
    % can present any amount of different ratings. First line in every cell
    % array corresponds to first rating, second line = second rating, and
    % so on; make sure that all cell arrays in the header section have the
    % same number of lines

    textVec = {'How unpleasant is this pattern to you?'; ...
               'How arousing is this pattern to you?'; ...
               'How fearful are you seeing this pattern?'; ...
               'How likely (in %) will this patern be followed by a shock?'};

    anchorsVec = {'very pleasant', 'very unpleasant'; ...
                  'not arousing at all', 'very arousing'; ...
                  'not fearful at all', 'very fearful'; ...
                  'never followed', 'always followed'};
    
    labelVec = {0:10; ...
                0:10; ...
                0:10; ...
                0:10:100};

    if length(unique([length(textVec), length(anchorsVec), length(labelVec)])) > 1
        error('Cell arrays in header section need same number of lines')
    end
    

    textCol = [0 0 0]; % color for text & scale
    noSelCol = [255 0 0]; % color for cursor before selection is confirmed
    selMadeCol = [0 255 0]; % color for cursor after selection is confirmed

    %% matrices that we need
    orderVec = randperm(length(csTextureVec));
    ratingMat = NaN(length(csTextureVec), length(labelVec), 3);
    
    % automatic computation of sizes relative to screen
    % global values / coordinates
    [xTotal, yTotal] = Screen('WindowSize', window);
    xCenter = xTotal / 2; %yCenter = yTotal ./ 2;

    % x coordinates
    widthLeft = round(xTotal * .2, 0);
    widthRight = round(xTotal * .8, 0);
    widthTicks = cell(length(labelVec),1);
    for i = 1:length(labelVec)
        widthTicks{i} = widthLeft : (widthRight-widthLeft)/(length(labelVec{i})-1) : widthRight;
    end

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
    coordsTicks = cell(length(labelVec),1);
    for i = 1:length(labelVec)
        coordsTicks{i} = [sort(repmat(widthTicks{i}, 1, 2)); ...
                          repmat([heightScale-.05*yTotal,  heightScale+.05*yTotal], ...
                          1, length(labelVec{i}))];
    end
    

    %% Rating Loop
    for RatI = 1:length(labelVec)
        for StimI = 1:length(csTextureVec)
            
            % initiate cursor properties
            pointCol = noSelCol;
            selMade = false;
            
            startTime = GetSecs;
    
            while selMade == false
                % get mouse coordinates
                [x,~,buttons] = GetMouse(window);
                % if mouse is outside the scale, set coordinates to scale min/max
                if x < widthLeft
                    x = widthLeft;
                elseif x > widthRight
                    x = widthRight;
                end
                % either put marker on mouse coordniates or on closest scale tick
                if continuousRatings == true
                    xPointer = x;
                elseif continuousRatings == false
                    xDist = abs(x - widthTicks{RatI});
                    xPointer = min(widthTicks{RatI}(xDist == min(xDist)));
                end
                % detect clicks; in case of click, log rating and end loop
                if sum(buttons(buttonsOK) > 0)
                    relPos = (xPointer-widthLeft) / (widthRight-widthLeft);
                    loggedRating = labelVec{RatI}(1) + relPos*(labelVec{RatI}(end)-labelVec{RatI}(1));
                    pointCol = selMadeCol;
                    selMade = true;
                elseif sum(buttons(buttonsOK)) == 0
                    pointCol = noSelCol;
                end
    
                % Draw Stimulus
                Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
                
                % Draw Question
                DrawFormattedText(window, textVec{RatI}, 'center', heightQuestion, textCol);
                
                % Draw Scale
                Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
                Screen('DrawLines', window, coordsTicks{RatI}, [], textCol);
                Screen('DrawLine', window, pointCol, ...
                                   xPointer, heightScale-.05*yTotal, ...
                                   xPointer, heightScale+.05*yTotal, 3);
                
                % Draw Labels
                for lab = 1:length(labelVec{RatI})
                    DrawFormattedText(window, int2str(labelVec{RatI}(lab)), ...
                                      'center', heightLabels, textCol, ...
                                      [], [], [], [], [], ...
                                      [widthTicks{RatI}(lab), heightLabels, widthTicks{RatI}(lab), heightLabels]);
                end
                
                % Draw Anchors
                DrawFormattedText(window, anchorsVec{RatI,1}, ...
                                  'center', heightAnchors, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicks{RatI}(1), heightAnchors, widthTicks{RatI}(1), heightAnchors]);
                DrawFormattedText(window, anchorsVec{RatI,2}, ...
                                  'center', heightAnchors, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicks{RatI}(end), heightAnchors, widthTicks{RatI}(end), heightAnchors]);
                
                % Flip it!
                Screen('Flip', window)
            end
            
            ratDur = GetSecs() - startTime;
    
            % log rating
            ratingMat(orderVec(StimI), RatI, 1) = loggedRating;
            ratingMat(orderVec(StimI), RatI, 2) = xPointer;
            ratingMat(orderVec(StimI), RatI, 3) = ratDur;
            
            % show the marked selected rating for a moment
            WaitSecs(.3);
    
            % draw fixation cross and wait for [pauseDur] seconds
            presFix(window, pauseDur);
        end % rating loop
    end

    % turn 3D ratingMat into 1 x [number all values] vector (i.e., a line)
    ratingMat = ratingMat(:)';
end