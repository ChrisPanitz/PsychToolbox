function gaborgen24_fMRI_Day2(subNo, csPerm)
%% SET TO ZERO FOR ACTUAL EXPERIMENT
Screen('Preference', 'SkipSyncTests', 0);


%% Header
TRdur = 2; % TR length in seconds

% trial parameters 
nrBlocksRecall = 2; % # of blocks in recall; ratings after each block
nrTrialsRecall = 18; % # of trials per stimulus and per block

minMaxItiSec = [5.5 15.4]; % min & max ITI in sec, taken from exponential distribution
meanItiSec = 7; % mean ITI duration in sec

% How many sequential presentations of the same CS/GS are allowed
% if no restrictions wanted ==> set to ridiculously high number
maxSeqCS = 2; 

% CS parameters
rotStim = [15 35 55 75]; % vector of orientation angles (clockwise from 0 = vertical)
cycPerDeg = 3.5; % spatial frequency of Gabors in cycles/degree visual angle
imSizeAng = [5 5]; % size of Gabor patches in degrees visual angle

% distance participant <-> screen
seatDistInch = 125/2.54;

% flicker parameters
% framPerCyc for 15 Hz stim: 8 frames at 120 Hz screen or 4 frames at 60 Hz screen
framPerCyc = 4; % duration of one cycle in frames
flickMod = 'box'; % how flicker stimulus is modulated. 'box' or 'sin'
flickDurFrames = 124; % duration of whole stimulus in frames; 124/60 = 2.067 s

% US parameters - for BIOPAC stimulator
% Not needed at the time of writing
% leaving it in, in case needed in the future
%nrPulses = 8; % number of single pulses
%secBtwPulses = 1/120; % time between single pulses in sec
%shockOnsetSecs = 2; % shock onset how many seconds after CS onset?
% analog output settings to BIOPAC signal generator (not the actual current stimulator)
% only change if you know what you are doing
%baselineVolt = -10; % baseline DC throughout experiment
%peakVolt = 10; % voltage level to trigger shock
%daqScanRate = 960000; % scan rate of the D/A box - max = 1 MHz

% rating & instruction parameters
pauseBtwRat = 1*TRdur; % pause duration between two ratings, in sec
ratingDur = 2*TRdur; % time out for ratings, in sec
instructDur = 2*TRdur; % fixed duration of instruction screens, in sec
buttonsLeft = KbName('b'); % DEC key codes for moving rating cursor left
buttonsRight = KbName('y'); % DEC key codes for moving rating cursor right
buttonsOK = KbName('g'); % DEC key codes for confirm rating selection

% aesthetics
backgroundCol = 127.5; % bakcground color; 127.5 = mid gray
gratingDark = 0; % dark grating stripes = 0 = black
gratingBright = 255; % bright grating stripes = 255 = white
fontSize = 32; % well... font size

% for logfiles and rating files
logFileFolder = 'C:\Users\dinglab.UFAD\Desktop\Gaborgen24 paradigm\';
logFilePrefix = 'gaborgen24_fMRI_Day2';
logFilePrefixDay1 = 'gaborgen24_fMRI_Day1';
ratFileFolder = 'C:\Users\dinglab.UFAD\Desktop\Gaborgen24 paradigm\';
ratFilePrefix = 'gaborgen24_fMRI_Day2'; 
ratAllFileFolder = 'C:\Users\dinglab.UFAD\Desktop\Gaborgen24 paradigm\';
ratAllFilePrefix = 'gaborgen24_fMRI_Day2'; 

% strings of instructions to participants          
welcomeMsg = ['Thank you for participanting in our experiment In the upcoming\n' ...
              'task you will be presented flickering pictures\n' ...
              'in the center of the screen. Please pay attention\n' ...
              'to them. An electric shock also may occur from time\n' ...
              'to time. Every now and then you will be asked to\n' ...
              'answer questions about how you experience the pictures.\n'...
              'The whole task will last about 50 minutes.\n\n\n' ...
              'If you have questions, please let the\n' ...
              'experimenter know now. Befor the task begins, please answer' ...
              'the questions about the pictures you know from yesterday.\n\n'];
preRatMsg = ['Please rate for each picture how unpleasant and arousing you find it,\n' ...
             'how fearful you are seeing it, and how likely you think it will be\n' ...
             'followed by a shock. You can move the cursor to the left or right\n' ...
             'by using the buttons in your hand. You have ' int2str(ratingDur) ' seconds\n'  ...
              'to move the cursor around before your rating is saved.'];
postRatMsg = ['Thank you for rating the pictures.\n' ...
              'The task will continue in a moment.'];
goodbyeMsg = ['You have completed the task. Thank you for your participation!\n\n' ...
              'The experimenter will be with you in a moment.'];


% ports & codes
addressesLPTHex = '3EFC'; % LPT port adress in HEX
firstTRcode = 99; % port code for first TR
% stimCodes = 2; % codes (DEC) sent via LPT when Gabors are presented -
% stimCodes currently is set for each trial
TRtriggerCodes = ['t','5%']; % serial port codes from fMRI, coded as keyboard strokes     


%% check if there is a Day 1 logfile for this participant and if they have
% the same CS-Gabor permutation
if nargin < 3
    forceSubAndPerm = false;
end
if forceSubAndPerm == false
    logFileNameDay1 = strcat(logFileFolder, logFilePrefixDay1, '_', num2str(subNo), '_logfile.dat');
    if exist(logFileNameDay1, 'file')
        permDay1 = dlmread(logFileNameDay1, ',', [1 1 1 1]);
    else
        error('No Day 1 logfile for this participant number!');
    end
    if csPerm ~= permDay1
        error('The specified CS-Gabor permutation is different from Day 1!');
    end
end


%% clearing and initializing stuff
IOPort('CloseAll');
clc;

% Setting rng seed; "old" rand function (instead of rng) was used for
% reasons of compatibility with lab computer
rand('state',sum(100*clock));

AssertOpenGL; 

% LPT port
lptPortObj = io64;
status = io64(lptPortObj);
adressesLPT = hex2dec(addressesLPTHex);


%% Prepare logging

% files to write trial & rating data + variable names for header line
% logfiles
logFileName = strcat(logFileFolder, logFilePrefix, '_', num2str(subNo), '_logfile.dat');
fIDLog = fopen(logFileName, 'w');
fprintf(fIDLog,'partNo,csPerm,trial,phase,block,stim,paired,timeSinceFirstTR,actFlickDur,itiDur\n');
fclose(fIDLog);
% ratings
ratingsHeader = ['partInd,ratInd,' ...
                 'val_csp,val_gs1,val_gs2,val_gs3,ar_csp,ar_gs1,ar_gs,ar_gs3,' ...
                 'fear_csp,fear_gs1,fear_gs2,fear_gs3,exp_csp,exp_gs1,exp_gs2,exp_gs3,' ...
                 'valDur_csp,valDur_gs1,valDur_gs2,valDur_gs3,arDur_csp,arDur_gs1,arDur_gs,arDur_gs3,' ... 
                 'fearDur_csp,fearDur_gs1,fearDur_gs2,fearDur_gs3,expDur_csp,expDur_gs1,expDur_gs2,expDur_gs3' ...
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
trialsPerBlockRecall = nrCS*nrTrialsRecall; % how many trials per HAB block across all CS

% initialize (empty) trialMat to store CS type, US [yes/no], ITI duration, phase, & block
trialMat = NaN(nrBlocksHab*trialsPerBlockRecall, 5);
% initialize (empty) ratingAfterTrial that contains indices for every
% block's last trial
ratingAfterTrial = NaN(nrBlocksRecall,1);

% some helping count variables
trialCount = 0;
ratingCount = 0;

% for recall trials
for bl = 1:nrBlocksRecall
    tempMat = NaN(trialsPerBlockRecall,5);
    tempMat(:,1:2) = pseudoShuffleWithBooster(conds, nrTrialsRecall, 0, ...
                     0, [], maxSeqCS, 0, true);
    % ITI, phase, & block columns
    tempMat(:,3) = createExpIntDurs(trialsPerBlockRecall, minMaxItiSec, meanItiSec*trialsPerBlockRecall); 
    tempMat(:,4) = 4;
    tempMat(:,5) = bl;
    
    % write values into trialMat and ratingAfterTrial
    trialCount = trialCount + trialsPerBlockRecall;
    ratingCount = ratingCount + 1;
    trialMat((trialCount-trialsPerBlockRecall+1):trialCount,:) = tempMat;
    ratingAfterTrial(ratingCount) = trialCount;
end

% flicker sequences (on/off frames)   
flickerVec = genFlickVec(framPerCyc, flickDurFrames, flickMod);
flickerVec = squeeze(flickerVec)';

% initialize Biopac device & output port for electrical stimulation
% in case it was already initialized 
%daqreset; 
% create object for data acquisition toolbox device, add output channel, 
% set scan rate, and set DC to baseline level
%d = daq('mcc');
%addoutput(d,'Board0','Ao0','Voltage');
%d.Rate = daqScanRate;
%write(d, baselineVolt); 
% create square function for shock delivery
%biopacTrigger = [baselineVolt*ones(round(shockOnsetSecs*d.Rate),1); ...
%                 repmat([peakVolt.*ones(round(secBtwPulses*d.Rate/2),1); ... 
%                        baselineVolt.*ones(round(secBtwPulses*d.Rate/2),1)], ...
%                        nrPulses, 1)]; 
%biopacNothing = [baselineVolt*ones(round(shockOnsetSecs*d.Rate),1); ...
%                 repmat(baselineVolt.*ones(round(secBtwPulses*d.Rate),1), ...
%                        nrPulses, 1)]; 


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
        
    % wait for the first MRI pulse, log system time and send code to EEG
    timeFirstTR = waitForScanTriggerKb(TRtriggerCodes);
    io64(lptPortObj, adressesLPT, firstTRcode);
    
    % welcome screen, terminated after fixed duration
    DrawFormattedText(w, welcomeMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);
    presFix(w, TRdur);

    % instruction screen & ratings, controlled by participant mouse click
    DrawFormattedText(w, preRatMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);
    presFix(w, TRDur);

    % run rating routine and write output into file (append mode)
    currRatings = rateCSMRI(w, buttonsLeft, buttonsRight, buttonsOK, ...
                            ratingDur, pauseBtwRat, textureVec);
    dlmwrite(ratFileName, [subNo, currentRatingInd, currRatings], '-append');
    dlmwrite(ratAllFileName, [subNo, currentRatingInd, currRatings], '-append');
    currentRatingInd = currentRatingInd+1;

    % another wait-for-TR to compensate for timing errors
    waitForScanTriggerKb(TRtriggerCodes);

    % instructions after ratings
    DrawFormattedText(w, postRatMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);
    
    % fix cross for the BOLD response to play out before ratings
    presFix(w, 3*TRdur);
    
    %% trials start

    for trial = 1 : size(trialMat,1)
        % create two/three-digit code for LPT: 
        % first digit: US yes (1) or no (-)... for now always 0/-
        % second digit: phase (4 = Recall)
        % third digit = stimulus (1 = CS, 2 = GS1, 3 = GS2, 4 = GS3)
        stimCodes = trialMat(trial,2)*100 + trialMat(trial,4)*10 + trialMat(trial,1);
        % present CS/GS, log stimulus duration and onset
        if trialMat(trial, 2) == 1 % if there is a US to presented (right now: never)
            [actFlickDur, onsetTime] = presFlickBiopacShock(w, ...
                                       flickerVec, [], textureVec(trialMat(trial,1)), ...
                                       biopacTrigger, d, ...
                                       lptPortObj, adressesLPT, stimCodes);
        else % no US
            [actFlickDur, onsetTime] = presFlickBiopacShock(w, ...
                                       flickerVec, [], textureVec(trialMat(trial,1)), ...
                                       biopacNothing, d, ...
                                       lptPortObj, adressesLPT, stimCodes);
        end % US vs NoUS conditional

        % present trial while logging actual flicker duration
        presFix(w, trialMat(trial,3));
        
        % write trial parameters to logfile
        % part ID, CS+ permutation, trial #, phase, block, stimulus type, paired [yes/no], 
        % stimulus onset (since first TR), simulus duration, iti duration;
        timeSinceFirstTR = onsetTime - timeFirstTR;
        trialOutVec = [subNo, csPerm, trial, trialMat(trial,4), trialMat(trial,5), ...
                       trialMat(trial,1), trialMat(trial,2), timeSinceFirstTR, ...
                       actFlickDur, trialMat(trial,3)];
        dlmwrite(logFileName, trialOutVec, '-append');

        % check if this is the last trial of a block (==> ratings)
        if ismember(trial, ratingAfterTrial)
            % wait for TR pulse to compensate for timing errors across participants
            waitForScanTriggerKb(TRtriggerCodes);
            
            % instructions for ratings
            DrawFormattedText(w, preRatMsg, 'center', 'center');
            Screen('Flip', w);
            WaitSecs(instructDur);

            % run rating routine and write output into file (append mode)
            currRatings = rateCSMRI(w, buttonsLeft, buttonsRight, buttonsOK, ...
                                    ratingDur, pauseBtwRat, textureVec);
            dlmwrite(ratFileName, [subNo, currentRatingInd, currRatings], '-append');
            dlmwrite(ratAllFileName, [subNo, currentRatingInd, currRatings], '-append');
            currentRatingInd = currentRatingInd+1;

            % another wait-for-TR to compensate for timing errors
            waitForScanTriggerKb(TRtriggerCodes);

            if trial < size(trialMat,1) % skip this at the very end of the paradigm
                % instructions after ratings
                DrawFormattedText(w, postRatMsg, 'center', 'center');
                Screen('Flip', w);
                WaitSecs(instructDur);
    
                % fix cross for the BOLD response to play out before ratings
                presFix(w, 3*TRdur);
            end
        end % rating conditional
        
    end % trial loop
        
    % final screen, disappears automatically
    DrawFormattedText(w, goodbyeMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);

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



function ratingMat = rateCSMRI(window, buttonsLeft, buttonsRight, buttonsOK, maxDurInSec, pauseDur, csTextureVec)
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

    textVec = {'How unpleasant was this pattern to you?'; ...
               'How arousing was this pattern to you?'; ...
               'How fearful did you feel when you saw this pattern?'; ...
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
    ratingMat = NaN(length(csTextureVec), length(labelVec), 2);
    
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
        % selection index variables
        selMin = 1;
        selMax = length(labelVec{RatI});
        selCent = ceil((selMax+selMin)/2);
        currTime = 0; % initialize currTime with some value that is smaller than any GetSecs
        
        for StimI = 1:length(csTextureVec)
            % initiate cursor properties
            currSelect = selCent;
            pointCol = noSelCol;
            
            startTime = GetSecs;
    
            while (currTime-startTime) < maxDurInSec
   
                % wait for button press, change selected rating value & cursor color
                [~, ~, buttonsPressed, ~] = KbCheck(-1);
                if sum(buttonsPressed(buttonsLeft)) > 0 && currSelect > selMin
                    currSelect = currSelect - 1;
                    pointCol = noSelCol;
                    WaitSecs(.100);
                elseif sum(buttonsPressed(buttonsRight)) > 0 && currSelect < selMax
                    currSelect = currSelect + 1;
                    pointCol = noSelCol;
                    WaitSecs(.100);
                elseif sum(buttonsPressed(buttonsOK)) > 0
                    pointCol = selMadeCol;
                end
                %when time is as good as up (within 20 ms), set cursor color to "selected"
                if (currTime-startTime) > (maxDurInSec-.020)
                    pointCol = selMadeCol;
                end

                % Draw Stimulus
                Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
                
                % Draw Question
                DrawFormattedText(window, textVec{RatI}, 'center', heightQuestion, textCol);
                
                % Draw Scale
                Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
                Screen('DrawLines', window, coordsTicks{RatI}, [], textCol);
                Screen('DrawLine', window, pointCol, ...
                                   widthTicks{RatI}(currSelect), heightScale-.05*yTotal, ...
                                   widthTicks{RatI}(currSelect), heightScale+.05*yTotal, 3);
                
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
            
                % log system time
                currTime = GetSecs;                
            end
            
            ratDur = GetSecs() - startTime;
    
            % log rating
            ratingMat(orderVec(StimI), RatI, 1) = labelVec{RatI}(currSelect);
            ratingMat(orderVec(StimI), RatI, 2) = ratDur;
    
            % show selected rating for 250 ms, draw fixation cross and show it for the remainder of [pauseDur] seconds
            % if [pauseDur] is shorter than 250 ms, selection will be shown for this amount of time, no fixation cross
            if pauseDur > .250
                WaitSecs(.250);
                presFix(window, pauseDur-.250);
            else
                WaitSecs(pauseDur);
            end
        end % rating loop
    end

    % turn 3D ratingMat into 1 x [number all values] vector (i.e., a line)
    ratingMat = ratingMat(:)';
end