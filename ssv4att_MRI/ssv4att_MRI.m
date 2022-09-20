function ssv4att_MRI(subNo)
%% SET TO ZERO FOR ACTUAL EXPERIMENT
Screen('Preference', 'SkipSyncTests', 1);


%% Header
% for logfiles
logFileFolder = 'C:\Users\dinglab.UFAD\Desktop\ssv4att_MRI\logfiles\';
logFilePrefix = 'ssv4att_MRI'; 

% scanner parameters
TRdur = 2; % TR length in seconds

% trial parameters 
breakAftTr = 76; % breaks after how many trials?
picsPerArray = 4; % there will be a colorpic and a grayscalepic array
shiftPix = 5; % how many pixels pictures shift per frame on the x and y axes, respectively 
trNoEvents = 116; % number of trials with no events, must be a multiple of 4
trEventsBL = 12; % number of trials with only baseline events, must be a multiple of 4
trEventsCue = 12; % number of trials with only cue task events, must be a multiple of 6
trEventsBoth = 12; % number of trials with events in both baseline and cue task, must be a multiple of 6

% flicker parameters
% frPerCyc for 120Hz screen: 14 = 8.57 Hz, 8 = 15 Hz
frPerCyc = [14 8]; % duration of one cycle in frames (for different frequencies)
flickMods = {'box'}; % how flicker stimulus is modulated. 'box' or 'sin' or both possible
frPerMetaCyc = 56; % meta cycle is the time window after which all frequencies are in phase again (here in frames)
metaCyclesTotal = 13; % trial length in meta cycles
metaCyclesBaseline = 3:5; % how many metacycle windows can the baseline last

% baseline part
blEventBuffer = 36; % in frames; 36 frames = 300 ms
blEventDur = 24; % in frames; 24 frames = 200 ms
%maxBlEventsPerMC = 1; % max amount of baseline events per metacycle
maxBlEventsPerTrial = 2; % max amount of baseline events per trial
standardTickLength = 50; % in pixels
eventTickRel = 1.10; % scaling factor for event
centTileCol = [127 127 127]; % color of the central tile

% cued attention part
cueEventBuffer = 96; % in frames; 96 frames = 800 ms
cueEventDur = 72; % in frames; 72 frames = 600 ms
maxCueEventsPerMC = 1/3; % max amount of cue task events per metacycle
%maxCueEventsPerTrial = 3; % max amount of cue task events per trial
standardLum = 1; % alpha value for flickering
eventLum = .25; % alpha value during cue task event
cueLetters = ['C','G']; % cues that are displayed in the center
textForRatings = {'color', 'gray'}; % for the ratings after trials

% ITI duration parameters
meanITI = 8; % in seconds
minMaxITItotal = [6.5 16.4]; % will be drawn fro man exponential distribution
minMaxITIprerat = [0.5 1.5]; % time between trial offset and rating, fixation cross
pauseDur = 20; % pause duration between blocks in seconds

% some aesthetics
backgroundCol = 0; 
textCol = 255;
fontSize = 36;

% rating & instruction parameters
ratingDur = 2*TRdur; % time out for ratings, in sec
scaleLabels = {0:6}; % response options
scaleAnchors = {' ',' '}; % don't want labels but function expects input
%ratingText = ... rating text will be set in trial loop
pauseBtwRat = 0.5; % pause after rating in seconds
instructDur = 3*TRdur; % fixed duration of instruction screens, in sec
buttonsLeft = KbName('b'); % DEC key codes for moving rating cursor left
buttonsRight = KbName('y'); % DEC key codes for moving rating cursor right
buttonsOK = KbName('g'); % DEC key codes for confirm rating selection

% image parameters
fieldSizeAng = 10.2; % diameter of circular field in which pictures are presented
imSizeAng = [1.5 1.5]; % picture size in degrees of visual angle
% distance participant <-> screen
seatDistInch = 70; 
% folder with images
imFolder = 'C:\Users\dinglab.UFAD\Desktop\ssv4att_MRI\imageStimuli\'; % path for image files
% lists of filenames for images - must have same length
imFileListColor = {'pic01c.jpg', 'pic02c.jpg', 'pic03c.jpg', 'pic04c.jpg', ...
                   'pic05c.jpg', 'pic06c.jpg', 'pic07c.jpg', 'pic08c.jpg', ...
                   'pic09c.jpg', 'pic10c.jpg', 'pic11c.jpg', 'pic12c.jpg', ...
                   'pic13c.jpg', 'pic14c.jpg', 'pic15c.jpg', 'pic16c.jpg', ...
                   'pic17c.jpg', 'pic18c.jpg', 'pic19c.jpg', 'pic20c.jpg'};
imFileListGray = {'pic01g.jpg', 'pic02g.jpg', 'pic03g.jpg', 'pic04g.jpg', ...
                  'pic05g.jpg', 'pic06g.jpg', 'pic07g.jpg', 'pic08g.jpg', ...
                  'pic09g.jpg', 'pic10g.jpg', 'pic11g.jpg', 'pic12g.jpg', ...
                  'pic13g.jpg', 'pic14g.jpg', 'pic15g.jpg', 'pic16g.jpg', ...
                  'pic17g.jpg', 'pic18g.jpg', 'pic19g.jpg', 'pic20g.jpg'};

% strings of instructions to participants          
welcomeMsg = ['Thank you for participanting in our experiment. We now start the task.\n\n' ...
              'Remember to count the events in the central cross and\n' ...
              'in the attended pictures. In each trial, there might be\n' ...
              'one, several, or no events at all.\n\n' ...
              'Spread your attention across the screen while you keep your\n' ...
              'eyes fixed on the central gray square. We recommend you to\n' ...
              'blink as little as possible while the flashing pictures are on.\n\n' ...
              'The complete task will take about 35 minutes.\n'...
              'There will a break halfway through the task.\n' ...
              'When you are ready, continue by pressing any button.'];
breakMsg = ['You can take a short break now.\n\n' ...
            'The task will continue after 20 seconds.\n'];
goodbyeMsg = ['You have completed the task. The scan is still in progress.\n\n' ...
              'Please remain still. The experimenter will be with you after the scan ends.'];

% ports & codes
addressLPTHex = '3EFC'; % LPT port adress in HEX
firstTRcode = 99; % port code for first TR
TRtriggerCodes = 't'; % serial port codes from fMRI, coded as keyboard strokes              

%% clearing and initializing stuff
Screen('CloseAll');
IOPort('CloseAll');
clc;

% Setting rng seed; "old" rand function (instead of rng) was used for
% reasons of compatibility with lab computer
rand('state',sum(100*clock));

AssertOpenGL; 

% LPT port
lptPortObj = io64;
status = io64(lptPortObj);
addressLPT = hex2dec(addressLPTHex);

%% Prepare logging
% check if logfile already exists
logFileName = strcat(logFileFolder, logFilePrefix, '_', num2str(subNo), '.dat');
if exist(logFileName, 'file')
    error(['There already is a logfile for this participant number. Please use different number or remove logfile from this folder: ' logFileFolder]);
end
% write header line of logfile
fID = fopen(logFileName, 'w');
fprintf(fID,['subNo,trial,attendCond,freqCol,freqAttend,blInFr,' ...
             'nrBLEvents,nrCueEventsCol,nrCueEventsGray,nrCueEventsAttend,nrCueEventsNonattend,eventsCounted,' ...
             'timeSinceFirstTR,itiDur,preRatDur,ratDur,actFlickDur,' ...
             'pCol1,pCol2,pCol3,pCol4,pGray1,pGray2,pGray3,pGray4\n']);
fclose(fID);

%% actual experiment is in try loop
try
    %% get all info for screen & window
    Screen('CloseAll');
    screens = Screen('Screens');
    screenNumber = max(screens);
    w = Screen('OpenWindow', screenNumber, backgroundCol);
    [winX,winY] = Screen('WindowSize',w);
    centX = winX/2; centY = winY/2;
        
    % set alpha blending mode for manipulating picture transparency
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % Set Font Size
    Screen('TextSize',w,fontSize);
    
    % Loading Screen
    DrawFormattedText(w, '[Preparing Stimuli...]', 'center', 'center', textCol);
    Screen('Flip', w);

    % Compute screen's PPI via window size (pix) and display size (mm)
    ppi = Screen('WindowSize', w) / Screen('DisplaySize', screenNumber) * 25.4;
    % translate visual angle into pixels for picture scaling
    fieldSizePix = pixFromAngle(fieldSizeAng, seatDistInch, ppi);
    imSizePix = pixFromAngle(imSizeAng, seatDistInch, ppi);

    % x and y coordinates of the central "fixation" tile
    centTileCoord = [centX-imSizePix(1)/3, centY-imSizePix(2)/3, ...
                     centX+imSizePix(1)/3, centY+imSizePix(2)/3];
    
    % radius for picture tiles and central tile
    picRad = sqrt(sum(imSizePix.^2)) / 2;
    centTileRad = sqrt(sum((imSizePix*2/3).^2)) / 2;

    % define allowed distances of pic centres from screen center (as [0,0])
    minRho = ceil(centTileRad + picRad);
    maxRho = floor(fieldSizePix/2 - picRad);

           
    %% create vectors & matrices for stimuli & trials
    % some intermediate computations
    flickDurFrames = frPerMetaCyc*metaCyclesTotal; % duration of whole stimulus in frames
    flickDurBaseline = metaCyclesBaseline*frPerMetaCyc; % possible durations of baseline in frames
    nrTrials = trNoEvents+trEventsBL+trEventsCue+trEventsBoth; % total number of trials
    
    % blLengthVec: length of baseline in frames
    blLengthVec = repmat(flickDurBaseline',ceil(nrTrials/length(flickDurBaseline)),1);
    
    %trialMat has all the trial info
    trialMat = zeros(nrTrials,11);
    
    % which cue (1/2)... see cueLetters for assignment
    % here: 1 = C, 2 = G
    trialMat(1:trNoEvents+trEventsBL,1) = repmat([1,1,2,2],1,(trNoEvents+trEventsBL)/4); 
    trialMat(trNoEvents+trEventsBL+1:end,1) = repmat([1,1,1,2,2,2],1,(trEventsCue+trEventsBoth)/6); 
    % driving frequency of color array (depends on frPerCyc)
    % here, 1: Color = 8.57 Hz, 2: Color = 15 Hz
    trialMat(:,2) = repmat([1,2],1,nrTrials/2);
        trialMat(trNoEvents+trEventsBL+1:6:end,2) = 2;
        trialMat(trNoEvents+trEventsBL+2:6:end,2) = 1;
        trialMat(trNoEvents+trEventsBL+3:6:end,2) = 2;
        trialMat(trNoEvents+trEventsBL+4:6:end,2) = 1;
        trialMat(trNoEvents+trEventsBL+5:6:end,2) = 2;
        trialMat(trNoEvents+trEventsBL+6:6:end,2) = 1;
    % recoding: driving frequency of Attended Pics (depends on frPerCyc)
    % here, 1 = attended at 8.57 Hz, 2 = attended at 15 Hz
    trialMat(:,3) = 2 - (trialMat(:,1) == trialMat(:,2)); 
    % length of baseline in frames
    trialMat(:,4) = blLengthVec(randperm(size(trialMat,1)));
    for i = trNoEvents+1 : trNoEvents+trEventsBL
        %trialMat(i,5) = randi([1, trialMat(i,4)/frPerMetaCyc*maxBlEventsPerMC]);
        trialMat(i,5) = randi([1, maxBlEventsPerTrial]);
    end
    % number of events in first array (here color col 6) and in second
    % array (here gray, col 7)
    trialMat(trNoEvents+trEventsBL+1:trNoEvents+trEventsBL+trEventsCue+trEventsBoth, 6:7) = ...
        repmat([1,0; 0,1; 1,1], (trEventsCue+trEventsBoth)/3, 1);
    for i = trNoEvents+trEventsBL+1 : trNoEvents+trEventsBL+trEventsCue
        maxCues = ceil((metaCyclesTotal - trialMat(i,4)/frPerMetaCyc) * maxCueEventsPerMC) - 1;
        trialMat(i,6:7) = trialMat(i,6:7) .* randi([1, maxCues], 1, 2);
        %trialMat(i,6:7) = trialMat(i,6:7) .* randi([1, (metaCyclesTotal - trialMat(i,4)/frPerMetaCyc) * maxCueEventsPerMC], 1, 2);
        %trialMat(i,6:7) = trialMat(i,6:7) .* randi([1, maxCueEventsPerTrial], 1, 2);
    end
    for i = trNoEvents+trEventsBL+trEventsCue+1 : trNoEvents+trEventsBL+trEventsCue+trEventsBoth
        %trialMat(i,5) = randi([1, trialMat(i,4)/frPerMetaCyc*maxBlEventsPerMC]);
        trialMat(i,5) = randi([1, maxBlEventsPerTrial]);
        maxCues = ceil((metaCyclesTotal - trialMat(i,4)/frPerMetaCyc) * maxCueEventsPerMC) - 1;
        trialMat(i,6:7) = trialMat(i,6:7) .* randi([1, maxCues], 1, 2);
        %trialMat(i,6:7) = trialMat(i,6:7) .* randi([1, (metaCyclesTotal - trialMat(i,4)/frPerMetaCyc) * maxCueEventsPerMC], 1, 2);
        %trialMat(i,6:7) = trialMat(i,6:7) .* randi([1, maxCueEventsPerTrial], 1, 2);
    end
    % recode number of events: attended array (col 8) and non-attended
    % array (col 9)
    trialMat(trialMat(:,1)==1, 8:9) = trialMat(trialMat(:,1)==1, 6:7);
    trialMat(trialMat(:,1)==2, 8:9) = trialMat(trialMat(:,1)==2, 7:-1:6);
    % total duration of ITI (including rating and fixation cross)
    trialMat(:,10) = createExpIntDurs(size(trialMat,1),minMaxITItotal,size(trialMat,1)*meanITI);
    % duartion of the ITI part between flicker offset and rating
    trialMat(:,11) = minMaxITIprerat(1) + rand(size(trialMat,1),1).*diff(minMaxITIprerat);
    trialMat = trialMat(randperm(size(trialMat,1)),:);
    
    % compute length and coordinates for ticks of central cross (for each frame)
    xyTick = pol2cart(0.25*pi,standardTickLength/2);
    xyEvent = pol2cart(0.25*pi,standardTickLength/2*eventTickRel);
    standardTickPos = [centX-xyTick, centX+xyTick, centX-xyTick, centX+xyTick; ...
                       centY-xyTick, centY+xyTick, centY+xyTick, centY-xyTick];
    eventTickPos = [centX-xyEvent, centX+xyEvent, centX-xyEvent, centX+xyEvent; ...
                    centY-xyEvent, centY+xyEvent, centY+xyEvent, centY-xyEvent];

    % select for each trial the frames during which events will be "active" 
    blEventFrames = prepareBlEventFrames(trialMat(:,5),trialMat(:,4),flickDurFrames,blEventBuffer,blEventDur);
    cueEventFrames = prepareCueEventFrames(trialMat(:,6),trialMat(:,7),trialMat(:,4),flickDurFrames,cueEventBuffer,cueEventDur,picsPerArray);
   
    % flicker sequences (on/off frames) for different tagging frequencies
    flickerVecs = NaN(length(frPerCyc), length(flickMods), flickDurFrames);
    for freq = 1:length(frPerCyc)
        for modul = 1:length(flickMods)
            flickerVecs(freq,modul,:) = genFlickVec(frPerCyc(freq), flickDurFrames, flickMods{modul});
        end % modul loop
    end % freq loop
    flickerVecs = squeeze(flickerVecs.*standardLum); % results in frequencies x samples
    relEventLum = eventLum/standardLum; % compute factor to multiply standardLum vector with later
    
    %% load and prepare images    
    imageMat = cell(length(imFileListColor),2);
    textureMat = NaN(length(imFileListColor),2);
    for im = 1:length(imFileListColor)
        imageMat{im,1} = imread([imFolder imFileListColor{im}]);
        textureMat(im,1) = Screen('MakeTexture', w, imageMat{im,1});
        imageMat{im,2} = imread([imFolder imFileListGray{im}]);
        textureMat(im,2) = Screen('MakeTexture', w, imageMat{im,2});
    end % filling imageMat and textureVec (im loop)

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   HERE COMES THE ACTUAL EXPERIMENT   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    io64(lptPortObj, addressLPT, 0); % just to be sure
    HideCursor();
    
    % welcome screen, terminated after participant button press
    DrawFormattedText(w, welcomeMsg, 'center', 'center');
    Screen('Flip', w);
    waitForScanTriggerKb([buttonsLeft, buttonsRight, buttonsOK]); % actually waits for any button press
        
    % wait for the first MRI pulse, log system time and send code to EEG
    Screen('Flip', w);
    timeFirstTR = waitForScanTriggerKb(TRtriggerCodes);
    io64(lptPortObj, addressLPT, firstTRcode);
    WaitSecs(.010);
    io64(lptPortObj, addressLPT, 0);
    presFix(w, 3*TRdur, textCol);
    
    
    %% trials start
    for trial = 1 : size(trialMat,1)
        % select pictures
        picListInd = randperm(length(imFileListColor));
        textureVec = [textureMat(picListInd(1:picsPerArray),1); ...
                       textureMat(picListInd(picsPerArray+1:2*picsPerArray),2)];

        % prepare vectors for central tile
        tickLengthVec = zeros(2,4,flickDurFrames);
        tickLengthVec(:,:,1:trialMat(trial,4)) = repmat(standardTickPos,1,1,trialMat(trial,4));
        tickLengthVec(:,:,logical(blEventFrames{trial,:})) = repmat(eventTickPos,1,1,sum(blEventFrames{trial,:}));
        cueVec(1,1:trialMat(trial,4)) = ' ';
        cueVec(1,trialMat(trial,4)+1:flickDurFrames) = cueLetters(trialMat(trial,1));
        
        % prepare luminance vectors
        lumVectors = NaN(2*picsPerArray, flickDurFrames);
        lumVectors(1:picsPerArray,:) = repmat(flickerVecs(trialMat(trial,2),:), picsPerArray,1);
        lumVectors(picsPerArray+1:end,:) = repmat(flickerVecs(3-trialMat(trial,2),:), picsPerArray,1);
        lumVectors(logical(cueEventFrames{trial})) = lumVectors(logical(cueEventFrames{trial})) .* relEventLum;
        
        % compute positions for pictures
        posVecCenter = createCircPositions (2*picsPerArray, flickDurFrames, ...
                                      minRho, maxRho, 2*picRad, shiftPix);
        posVecCenter = permute(posVecCenter,[2 1 3]);
        posVecCoord = cat(1, ...
                        centX + posVecCenter(1,:,:) - imSizePix(1)/2, ...
                        centY + posVecCenter(2,:,:) - imSizePix(2)/2, ...
                        centX + posVecCenter(1,:,:) + imSizePix(1)/2, ...
                        centY + posVecCenter(2,:,:) + imSizePix(2)/2);
        
        % set text for rating
        ratingText = {['How many events in the central cross and the flashing ' ...
                       textForRatings{trialMat(trial,1)} ' pictures\n' ...
                       'did you count in total?']};

        % marker code to send to EEG data
        % tens: 1 = attend to color; 2 = attend to gray
        % ones/units: 1 = attend to 8.57 Hz; 2 = attend to 15 Hz
        % i.e.: 11 = attend to color @ 8.57 Hz; 12 = attend to color @ 15 Hz
        % 21 = attend to gray @ 8.57 Hz; 22 = attend to gray @ 15 Hz
        lptCode = 10*trialMat(trial,1) + 1*trialMat(trial,3);

        % present trial while logging actual flicker duration
        [actFlickDur, onsetTime] = presTrial_ssv4att_MRI(w, ...
                                   textureVec, posVecCoord, lumVectors, trialMat(trial,4), ...
                                   centTileCol, centTileCoord, tickLengthVec, cueVec, ...
                                   lptPortObj, addressLPT, lptCode);

        %fixation
        presFix(w, trialMat(trial,11), textCol);

        % rating
        ratingMat = rateMRInopic(w, buttonsLeft, buttonsRight, buttonsOK, ...
                                 'lowest', pauseBtwRat, false, ratingDur, ...
                                 scaleLabels, scaleAnchors, ratingText, textCol);

        % fixation
        presFix(w, trialMat(trial,10)-trialMat(trial,11)-ratingDur-pauseBtwRat, textCol);

        % write parameters to data file
        timeSinceFirstTR = onsetTime - timeFirstTR;
        trialOutVec = [subNo trial trialMat(trial,1:9) ratingMat(1,1) ...
                       timeSinceFirstTR trialMat(trial,10:11) ratingMat(1,2) actFlickDur ...
                       picListInd(1:2*picsPerArray)];
        dlmwrite(logFileName, trialOutVec, '-append');
        
        % Short break for participants (after [pauseAftTr] trials but not after last one)
        if trial/breakAftTr == round(trial/breakAftTr, 0) && trial ~= size(trialMat, 1)
           DrawFormattedText(w, breakMsg, 'center', 'center');
           Screen('Flip', w);
           WaitSecs(pauseDur);
           presFix(w, 3*TRdur, textCol);
        end % break conditional
    end % trial loop
        
    % final screen, disappears automatically
    DrawFormattedText(w, goodbyeMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);
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

end % ssv4 function

%% custom functions

% computes # of pixels required for picture based on the required visual
% angle, seating distance in inches, and the PPI of the screen; can take
% multiple and compute vector of values in one call
function nrPix = pixFromAngle(visAngle, distInch, ppi)
    nrPix = round(tan(visAngle/360*pi) * ppi * 2 * distInch, 0);
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



function startTime = waitForScanTriggerKb(keys,devNum)
% Written to wait for and detect MRI triggers sent via serial port; uses
% keyboard logging commands. Logs system time when trigger is detected.
% keys = vector of acceptable keys as chars (e.g. ['a', 'b']) or ints ([4, 5])
% devNum = keyboard device number; optional, defaults to -1 (all connected keyboards)
    if nargin < 2
        devNum = -1;
    end

    % create keyCode object with 256 zeros (length of KbCheck keyCode output)
    keyCode = zeros(1,256);
    
    if ischar(keys)
        % translate key strings into DEC format
        keyInd = NaN(length(keys),1);
        for key = 1:length(keys)
            keyInd(key) = KbName(keys(key));
        end
    else
        % leave integers as is
        keyInd = keys;
    end
    
    % check for "key strokes" and log them in keyCode
    while sum(keyCode(keyInd)) == 0
        [~, ~, keyCode] = KbCheck(devNum);
        WaitSecs(0.001);
    end

    % When key stroke (aka trigger) has been detected, log system time
    startTime = GetSecs;

end % function



% present fixation cross; cross size is hardcoded and relative to window
% size (pixels), duration is given in seconds
function presFix(window, fixDurSec, crossColor)
    if nargin < 3
        crossColor = 0; % default color is black
    end

    [winSize(1), winSize(2)] = Screen('WindowSize', window);
    winCenter = winSize ./ 2;
    
    % hard-coded length of fixation cross arme relative to screen height
    fixL = round(winSize(2)*.02, 0);
    
    % Present fixation cross
    Screen('DrawLines', window, [-fixL fixL 0 0; 0 0 -fixL fixL], 7, crossColor, winCenter);
    Screen('Flip', window);
    WaitSecs(fixDurSec);
end



function [positions] = createCircPositions(numPos, nrSamples, minRho, maxRho, minDistAllowed, maxPixDelta)
% creates x-y coordinates with minimal distances between each other and in 
% a circlular field with inside/outside restrictions 
% numPos = number of coordinate pairs
% nrSamples = number of samples in the time series
% minRho = minimal distance allowed from center ([0,0])
% maxRho = maximal distance allowed from center ([0,0])
% minDistAllowed = minimal direct distance between coordinate pairs
% maxPixDelta = step size in pixels between two samples on x and y axis,
% can be -maxPixDelta, 0, or maxPixDelta
    positions = NaN(numPos,2,nrSamples);
    maxNrTries = 1000;
    nrPosOkay = 0;
    
    while nrPosOkay < numPos
        [positions(1,1,1),positions(1,2,1)] = pol2cart(2*pi*rand, minRho+(maxRho-minRho)*rand);
        nrPosOkay = 1;
        for posI = 2:numPos
            minDist = 0;
            nrTries = 0;
            while minDist < minDistAllowed && nrTries < maxNrTries+1
                % randomly draw position within the min/max radius restriction
                [positions(posI,1,1),positions(posI,2,1)] = pol2cart(2*pi*rand, minRho+(maxRho-minRho)*rand);
                % compute distance to closest existion position (Pythagoras)
                minDist = min(sqrt(sum((positions(posI,:,1) - positions(1:posI-1,:,1)).^2, 2)));
                nrTries = nrTries + 1;
            end % minDist loop
            if minDist >= minDistAllowed
                nrPosOkay = nrPosOkay + 1;
            end
        end % posI loop
    end % startPosOkay loop

    for sampI = 2:nrSamples
        allRhoDist = 0;
        % loop until all distances are at least the allowed minimum
        while min(allRhoDist) < minDistAllowed
            % add x and y shift to previous positions
            positions(:,:,sampI) = positions(:,:,sampI-1) + maxPixDelta*randi([-1,1],numPos,2);
            % compute distances from center 
            [~, rho] = cart2pol(positions(:,1,sampI), positions(:,2,sampI));
            % positions inside fixation are pushed outside
            positions(rho < minRho,:,sampI) = positions(rho < minRho,:,sampI) + ...
                positions(rho < minRho,:,sampI)./abs(positions(rho < minRho,:,sampI))*maxPixDelta;
            % positions outside stimulus field are pushed inside
            positions(rho > maxRho,:,sampI) = positions(rho > maxRho,:,sampI) - ...
                positions(rho > maxRho,:,sampI)./abs(positions(rho > maxRho,:,sampI))*maxPixDelta;
            % compute all distances between positions (in pixels on x and y) 
            allXDist = positions(:,1,sampI) - positions(:,1,sampI)';
            allXDist = allXDist(logical(1-eye(numPos)));
            allYDist = positions(:,2,sampI) - positions(:,2,sampI)';
            allYDist = allYDist(logical(1-eye(numPos)));
            % compute direct distance
            [~,allRhoDist] = cart2pol(allXDist,allYDist);
        end
    end

end



function blEventFrames = prepareBlEventFrames(nrOfEventsVec,blLengthVec,durInFrames,blEventBuffer,blEventDur)
% creates on/off vectors the length of a trial in frames for each trial
% nrOfEventsVec = vector with a number of events to be presented for each trail
% blLengthVec = vector with length of baseline in frames for each trial
% durInFrames = duration of trial in frames (currently a single value)
% blEventBuffer = number of frames in the beginning and after previous events when no events can take place 
% blEventDur = duration of an event in frames
    blEventFrames = cell(length(nrOfEventsVec),1);
    for trialI = 1:length(blEventFrames)
        if nrOfEventsVec(trialI) > 0
            eventsValid = false;
            while eventsValid == false
                blEventFrames{trialI} = zeros(1,durInFrames);
                bufferCheck = zeros(1,durInFrames);
                eventLats = randi([blEventBuffer+1, blLengthVec(trialI)-blEventDur], nrOfEventsVec(trialI), 1);
                for latI = 1:length(eventLats)
                    blEventFrames{trialI}(eventLats(latI):eventLats(latI)+blEventDur) = 1;
                    bufferCheck(eventLats(latI):eventLats(latI)+blEventDur) = ...
                        bufferCheck(eventLats(latI):eventLats(latI)+blEventDur)+1;
                end
                if sum(bufferCheck > 1) == 0
                    eventsValid = true;
                end
            end
        end
    end
end



function cueEventFrames = prepareCueEventFrames(nrOfEventsVecColor,nrOfEventsVecGray,blLengthVec,durInFrames,cueEventBuffer,cueEventDur,picsPerArray)
% creates on/off vectors the length of a trial in frames for each trial
% nrOfEventsVecColor = vector with a number of events in the color pic array to be presented for each trail
% nrOfEventsVecGray = vector with a number of events in the gray pic array to be presented for each trail
% blLengthVec = vector with length of baseline in frames for each trial
% durInFrames = duration of trial in frames (currently a single value)
% cueEventBuffer = number of frames in the beginning and after previous events when no events can take place 
% cueEventDur = duration of an event in frames
% picsPerArray = number of pics per array
    cueEventFrames = cell(length(nrOfEventsVecColor),1);
    for trialI = 1:length(cueEventFrames)
        if nrOfEventsVecColor(trialI)+nrOfEventsVecGray(trialI) > 0
            cueEventFrames{trialI} = zeros(2*picsPerArray,durInFrames);
            % Color
            eventsValid = false;
            while eventsValid == false
                cueEventFrames{trialI}(1:picsPerArray,:) = 0;
                bufferCheck = zeros(1,durInFrames+cueEventBuffer);
                eventLats = randi([blLengthVec(trialI)+cueEventBuffer+1, durInFrames-cueEventDur], nrOfEventsVecColor(trialI), 1);
                for latI = 1:length(eventLats)
                    cueEventFrames{trialI}(randi([1, picsPerArray]),eventLats(latI):eventLats(latI)+cueEventDur) = 1;
                    bufferCheck(eventLats(latI):eventLats(latI)+cueEventDur+cueEventBuffer) = ...
                        bufferCheck(eventLats(latI):eventLats(latI)+cueEventDur+cueEventBuffer)+1;
                end
                if sum(bufferCheck > 1) == 0
                    eventsValid = true;
                end
            end
            % Gray
            eventsValid = false;
            while eventsValid == false
                cueEventFrames{trialI}(picsPerArray+1:end,:) = 0;
                bufferCheck = zeros(1,durInFrames+cueEventBuffer);
                eventLats = randi([blLengthVec(trialI)+1, durInFrames-cueEventDur], nrOfEventsVecGray(trialI), 1);
                for latI = 1:length(eventLats)
                    cueEventFrames{trialI}(randi([picsPerArray+1, 2*picsPerArray]),eventLats(latI):eventLats(latI)+cueEventDur) = 1;
                    bufferCheck(eventLats(latI):eventLats(latI)+cueEventDur+cueEventBuffer) = ...
                        bufferCheck(eventLats(latI):eventLats(latI)+cueEventDur+cueEventBuffer)+1;
                end
                if sum(bufferCheck > 1) == 0
                    eventsValid = true;
                end
            end            
        end
    end
end
    
    

function [actFlickDur, startTime] = presTrial_ssv4att_MRI(window, textureVec, ...
                                                          posVecs, flickVecs, blFrames, ...
                                                          centTileCol, centTileCoord, ...
                                                          fixVec, cueVec, ...
                                                          lptObject, lptAdress, lptCode)
% presents a flickering image texture; flickVec is a one-dimensional vector
% with opaqueness values from 0 to 1 for each frame; imageSize can be empty
% (no scaling of original image) or have two values for size in pixels ([x
% y]); texture is texture object; returns actual flicker stimulus duration
% in seconds obtained with tic & toc commands
    [winSize(1), winSize(2)] = Screen('WindowSize', window);
    winCenter = winSize ./ 2;

    tSize = Screen('TextSize', window);

    startTime = GetSecs();

    % Present Flicker Stimulus (Baseline)
    for frame = 1:blFrames
        Screen('DrawTextures', window, textureVec, [], posVecs(:,:,frame), 0, [], flickVecs(:,frame));
        Screen('FillRect', window, centTileCol, centTileCoord);
        Screen('DrawLines', window, fixVec(:,:,frame), 5, 0);
        Screen('DrawText', window, cueVec(1,frame), winCenter(1)-tSize*.55, winCenter(2)-tSize*.75 ,0);
        Screen('Flip', window);
    end

    % send marker via LPT Port
    io64(lptObject, lptAdress, lptCode);

    % Present Flicker Stimulus (Cue)
    for frame = blFrames+1:size(flickVecs,2)
        Screen('DrawTextures', window, textureVec, [], posVecs(:,:,frame), 0, [], flickVecs(:,frame));
        Screen('FillRect', window, centTileCol, centTileCoord);
        Screen('DrawLines', window, fixVec(:,:,frame), 5, 0);
        Screen('DrawText', window, cueVec(1,frame), winCenter(1)-tSize*.55, winCenter(2)-tSize*.75 ,0);
        Screen('Flip', window);
    end
    
    % set LPT Port to 0
    io64(lptObject, lptAdress, 0);

    % check time passed since before flicker stim
    actFlickDur = GetSecs() - startTime;
end


function ratingMat = rateMRInopic(window, buttonsLeft, buttonsRight, buttonsOK, ...
                                  startSelect, pauseDur, termByRater, maxDurInSec, ...
                                  labelVec, anchorsVec, textVec, textCol)
% PTB rating routine for Likert scales with mouse in put & slider response mode;
% The routine will first acquire rating # 1 for all stimuli, then rating 2, 
% and so on; order of stimuli is randomized at beginning of rating
% turn (but constant across ratings within turn);
% item texts, scale ticks & labels can be easily adapted (see "header") within this function;
% scales everything relative to screen (except for font size)
% window: Psychtoolbox window pointer
% buttonsLeft: indices of keyboard buttons to move cursor to left
% buttonsRight: indices of keyboard buttons to move cursor to right
% buttonsOK: indices of keyboard buttons to select current rating
% startSelect: starting level of cursor on Likert scale (integer);
% special values (strings) = 'lowest', 'center', and 'highest'
% pauseDur: fixed duration of fixation cross after each rating in seconds
% termByRater: true = rating ends upon okay-button press; false = fixed
% duration (maxDurInSec)
% labelVec: [# ratings]-element cell array; within cells = numbers for scake levels (integers)
% anchorsVec: [# ratings]-element cell array; within cells = strings for lowest and highest levels
% textVec: [# ratings]-element cell array; within cells = item text (string)
% output: ratingMat = [# ratings] x 2 matrix; rows = ratings; 1st col =
% rating; 2nd col = rating duration


    %% header with default parameters

    if length(unique([size(textVec,1), size(anchorsVec,1), size(labelVec,1)])) > 1
        error('Cell arrays in header section need same number of lines')
    end
    
    noSelCol = [255 0 0]; % color for cursor before selection is confirmed
    selMadeCol = [0 255 0]; % color for cursor after selection is confirmed

    %% matrices that we need
    ratingMat = NaN(length(labelVec), 3);
    
    % automatic computation of sizes relative to screen
    % global values / coordinates
    [xTotal, yTotal] = Screen('WindowSize', window);
    %xCenter = xTotal / 2; %yCenter = yTotal ./ 2;

    % x coordinates
    widthLeft = round(xTotal * .2, 0);
    widthRight = round(xTotal * .8, 0);
    widthTicks = cell(length(labelVec),1);
    for i = 1:length(labelVec)
        widthTicks{i} = widthLeft : (widthRight-widthLeft)/(length(labelVec{i})-1) : widthRight;
    end

    % y coordinates
    heightQuestion = round(yTotal * .3, 0);
    heightScale = round(yTotal * .5, 0);
    heightLabels = round(yTotal * .6, 0);
    heightAnchors = round(yTotal * .65, 0);
    
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
        if strcmp(startSelect, 'lowest')
            currSelect = selMin;
        elseif strcmp(startSelect, 'center')
            currSelect = ceil((selMax+selMin)/2);
        elseif strcmp(startSelect, 'highest')
            currSelect = selMax;
        else
            currSelect = startSelect;
        end

        % initiate cursor properties
        pointCol = noSelCol;

        startTime = GetSecs;
        keepRunning = true;

        while keepRunning == true
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
                if termByRater == true
                    keepRunning = false;
                end                
            end
            %when time is as good as up (within 20 ms), set cursor color to "selected"
            if (GetSecs()-startTime) > (maxDurInSec-.020)
                pointCol = selMadeCol;
                WaitSecs(.020);
                keepRunning = false;
            end

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
            Screen('Flip', window);
        end

        % logging rating duration
        ratDur = GetSecs() - startTime;

        % log rating
        ratingMat(RatI, 1) = labelVec{RatI}(currSelect);
        ratingMat(RatI, 2) = ratDur;

        % show selected rating for 250 ms, draw fixation cross and show it for the remainder of [pauseDur] seconds
        % if [pauseDur] is shorter than 250 ms, selection will be shown for this amount of time, no fixation cross
        if pauseDur > .250
            WaitSecs(.250);
            presFix(window, pauseDur-.250, textCol);
        else
            WaitSecs(pauseDur);
        end
    end
end