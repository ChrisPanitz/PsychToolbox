function ssv4(subNo)
%% SET TO ZERO FOR ACTUAL EXPERIMENT
Screen('Preference', 'SkipSyncTests', 0);


%% Header
% for logfiles
logFileFolder = '/home/andreaskeil/Desktop/As_Exps/ssv4/logfiles/';
logFilePrefix = 'ssv4'; % change for new study

% trial parameters 
trPerPic = 5; % # of picture presentations per pic & condition
breakAftTr = 100; % breaks after how many trials?
minMaxItiSec = [.8 3]; % min & max ITI in sec, taken from exponential distribution

% flicker parameters
% frPerCyc for 120Hz screen: 20 = 6 Hz, 14 = 8.57 Hz, 8 = 15 Hz
frPerCyc = [20 14 8]; % duration of one cycle in frames (for different frequencies)
flickMods = {'box'; 'sin'}; % how flicker stimulus is modulated. 'box' or 'sin' or both possible
flickDurFrames = 280; % duration of whole stimulus in frames

% image parameters
% no rescaling, original images scaled for 4x4 ° visual angle at 125 cm
% seating distance and 68.8 PPI
imSizeAng = []; 
% distance participant <-> screen
seatDistInch = 125/2.54; 
% folder with images
imFolder = '/home/andreaskeil/Desktop/As_Exps/ssv4/imageStimuli/'; % path for image files
% list of filenames for images
imFileList = {'pic01.jpg', 'pic02.jpg', 'pic03.jpg', 'pic04.jpg', ...
              'pic05.jpg', 'pic06.jpg', 'pic07.jpg', 'pic08.jpg', ...
              'pic09.jpg', 'pic10.jpg', 'pic11.jpg', 'pic12.jpg', ...
              'pic13.jpg', 'pic14.jpg', 'pic15.jpg', 'pic16.jpg', ...
              'pic17.jpg', 'pic18.jpg', 'pic19.jpg', 'pic20.jpg'};

% strings of instructions to participants          
welcomeMsg_expStart = ['Welcome and thank you for participating in our\n' ...
                       'experiment. The experimenter will start the experiment soon.'];
welcomeMsg_partStart = ['Everything has been set up now. In the upcoming\n' ...
                        'task you will be presented flickering pictures\n' ...
                        'in the center of the screen. Your task is to\n' ...
                        'pay close attention to these pictures. The task\n' ...
                        'will take about 35 to 40 minutes with several\n' ...
                        'short breaks in between.\n\n\n' ...
                        'If you have questions, please let the\n' ...
                        'experimenter know now. When you are ready\n' ...
                        'you can start the task by pressing any mouse button.'];                   
breakMsg = ['You can take a short break now.\n\n' ...
            'When you are ready, please continue the experiment\n' ...
            'by pressing any mouse button.'];
goodbyeMsg = ['You have completed the task. Thank you for your participation!\n\n' ...
              'Please remain seated. The experimenter will be with you in a moment.'];

          

%% clearing and initializing stuff
IOPort('CloseAll');
clc;

% Setting rng seed; "old" rand function (instead of rng) was used for
% reasons of compatibility with lab computer
rand('state',sum(100*clock));

%Initialize sound driver & push for low latencies (1)
%InitializePsychSound(1) 

AssertOpenGL; 

%% Prepare logging
% open serial port
[s3, ~] = IOPort('OpenSerialPort', '/dev/ttyS0', 'Baudrate = 9600');
%fopen(s3);

% file to write trial data in + variable names for header line
logFileName = strcat(logFileFolder, logFilePrefix, '_', num2str(subNo), '.dat');
fID = fopen(logFileName, 'w');
fprintf(fID,'subNo,trial,cond,pic,actFlickDur,fixDur\n');
fclose(fID);

%% actual experiment is in try loop
try
    %% get all info for screen & window
    screens = Screen('Screens');
    screenNumber = max(screens);
    w = Screen('OpenWindow', screenNumber, [122 122 122]);
        
    % set alpha blending mode for manipulating picture transparency
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % Compute screen's PPI via window size (pix) and display size (mm)
    ppi = Screen('WindowSize', w) / Screen('DisplaySize', screenNumber) * 25.4;
    % translate visual angle into pixels for picture scaling
    imSizePix = pixFromAngle(imSizeAng, seatDistInch, ppi);

    %% check header parameters for plausibility
    paramCheck(frPerCyc, flickDurFrames, flickMods, breakAftTr, ...
               length(frPerCyc)*length(flickMods)*length(imFileList)*trPerPic, ...
               imSizeAng, w);
           
    %% create vectors & matrices for stimuli & trials
    % Dynamically creates matrix containing indices to 'look up' frequency
    % (in frPerCyc), flicker modulation (in flickMods), and picture (in
    % textureVec)... for each trial
    % #rows = total trial number
    % 1st col = frequency, 2nd col = modulation, 3rd col = picture
    conds = [sort(repmat(1:length(frPerCyc), 1, length(flickMods)*length(imFileList)))', ...
             repmat(sort(repmat(1:length(flickMods), 1, length(imFileList))), 1, length(frPerCyc))',  ...
             repmat(1:length(imFileList), 1, length(frPerCyc)*length(flickMods))'];
    condVec = repmat(conds, trPerPic, 1);
    condVec = condVec(randperm(size(condVec,1)),:);
   
    % flicker sequences (on/off frames) for different tagging frequencies
    flickerVecs = NaN(length(frPerCyc), length(flickMods), flickDurFrames);
    for freq = 1:length(frPerCyc)
        for modul = 1:length(flickMods)
            flickerVecs(freq,modul,:) = genFlickVec(frPerCyc(freq), flickDurFrames, flickMods{modul});
        end % modul loop
    end % freq loop

    %% load an prepare images    
    imageMat = cell(length(imFileList),1);
    textureVec = NaN(length(imFileList),1);
    for im = 1:length(imFileList)
        imageMat{im} = imread([imFolder imFileList{im}]);
        textureVec(im) = Screen('MakeTexture', w, imageMat{im});
    end % filling imageMat and textureVec (im loop)

    %%    BEGIN  MAIN EXPERIMENT %%%%%%%%
    HideCursor
    
    % welcome screen, initiated by experimenter terminated by participant
    DrawFormattedText(w, welcomeMsg_expStart, 'center', 'center');
    Screen('Flip', w);
    KbStrokeWait;
    DrawFormattedText(w, welcomeMsg_partStart, 'center', 'center');
    Screen('Flip', w);
    waitForClick;
    presFix(w, 5);
        
    %% trials start
    for trial = 1 : size(condVec,1)
        % set fix cross duration (sec) from exponential distribution with mu = .1
        fixDurThisTime = randExpoInt(minMaxItiSec); 
        
        % present trial while logging actual flicker duration
        presFix(w, fixDurThisTime);
        actFlickDur = presFlick(w, flickerVecs(condVec(trial,1),condVec(trial,2),:), ...
                                imSizePix, textureVec(condVec(trial,3)), s3);
        
        % write parameters to data file
        combCond = condVec(trial,1)*10 + condVec(trial,2);
        trialOutVec = [subNo trial combCond condVec(trial,3) actFlickDur fixDurThisTime];
        dlmwrite(logFileName, trialOutVec, '-append');
        
        % Short break for participants (after [pauseAftTr] trials but not after last one)
        if trial/breakAftTr == round(trial/breakAftTr, 0) && trial ~= size(condVec, 1)
           DrawFormattedText(w, breakMsg, 'center', 'center');
           Screen('Flip', w);
           waitForClick;
           presFix(w, 5);
        end % break conditional
    end % trial loop
        
    % final screen, disappears automatically after 10 sec
    DrawFormattedText(w, goodbyeMsg, 'center', 'center');
    Screen('Flip', w);
    KbStrokeWait;

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

% checks the validity of some some parameters and gives error messages or
% warnings. If parameters would result in error, terminate script and give
% explanation. If choice of parameters is problematic but does not result
% in error, give warning and allow user to continue with button press.
function paramCheck(framesPerCycle, framesTotal, flickMods, breakAfterTrials, totalTrials, imageSize, window)
    % check whether framesTotal is a multiple of all framesPerCycle (error)
    probFreqs = NaN(1,length(framesPerCycle));
    for i = length(framesPerCycle):-1:1
        if framesTotal/framesPerCycle(i) ~= round(framesTotal/framesPerCycle(i),0)
            probFreqs(i) = framesPerCycle(i);
        else
            probFreqs(i) = [];
        end
    end
    if ~isempty(probFreqs)
        error(['Total number of frames must be an integer multiple of ' ...
               'each value of frames per cycle. The following values ' ...
               'do not fit into the total number of frames: ' ...
               num2str(probFreqs)]);
    end
    
    % check whether flickMods has valid arguments (error)
    if sum(ismember(flickMods, {'box','sin'})) < length(flickMods)
        error('Only "box" and "sin" are allowed for genFlickVec modulation input.');
    end
    
    % check whether image size is empty or has two values (error)
    if ~isempty(imageSize) && length(imageSize) ~= 2
        error('imageSize parameter must be empty or contain both x and y values.');
    end
    
    % check whether breakAfterTrials leads to equal-sized blocks of trials
    % (warning)
    if totalTrials/breakAfterTrials ~= round(totalTrials/breakAfterTrials, 0)
        fprintf(['the total number of trials is not an integer multiple of ' ...
                 'breakAfterTrials. This means that the last block will be' ...
                 'shorter than the others. If you are sure about this, continue' ...
                 'with a key press\n']);
        KbStrokeWait;
    end
    
    % check whether number of frames for flicker stimulus is even (warning)
    if sum(framesPerCycle./2 == round(framesPerCycle./2,0)) < length(framesPerCycle)
        fprintf(['It is recommended to use only even numbers for frames per ' ...
                 'cycle. If you are sure about this, continue with a key press.']);
        KbStrokeWait;
    end

    % check whether there are overlaps between frequencies and harmonics
    % (warning)
    fRate = Screen('NominalFrameRate', window);
    frequencies = fRate./framesPerCycle;
    harmonics = repmat(frequencies,1,9) .* sort(repmat(2:10,1,length(framesPerCycle)));
    critFreqInd = sum(ismember(frequencies, harmonics));
    if sum(critFreqInd)
        fprintf(['It is recommend not to use frequencies that have ' ...
                 'a harmonic that is another used frequency (i.e., f1 * x). ' ...
                 'Critical frequencies: ' num2str(frequencies(critFreqInd)) ...
                 ' Hz (i.e., ' num2str(framesPerCycle(critFreqInd)) ' frames ' ...
                 'per cycle). If you are sure about this, continue with a key press.']);
        KbStrokeWait;
    end
end % function


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
function actFlickDur = presFlick(window, flickVec, imageSize, texture, portOut)
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
    tic;
    
    % Present Flicker Stimulus
    for frame = 1:size(flickVec, 3)
        Screen('DrawTexture', window, texture, [], imageCoordinates, 0, [], flickVec(frame));
        Screen('Flip', window);
    end
    
    % check time passed since before flicker stim
    actFlickDur = toc;
end
