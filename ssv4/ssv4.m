function ssv4(subNo)
%% TAKE OUT FOR ACTUAL EXPERIMENT OR SET TO ZERO
Screen('Preference', 'SkipSyncTests', 1);


%% Header
% for logfiles
logFileFolder = [pwd '/logfiles/'];
logFilePrefix = 'ssv4'; % change for new study

% trial parameters 
trPerPic = 10; % # of picture presentations per pic & condition
breakAftTr = 100; % breaks after how many trials?
minMaxItiSec = [.8 3]; % min & max ITI in sec, taken from exponential distribution

% flicker parameters
% frPerCyc for 120Hz screen: 20 = 6 Hz, 14 = 8.57 Hz, 8 = 15 Hz
frPerCyc = [20 14 8]; % duration of one cycle in frames (for different frequencies)
flickMods = {'box'; 'sin'}; % how flicker stimulus is modulated. 'box' or 'sin' or both possible
flickDurFrames = 280; % duration of whole stimulus in frames

% image parameters
imSizeAng = [];%[4 4]; % visual angle, horizontal & vertical
% distance participant <-> screen
seatDistInch = 80/2.54; 
% folder with images
imFolder = [pwd '/imageStimuli/']; % path for image files
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

% Setting some form of seed
rng(sum(100*clock));

%Initialize sound driver & push for low latencies (1)
%InitializePsychSound(1) 

AssertOpenGL; 

%% Prepare logging
% open serial port
%s3 = serialport('com3', 9600);
%fopen(s3);

% file to write trial data in + variable names for header line
logFileName = strcat(logFileFolder, logFilePrefix, '_', num2str(subNo), '.dat');
fID = fopen(logFileName, 'w');
fprintf(fID,'subNo,trial,freq,mod,pic,actFlickDur,fixDur\n');
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
        fixDurThisTime = randExpoInt(.1, minMaxItiSec); 
        % present trial while logging actual flicker duration
        presFix(w, fixDurThisTime);
        actFlickDur = presFlick(w, condVec(trial,:), flickerVecs, imSizePix, textureVec);
        
        % write parameters to data file
        trialOutVec = [subNo trial condVec(trial,:) actFlickDur fixDurThisTime];
        dlmwrite(logFileName, trialOutVec, '-append');
        
        % Short break for participants (after [pauseAftTr] trials)
        if trial/breakAftTr == round(trial/breakAftTr, 0)
           DrawFormattedText(w, breakMsg, 'center', 'center');
           Screen('Flip', w);
           waitForClick;
           presFix(w, 5);
        end % break conditional
    end % trial loop
        
    % final screen, disappears automatically after 10 sec
    DrawFormattedText(w, goodbyeMsg, 'center', 'center');
    Screen('Flip', w);
    presFix(w, 10);

    % tidy up
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
