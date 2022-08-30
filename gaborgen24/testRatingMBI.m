function testRatingMBI()
%% SET TO ZERO FOR ACTUAL EXPERIMENT
Screen('Preference', 'SkipSyncTests', 1);


%% Header

% rating & instruction parameters
TRdur = 2;
pauseBtwRat = 1*TRdur; % pause duration between two ratings, in sec
ratingDur = 3*TRdur; % time out for ratings, in sec
instructDur = 3*TRdur; % fixed duration of instruction screens, in sec
buttonsLeft = KbName('LeftArrow'); % DEC key codes for moving rating cursor left
buttonsRight = KbName('RightArrow'); % DEC key codes for moving rating cursor right
buttonsOK = KbName('q');

% aesthetics
backgroundCol = 127.5; % bakcground color; 127.5 = mid gray
fontSize = 24; % well... font size

preRatMsg = ['Please rate for each picture how unpleasant and arousing you find it,\n' ...
             'and how likely you think it will be followed by a shock.\n' ...
             'You can move the cursor to the left or right\n' ...
             'by using the buttons in your hand. You have ' int2str(ratingDur) ' seconds\n'  ...
              'to move the cursor around before your rating is saved.'];
postRatMsg = 'Thank you for rating the pictures.';


%% clearing and initializing stuff
IOPort('CloseAll');
clc;

% Setting rng seed; "old" rand function (instead of rng) was used for
% reasons of compatibility with lab computer
rand('state',sum(100*clock));

AssertOpenGL; 





%% PTB code is in try loop
try
    %% get all info for screen & window
    Screen('CloseAll');
    screens = Screen('Screens');
    screenNumber = max(screens);
    w = Screen('OpenWindow', screenNumber, backgroundCol);
        
    % set alpha blending mode for manipulating picture transparency
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % Compute screen's PPI via window size (pix) and display size (mm)
    %ppi = Screen('WindowSize', w) / Screen('DisplaySize', screenNumber) * 25.4;
    % translate visual angle into pixels for picture scaling
    %imSizePix = pixFromAngle(imSizeAng, seatDistInch, ppi); % size of Gabors

    % set font size
    Screen('TextSize', w, fontSize);    

    %% prepare geometric shapes


    squarePixRed = cat(3, 255*ones(300,300), zeros(300,300), zeros(300,300));
    squarePixBlue = cat(3, zeros(300,300), zeros(300,300), 255*ones(300,300));
    
    % create texture object
    textureVec = NaN(2,1);
    % texture of red circle...
    textureVec(1) = Screen('MakeTexture', w, squarePixRed);
    % ... and of blue circle
    textureVec(2) = Screen('MakeTexture', w, squarePixBlue);

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%   HERE COMES THE ACTUAL EXPERIMENT   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HideCursor();
        
    presFix(w, TRdur);

    % instructions for ratings
    DrawFormattedText(w, preRatMsg, 'center', 'center');
    Screen('Flip', w);
    WaitSecs(instructDur);

    % run rating routine
    rateCSMRI(w, buttonsLeft, buttonsRight, buttonsOK, ...
              ratingDur, pauseBtwRat, textureVec);

    % instructions after ratings
    DrawFormattedText(w, postRatMsg, 'center', 'center');
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


% computes # of pixels required for picture based on the required visual
% angle, seating distance in inches, and the PPI of the screen; can take
% multiple and compute vector of values in one call
function nrPix = pixFromAngle(visAngle, distInch, ppi)
    nrPix = round(tan(visAngle/360*pi) * ppi * 2 * distInch, 0);
end


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
               %'How fearful did you feel when you saw this pattern?'; ...
               'How likely (in %) will this pattern be followed by a shock?'};

    anchorsVec = {'very pleasant', 'very unpleasant'; ...
                  'not arousing at all', 'very arousing'; ...
                  %'not fearful at all', 'very fearful'; ...
                  'never followed', 'always followed'};
    
    labelVec = {0:10; ...
                0:10; ...
                %0:10; ...
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