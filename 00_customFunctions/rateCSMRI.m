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