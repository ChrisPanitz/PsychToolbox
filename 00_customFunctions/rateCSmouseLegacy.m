function ratingMat = rateCSmouseLegacy(window, buttonsOK, pauseDur, continuousRatings, csTextureVec)
% PTB rating routine for a threat conditioning paradigm with habituation, 
% acquistion, and extinction; items have partricipants rate stimulus
% valence, arousal, their fear, and US expectancy; Likert scales with
% slider response mode;
% The routine will first acquire valence ratings for all CS, then arousal,
% fear, and expectancy; order of CS is randomized at beginning of rating
% trun (but constant within turn);
% item texts, scale ticks & labels can be easily adapted (see "header");
% scales everything relative to screen (except for font size)
% buttonsOK: indices of mouse buttons to select current rating
% pauseDur: fixed duration of fixation cross after each rating in seconds
% continuousRatings: false = Likert scale, only ticks can be selected; true
% = continuous ratings, i.e., participants can select values between ticks
% csTextureVec: vector with textures for all CS to be rated
% output: ratingMat = 1 x [nrRatings*3] vector with all ratings (valence_stim1, 
% valence_stim2 ... arousal_stim1 ... fear_stim1 ... expectancy_last_stimulus);
% after ratings, x coordinates of mouse cursor & rating durations are
% logged (same order as ratings)


    %% header with parameters
    textValence = 'How unpleasant was this pattern to you?';
    anchorsValence = {'very pleasant', 'very unpleasant'};
    textArousal = 'How arousing was this pattern to you?';
    anchorsArousal = {'not arousing at all', 'very arousing'};
    textFear = 'How fearful did you feel when you saw this pattern?';
    anchorsFear = {'not fearful at all', 'very fearful'};
    textExpectancy = 'How likely (in %) will this patern be followed by a shock?';
    anchorsExpectancy = {'never followed', 'always followed'};
    
    labelsUnipolar = 0:10;
    labelsBipolar = 0:10;
    labelsExpectancy = 0:10:100;

    textCol = [0 0 0]; % color for text & scale
    noSelCol = [255 0 0]; % color for cursor before selection is confirmed
    selMadeCol = [0 255 0]; % color for cursor after selection is confirmed

    %% matrices that we need
    orderVec = randperm(length(csTextureVec));
    valenceMat = NaN(1, length(csTextureVec));
    valenceMatX = NaN(1, length(csTextureVec));
    valenceMatDur = NaN(1, length(csTextureVec));
    arousalMat = NaN(1, length(csTextureVec));
    arousalMatX = NaN(1, length(csTextureVec));
    arousalMatDur = NaN(1, length(csTextureVec));
    fearMat = NaN(1, length(csTextureVec));
    fearMatX = NaN(1, length(csTextureVec));
    fearMatDur = NaN(1, length(csTextureVec));
    expectancyMat = NaN(1, length(csTextureVec));
    expectancyMatX = NaN(1, length(csTextureVec));
    expectancyMatDur = NaN(1, length(csTextureVec));

    % automatic computation of sizes relative to screen
    % global values / coordinates
    [xTotal, yTotal] = Screen('WindowSize', window);
    xCenter = xTotal / 2; %yCenter = yTotal ./ 2;

    % x coordinates
    widthLeft = round(xTotal * .2, 0);
    widthRight = round(xTotal * .8, 0);
    widthTicksUnipolar = widthLeft : (widthRight-widthLeft)/(length(labelsUnipolar)-1) : widthRight;
    widthTicksBipolar = widthLeft : (widthRight-widthLeft)/(length(labelsBipolar)-1) : widthRight;
    widthTicksExpectancy = widthLeft : (widthRight-widthLeft)/(length(labelsExpectancy)-1) : widthRight;

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
                xDist = abs(x - widthTicksBipolar);
                xPointer = min(widthTicksBipolar(xDist == min(xDist)));
            end
            % detect clicks; in case of click, log rating and end loop
            if sum(buttons(buttonsOK) > 0)
                relPos = (xPointer-widthLeft) / (widthRight-widthLeft);
                loggedRating = labelsBipolar(1) + relPos*(labelsBipolar(end)-labelsBipolar(1));
                pointCol = selMadeCol;
                selMade = true;
            elseif sum(buttons(buttonsOK)) == 0
                pointCol = noSelCol;
            end

            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textValence, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksBipolar, [], textCol);
            Screen('DrawLine', window, pointCol, ...
                               xPointer, heightScale-.05*yTotal, ...
                               xPointer, heightScale+.05*yTotal, 3);
            
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
        end
        
        ratDur = GetSecs() - startTime;

        % log rating
        valenceMat(orderVec(StimI)) = loggedRating;
        valenceMatX(orderVec(StimI)) = xPointer;
        valenceMatDur(orderVec(StimI)) = ratDur;
        
        % show the marked selected rating for a moment
        WaitSecs(.3);

        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end % valence rating


    %% Arousal Rating
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
                xDist = abs(x - widthTicksUnipolar);
                xPointer = min(widthTicksUnipolar(xDist == min(xDist)));
            end
            % detect clicks; in case of click, log rating and end loop
            if sum(buttons(buttonsOK) > 0)
                relPos = (xPointer-widthLeft) / (widthRight-widthLeft);
                loggedRating = labelsUnipolar(1) + relPos*(labelsUnipolar(end)-labelsUnipolar(1));
                pointCol = selMadeCol;
                selMade = true;
            elseif sum(buttons(buttonsOK)) == 0
                pointCol = noSelCol;
            end

            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textArousal, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksUnipolar, [], textCol);
            Screen('DrawLine', window, pointCol, ...
                               xPointer, heightScale-.05*yTotal, ...
                               xPointer, heightScale+.05*yTotal, 3);
            
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
        end
        
        ratDur = GetSecs() - startTime;

        % log rating
        arousalMat(orderVec(StimI)) = loggedRating;
        arousalMatX(orderVec(StimI)) = xPointer;
        arousalMatDur(orderVec(StimI)) = ratDur;
        
        % show the marked selected rating for a moment
        WaitSecs(.3);

        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end % arousal rating


    %% Fear Rating
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
                xDist = abs(x - widthTicksUnipolar);
                xPointer = min(widthTicksUnipolar(xDist == min(xDist)));
            end
            % detect clicks; in case of click, log rating and end loop
            if sum(buttons(buttonsOK) > 0)
                relPos = (xPointer-widthLeft) / (widthRight-widthLeft);
                loggedRating = labelsUnipolar(1) + relPos*(labelsUnipolar(end)-labelsUnipolar(1));
                pointCol = selMadeCol;
                selMade = true;
            elseif sum(buttons(buttonsOK)) == 0
                pointCol = noSelCol;
            end

            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textFear, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksUnipolar, [], textCol);
            Screen('DrawLine', window, pointCol, ...
                               xPointer, heightScale-.05*yTotal, ...
                               xPointer, heightScale+.05*yTotal, 3);
            
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
        end
        
        ratDur = GetSecs() - startTime;

        % log rating
        fearMat(orderVec(StimI)) = loggedRating;
        fearMatX(orderVec(StimI)) = xPointer;
        fearMatDur(orderVec(StimI)) = ratDur;
        
        % show the marked selected rating for a moment
        WaitSecs(.3);

        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end % fear rating


    %% Expectancy Rating
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
                xDist = abs(x - widthTicksExpectancy);
                xPointer = min(widthTicksExpectancy(xDist == min(xDist)));
            end
            % detect clicks; in case of click, log rating and end loop
            if sum(buttons(buttonsOK) > 0)
                relPos = (xPointer-widthLeft) / (widthRight-widthLeft);
                loggedRating = labelsExpectancy(1) + relPos*(labelsExpectancy(end)-labelsExpectancy(1));
                pointCol = selMadeCol;
                selMade = true;
            elseif sum(buttons(buttonsOK)) == 0
                pointCol = noSelCol;
            end

            % Draw Stimulus
            Screen('DrawTexture', window, csTextureVec(orderVec(StimI)),[], destRectCs);
            
            % Draw Question
            DrawFormattedText(window, textExpectancy, 'center', heightQuestion, textCol);
            
            % Draw Scale
            Screen('DrawLine', window, textCol, widthLeft, heightScale, widthRight, heightScale);
            Screen('DrawLines', window, coordsTicksExpectancy, [], textCol);
            Screen('DrawLine', window, pointCol, ...
                               xPointer, heightScale-.05*yTotal, ...
                               xPointer, heightScale+.05*yTotal, 3);
            
            % Draw Labels
            for lab = 1:length(labelsExpectancy)
                DrawFormattedText(window, int2str(labelsExpectancy(lab)), ...
                                  'center', heightLabels, textCol, ...
                                  [], [], [], [], [], ...
                                  [widthTicksExpectancy(lab), heightLabels, widthTicksExpectancy(lab), heightLabels]);
            end
            
            % Draw Anchors
            DrawFormattedText(window, anchorsExpectancy{1}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksExpectancy(1), heightAnchors, widthTicksExpectancy(1), heightAnchors]);
            DrawFormattedText(window, anchorsArousal{2}, ...
                              'center', heightAnchors, textCol, ...
                              [], [], [], [], [], ...
                              [widthTicksExpectancy(end), heightAnchors, widthTicksExpectancy(end), heightAnchors]);
            
            % Flip it!
            Screen('Flip', window)
        end
        
        ratDur = GetSecs() - startTime;

        % log rating
        expectancyMat(orderVec(StimI)) = loggedRating;
        expectancyMatX (orderVec(StimI)) = xPointer;
        expectancyMatDur(orderVec(StimI)) = ratDur;
        
        % show the marked selected rating for a moment
        WaitSecs(.3);

        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur);
    end % arousal rating


    % concatenate ratings for the different items; result is a 1 x nrRatings vector
    ratingMat = cat(2, valenceMat, arousalMat, fearMat, expectancyMat, ...
                       valenceMatX, arousalMatX, fearMatX, expectancyMatX, ...
                       valenceMatDur, arousalMatDur, fearMatDur, expectancyMatDur);
end