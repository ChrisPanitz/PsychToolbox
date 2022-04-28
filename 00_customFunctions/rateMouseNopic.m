function ratingMat = rateMouseNopic(window, buttonsOK, pauseDur, continuousRatings, termByRater, maxDurInSec, labelVec, anchorsVec, textVec, textCol)
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

        % initiate cursor properties
        pointCol = noSelCol;
        %selMade = false;

        startTime = GetSecs;
        loggedRating = -999;
        keepRunning = true;

        while keepRunning == true
            % get mouse coordinates
            if loggedRating == -999
                [x,~,buttons] = GetMouse(window);
            end
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
                if termByRater == true
                    keepRunning = false;
                end
            %elseif sum(buttons(buttonsOK)) == 0
            %    pointCol = noSelCol;
            end

            % log rating duration and check for time out
            ratDur = GetSecs() - startTime;
            if ratDur > maxDurInSec
                %pointCol = selMadeCol;
                keepRunning = false;
            end

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
            Screen('Flip', window);
        end

        % log rating
        ratingMat(RatI, 1) = loggedRating;
        ratingMat(RatI, 2) = xPointer;
        ratingMat(RatI, 3) = ratDur;

        % show the marked selected rating for a moment
        %WaitSecs(.3);

        % draw fixation cross and wait for [pauseDur] seconds
        presFix(window, pauseDur, textCol);
    end

    % turn 3D ratingMat into 1 x [number all values] vector (i.e., a line)
    ratingMat = ratingMat(:)';
end