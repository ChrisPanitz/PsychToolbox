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