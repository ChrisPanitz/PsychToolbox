function paramCheck(framesPerCycle, framesTotal, flickMods, breakAfterTrials, imFileList, imageSize, window)
    % check whether framesTotal is a multiple of all framesPerCycle
    probFreqs = NaN(1,length(framesPerCycle));
    for i = length(framesPerCycle):-1:1
        if framesPerCycle(i)/framesTotal ~= round(framesPerCycle(i)/framesTotal,0)
            probFreqs(i) = framesPerCycle(i);
        else
            probFreqs(i) = [];
        end
        if ~isempty(probFreqs)
            error(["total number of frames must be an integer multiple of " ...
                   "each value of frames per cycle. The following values " ...
                   "do not fit into the total number of frames:\n" ...
                   num2str(probFreqs) "\n"]);
        end
    end
    
    % check whether flickMods has valid arguments
    if sum(ismember(flickMods, {'box','sin'})) < length(flickMods)
        error("only 'box' and 'sin' are allowed for genFlickVec modulation input\n");
    end
    
    % check whether image size is empty or has two values
    if ~isempty(imageSize) && length(imageSize) ~= 2
        error("imageSize parameter must be empty or contain both x and y values\n");
    end
    
    % check whether breakAfterTrials leads to equal-sized blocks of trials
    if length(imFileList)/breakAfterTrials ~= round(length(imFileList)/breakAfterTrials, 0)
        fprintf(["the total number of trials is not an integer multiple of " ...
                 "breakAfterTrials. This means that the last block will be" ...
                 "shorter than the others. If you're sure about this, continue" ...
                 "with a key press\n"]);
        KbStrokeWait;
    end
    
    % check whether number of frames for flicker stimulus is even
    if sum(framesPerCycle./2 == round(framesPerCycle./2,0)) < length(framesPerCycle)
        fprintf(["it is recommended to use only even numbers for frames per " ...
                 "cycle. If you're sure about this, continue with a key press\n"]);
        KbStrokeWait;
    end

    % check whether there are overlaps between frequencies and harmonics
    fRate = Screen('NominalFrameRate', window);
        frequencies = fRate./framesPerCycle;
        critFreqInd = sum(ismember(frequencies, frequencies.*2));
        if sum(critFreqInd)
            fprintf(["it is recommend not to use frequencies that are the " ...
                     "a harmonic of another used frequencies (i.e., f1 * x). " ...
                     "Critical frequencies: " num2str(frequencies(critFreqInd)) ...
                     " Hz (i.e., " num2str(frequencies(framesPerCycle)) " frames " ...
                     "per cycle). If you're sure about this, continue with a key press\n"]);
            KbStrokeWait;
        end