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
        error('only "box" and "sin" are allowed for genFlickVec modulation input\n');
    end
    
    % check whether image size is empty or has two values (error)
    if ~isempty(imageSize) && length(imageSize) ~= 2
        error('imageSize parameter must be empty or contain both x and y values\n');
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
        fprintf(['it is recommended to use only even numbers for frames per ' ...
                 'cycle. If you are sure about this, continue with a key press\n']);
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