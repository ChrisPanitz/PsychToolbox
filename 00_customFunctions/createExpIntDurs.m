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