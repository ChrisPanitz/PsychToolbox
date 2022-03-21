function shuffledMat = pseudoShuffleWithBooster(condVec, trPerCS, nrBoostersMassed, nrBoostersRandom, maxTrBtwBooster, maxSeqCS, percPairing, pairingProbLocal) 
    
    clusteredCS = 1; % any value > 0
    while clusteredCS > 0
        % trials for the booster "phase"
        boostMassInd = NaN(nrBoostersMassed,1);
        if sum(nrBoostersMassed + nrBoostersRandom) > 0
            boostMassInd(1) = randi(maxTrBtwBooster+1);
        end
        for i = 2:length(boostMassInd)
            boostMassInd(i) = boostMassInd(i-1) + randi(maxTrBtwBooster+1);
        end
        boostMat = zeros(length(boostMassInd),2);
        boostMat(boostMassInd,1) = 1;
        boostMat(boostMat(:,1) == 0,1) = randi([2,max(condVec)], sum(boostMat(:,1) == 0), 1);

        remainTrials = repmat(trPerCS, 1, length(condVec)) - histcounts(boostMat(:,1), 1:condVec(end)+1);

        afterBoostMat = zeros(sum(remainTrials),2);
        for condI = condVec
            afterBoostMat(sum(remainTrials(1:condI-1))+1 : sum(remainTrials(1:condI)), 1) = condI;
        end
    
        % shuffle after-massed-booster trial order and concatenate all trials
        afterBoostMat = afterBoostMat(randperm(sum(remainTrials)),:); % shuffle trial order
        shuffledMat = cat(1, boostMat, afterBoostMat);

         % 0 if the same CS occurs twice in sequence
        diffTrials = diff(shuffledMat(:,1));
        % adds neighboring values in groups of [maxSeqCS]
        % => [maxSeqCS] times no change is too much
        % "abs" is to avoid sums of zero from other values than only zeros
        seqVec = filter(ones(maxSeqCS,1), 1, abs(diffTrials)); 
        % first [maxSeqCS-1] values are not based on the full set of [maxSeqCS] elements
        seqVec = seqVec(maxSeqCS:end); 
        clusteredCS = sum(seqVec == 0); % if there are zeros ==> keep shuffling
    end

    % US [yes/no] column
    if pairingProbLocal == true
        nrUSafterBooster = floor((trPerCS-nrBoostersMassed-nrBoostersRandom) * percPairing);
    elseif pairingProbLocal == false
        nrUSafterBooster = floor(trPerCS*percPairing - nrBoostersMassed - nrBoostersRandom);
    end
    usVecAftBoost = zeros(remainTrials(1)-nrBoostersRandom,1);
    usVecAftBoost(1:nrUSafterBooster) = 1;
    usVecAftBoost = usVecAftBoost(randperm(length(usVecAftBoost)));
    shuffledMat(shuffledMat(:,1) == 1, 2) = cat(1, ones(nrBoostersMassed+nrBoostersRandom,1), usVecAftBoost);
end