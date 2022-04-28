function [positions] = createCircPositions(numPos, nrSamples, minRho, maxRho, minDistAllowed, maxPixDelta)
    positions = NaN(numPos,2,nrSamples);
    
    for posI = 1:numPos
        minDist = 0;
        while minDist < minDistAllowed
            % randomly draw position within the min/max radius restriction
            [positions(posI,1,1),positions(posI,2,1)] = pol2cart(2*pi*rand, minRho+(maxRho-minRho)*rand);
            % compute distance to closest existion position (Pythagoras)
            minDist = min(sqrt(sum((positions(posI,:,1) - positions(1:posI-1,:,1)).^2, 2)));
        end % minDist loop
    end % posI loop

    for sampI = 2:nrSamples
        allRhoDist = 0;
        % loop until all distances are at least the allowed minimum
        while min(allRhoDist) < minDistAllowed
            % add x and y shift to previous positions
            positions(:,:,sampI) = positions(:,:,sampI-1) + maxPixDelta*randi([-1,1],numPos,2);
            % compute distances from center 
            [~, rho] = cart2pol(positions(:,1,sampI), positions(:,2,sampI));
            % positions inside fixation are pushed outside
            positions(rho < minRho,:,sampI) = positions(rho < minRho,:,sampI) + ...
                positions(rho < minRho,:,sampI)./abs(positions(rho < minRho,:,sampI))*maxPixDelta;
            % positions outside stimulus field are pushed inside
            positions(rho > maxRho,:,sampI) = positions(rho > maxRho,:,sampI) - ...
                positions(rho > maxRho,:,sampI)./abs(positions(rho > maxRho,:,sampI))*maxPixDelta;
            % compute all distances between positions (in pixels on x and y) 
            allXDist = positions(:,1,sampI) - positions(:,1,sampI)';
            allXDist = allXDist(logical(1-eye(numPos)));
            allYDist = positions(:,2,sampI) - positions(:,2,sampI)';
            allYDist = allYDist(logical(1-eye(numPos)));
            % compute direct distance
            [~,allRhoDist] = cart2pol(allXDist,allYDist);
        end
    end

end


