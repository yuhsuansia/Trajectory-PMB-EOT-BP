function [estimatedTracks,estimatedExtents] = allTrajectoryFormation(multiBernoulli,step,parameters)

detectionThreshold = parameters.detectionThreshold;

isNotRedundant = [multiBernoulli.existence] >= detectionThreshold;
multiBernoulli = multiBernoulli(isNotRedundant);

numTracks = length(multiBernoulli);
estimatedTracks = nan(4,step,numTracks);
estimatedExtents = nan(2,2,step,numTracks);

for target = 1:numTracks
    particleStartTime = multiBernoulli(target).particleStartTime;
    uniqueStartTime = unique(particleStartTime);
    numUniqueStartTime = length(uniqueStartTime);
    probabilityMassStartTime = zeros(numUniqueStartTime,1);
    for startTime = 1:numUniqueStartTime
        probabilityMassStartTime(startTime) = sum(multiBernoulli(target).particleWeights(particleStartTime==uniqueStartTime(startTime)));
    end
    [~,indexMapStartTime] = max(probabilityMassStartTime);
    mapStartTime = uniqueStartTime(indexMapStartTime);

    particleEndTime = multiBernoulli(target).particleEndTime;
    uniqueEndTime = unique(particleEndTime);
    numUniqueEndTime = length(uniqueEndTime);
    probabilityMassEndTime = zeros(numUniqueEndTime,1);
    for endTime = 1:numUniqueEndTime
        probabilityMassEndTime(endTime) = sum(multiBernoulli(target).particleWeights(particleEndTime==uniqueEndTime(endTime)));
    end
    [~,indexMapEndTime] = max(probabilityMassEndTime);
    mapEndTime = uniqueEndTime(indexMapEndTime);

    for time = mapStartTime:mapEndTime
        wp = multiBernoulli(target).weightHistory{time};
        xp = multiBernoulli(target).kinematicHistory{time};
        ep = multiBernoulli(target).extentHistory{time};

        estimatedTracks(:,time,target) = sum(xp.*wp',2);
        weightReshape = ones(1,1,length(wp));
        weightReshape(1,1,:) = wp;
        estimatedExtents(:,:,time,target) = sum(ep.*weightReshape,3);
    end
end

end