function particles = extractParticles(multiBernoulli,numSteps,parameters)

detectionThreshold = parameters.detectionThreshold;

isNotRedundant = [multiBernoulli.existence] >= detectionThreshold;
multiBernoulli = multiBernoulli(isNotRedundant);

numTracks = length(multiBernoulli);
particles = cell(numSteps,numTracks);
for target = 1:numTracks
    % find map estimate of start time
    particleStartTime = multiBernoulli(target).particleStartTime;
    uniqueStartTime = unique(particleStartTime);
    numUniqueStartTime = length(uniqueStartTime);
    probabilityMassStartTime = zeros(numUniqueStartTime,1);
    for startTime = 1:numUniqueStartTime
        probabilityMassStartTime(startTime) = sum(multiBernoulli(target).particleWeights(particleStartTime==uniqueStartTime(startTime)));
    end
    [~,indexMapStartTime] = max(probabilityMassStartTime);
    mapStartTime = uniqueStartTime(indexMapStartTime);

    % find map estimate of end time
    particleEndTime = multiBernoulli(target).particleEndTime;
    uniqueEndTime = unique(particleEndTime);
    numUniqueEndTime = length(uniqueEndTime);
    probabilityMassEndTime = zeros(numUniqueEndTime,1);
    for endTime = 1:numUniqueEndTime
        probabilityMassEndTime(endTime) = sum(multiBernoulli(target).particleWeights(particleEndTime==uniqueEndTime(endTime)));
    end
    [~,indexMapEndTime] = max(probabilityMassEndTime);
    mapEndTime = uniqueEndTime(indexMapEndTime);

    for step = mapStartTime:mapEndTime
        particles{step,target}.wp = multiBernoulli(target).weightHistory{step};
        particles{step,target}.xp = multiBernoulli(target).kinematicHistory{step};
        particles{step,target}.ep = multiBernoulli(target).extentHistory{step};
    end
end

end