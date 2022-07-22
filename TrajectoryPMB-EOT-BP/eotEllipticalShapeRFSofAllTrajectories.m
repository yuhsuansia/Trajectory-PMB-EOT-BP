function [ estimatedTracks, estimatedExtents, estimatedParticles, trajectoryMetric, runTime] = eotEllipticalShapeRFSofAllTrajectories( measurementsCell,targetTracks,targetExtents,parameters )

numParticles = parameters.numParticles;
numNewBornParticles = parameters.numNewBornParticles;

meanClutter = parameters.meanClutter;
meanMeasurements = parameters.meanMeasurements;
scanTime = parameters.scanTime;
meanBirths = parameters.meanBirths;
survivalProbability = parameters.survivalProbability;

priorVelocityCovariance = parameters.priorVelocityCovariance;
surveillanceRegion = parameters.surveillanceRegion;
areaSize =  (surveillanceRegion(2,1)-surveillanceRegion(1,1)) * (surveillanceRegion(2,2)-surveillanceRegion(1,2));
measurementsCovariance = parameters.measurementVariance * eye(2);

priorExtent1 = parameters.priorExtent1;
priorExtent2 = parameters.priorExtent2;
meanExtentPrior = priorExtent1/(priorExtent2-3);
totalCovariance = meanExtentPrior+measurementsCovariance;

thresholdPruning = parameters.thresholdPruning;
thresholdTruncate = parameters.thresholdTruncate;
numOuterIterations = parameters.numOuterIterations;

[numSteps, ~] = size(measurementsCell);
constantFactor = areaSize*(meanMeasurements/meanClutter);
probabilityEffectiveDetect = 1-exp(-meanMeasurements);
uniformWeight = log(1/areaSize);

multiBernoulli = repmat(struct('existence',[],'particleWeights',[],'particleStartTime',[],'particleEndTime',[],'particlesKinematic',[],'particlesExtent',[],'weightHistory',[],'kinematicHistory',[],'extentHistory',[]),[0,1]);
PoissonPointProcess.meanNumOfUndetectedObject = 0;

trajectoryMetric.total = zeros(numSteps,1);
trajectoryMetric.loc = zeros(numSteps,1);
trajectoryMetric.miss = zeros(numSteps,1);
trajectoryMetric.false = zeros(numSteps,1);
trajectoryMetric.switch = zeros(numSteps,1);

fprintf('Time step: ')

elapsedTime = zeros(numSteps,1);
for step = 1:numSteps

    tic

    fprintf('%d ',step)

    % load current measurements
    measurements = measurementsCell{step};
    numMeasurements = size(measurements,2);

    % perform prediction step
    multiBernoulli = performBernoulliAllTrajectoryPrediction(multiBernoulli,scanTime,step,parameters);
    numLegacy = length(multiBernoulli);
    indexAliveParticle = cell(numLegacy,1);
    for target = 1:numLegacy
        indexTmp = find(multiBernoulli(target).particleEndTime==step);
        aliveProbability = sum(multiBernoulli(target).particleWeights(indexTmp));
        multiBernoulli(target).existence = multiBernoulli(target).existence*(1-aliveProbability*probabilityEffectiveDetect)/(1-multiBernoulli(target).existence*aliveProbability*probabilityEffectiveDetect);
        multiBernoulli(target).particleWeights(indexTmp) = multiBernoulli(target).particleWeights(indexTmp)*(1-probabilityEffectiveDetect);
        multiBernoulli(target).particleWeights = multiBernoulli(target).particleWeights/sum(multiBernoulli(target).particleWeights);
        indexAliveParticle{target} = indexTmp;
    end
    PoissonPointProcess.meanNumOfUndetectedObject = survivalProbability*PoissonPointProcess.meanNumOfUndetectedObject+meanBirths;

    % get indexes of promising new objects
    [newIndexes,measurements,likelihood1] = getPromisingNewMeasurementAllTrajectory(multiBernoulli,measurements,step,indexAliveParticle,parameters);
    numNew = size(newIndexes,1);
    newMultiBernoulli = repmat(struct('existence',[],'particleWeights',[],'particleStartTime',[],'particleEndTime',[],'particlesKinematic',[],'particlesExtent',[],'weightHistory',[],'kinematicHistory',[],'extentHistory',[]),[numNew,1]);
    for target = 1:numNew
        newMultiBernoulli(target).particleStartTime = step*ones(numNewBornParticles,1);
        newMultiBernoulli(target).particleEndTime = newMultiBernoulli(target).particleStartTime;
        proposalMean = measurements(:,newIndexes(target));
        proposalCovariance = 2 * totalCovariance; % stretch covariance matrix to make proposal distribution heavier-tailed then target distribution
        newMultiBernoulli(target).particlesKinematic = zeros(4,numSteps,numNewBornParticles);
        newMultiBernoulli(target).particlesKinematic(1:2,step,:) = proposalMean + sqrtm(proposalCovariance) * randn(2,numNewBornParticles);
        newMultiBernoulli(target).particlesKinematic(3:4,step,:) = mvnrnd([0;0],priorVelocityCovariance,numNewBornParticles)';
        newMultiBernoulli(target).particlesExtent = zeros(2,2,numSteps,numNewBornParticles);
        newMultiBernoulli(target).particlesExtent(:,:,step,:) = iwishrndFastVector(priorExtent1,priorExtent2,numNewBornParticles);
        newMultiBernoulli(target).particleWeights = exp(normalizeLogWeights(uniformWeight-log_mvnpdf(reshape(newMultiBernoulli(target).particlesKinematic(1:2,step,:),2,numNewBornParticles)', proposalMean', proposalCovariance')));
        aliveProbability = PoissonPointProcess.meanNumOfUndetectedObject*(1-probabilityEffectiveDetect);
        newMultiBernoulli(target).existence = aliveProbability/(aliveProbability+1);
    end

    % initialize belief propagation (BP) message passing
    multiBernoulli = [multiBernoulli;newMultiBernoulli];
    currentExistenceExtrinsic = repmat([multiBernoulli.existence]',[1,numMeasurements]);
    numParticlesBernoulli = arrayfun(@(x) length(x.particleWeights),multiBernoulli);
    weightsExtrinsic = cell(numLegacy,1);
    for target = 1:numLegacy
        weightsExtrinsic{target} = zeros(numParticlesBernoulli(target),numMeasurements);
    end
    weightsExtrinsicNew = cell(numNew,1);
    for target = 1:numNew
        weightsExtrinsicNew{target} = zeros(numParticlesBernoulli(numLegacy+target),numMeasurements);
    end
    likelihoodNew1 = weightsExtrinsicNew;

    for outer = 1:numOuterIterations
        % perform one BP message passing iteration for each measurement
        outputDA = cell(numMeasurements,1);
        targetIndexes = cell(numMeasurements,1);
        for measurement = numMeasurements:-1:1
            inputDA = ones(2,numLegacy);
            for target = 1:numLegacy
                if(outer == 1)
                    inputDA(2,target) = currentExistenceExtrinsic(target,measurement) * multiBernoulli(target).particleWeights'*likelihood1{target}(:,measurement);
                else
                    inputDA(2,target) = currentExistenceExtrinsic(target,measurement) * (weightsExtrinsic{target}(:,measurement)'*likelihood1{target}(:,measurement));
                end
                inputDA(1,target) = 1;
            end
            targetIndex = numLegacy;
            targetIndexesCurrent = nan(numLegacy,1);

            % only new targets with index >= measurement index are connected to measurement
            for target = numMeasurements:-1:measurement
                if(any(target==newIndexes))
                    targetIndex = targetIndex + 1;
                    targetIndexesCurrent = [targetIndexesCurrent;target];
                    if(outer == 1)
                        likelihoodNew1{targetIndex-numLegacy}(:,measurement) = constantFactor * exp(getLogWeightsFast(measurements(:,measurement),reshape(multiBernoulli(targetIndex).particlesKinematic(1:2,step,:),2,numParticlesBernoulli(targetIndex)),reshape(multiBernoulli(targetIndex).particlesExtent(:,:,step,:),2,2,numParticlesBernoulli(targetIndex)) + repmat(measurementsCovariance,[1,1,numParticlesBernoulli(targetIndex)])));
                        inputDA(2,targetIndex) = currentExistenceExtrinsic(targetIndex,measurement) * (multiBernoulli(targetIndex).particleWeights'*likelihoodNew1{targetIndex-numLegacy}(:,measurement));
                    else
                        inputDA(2,targetIndex) = currentExistenceExtrinsic(targetIndex,measurement) * (weightsExtrinsicNew{targetIndex-numLegacy}(:,measurement)'*likelihoodNew1{targetIndex-numLegacy}(:,measurement));
                    end
                    inputDA(1,targetIndex) = 1;
                    if(target == measurement)
                        inputDA(1,targetIndex) = 1 - currentExistenceExtrinsic(targetIndex,measurement);
                    end
                end
            end
            targetIndexes{measurement} = targetIndexesCurrent;
            outputDA{measurement} = dataAssociationBP(inputDA);
        end

        % perform update step for legacy targets
        for target = 1:numLegacy
            weights = zeros(numParticlesBernoulli(target),numMeasurements+1);
            weights(:,numMeasurements+1) = log(multiBernoulli(target).particleWeights);
            for measurement = 1:numMeasurements
                currentWeights = 1 + likelihood1{target}(:,measurement) * outputDA{measurement}(1,target);
                currentWeights = log(currentWeights);
                weights(:,measurement) = currentWeights;
            end
            % calculate extrinsic information for legacy targets (at all except last iteration) and belief (at last iteration)
            if(outer ~= numOuterIterations)
                for measurement = 1:numMeasurements
                    [weightsExtrinsic{target}(:,measurement),currentExistenceExtrinsic(target,measurement)] = getWeightsUnknownNewAll(weights,multiBernoulli(target).existence,measurement);
                end
            else
                multiBernoulli(target) = updateParticlesNewBernoulliAllTrajectory(multiBernoulli(target),weights);
                % compute probability mass function of end time
                uniqueTrajectoryEndTime = unique(multiBernoulli(target).particleEndTime);
                cardinalityOfSupportOfEndTime = length(uniqueTrajectoryEndTime);
                probabilityMassOfEndTime = zeros(cardinalityOfSupportOfEndTime,1);
                indexTruncate = true(numParticlesBernoulli(target),1);
                for endTime = 1:cardinalityOfSupportOfEndTime
                    indexTrajectoryEndTime = multiBernoulli(target).particleEndTime == uniqueTrajectoryEndTime(endTime);
                    probabilityMassOfEndTime(endTime) = sum(multiBernoulli(target).particleWeights(indexTrajectoryEndTime));
                    if probabilityMassOfEndTime(endTime) >= thresholdTruncate
                        indexTruncate(indexTrajectoryEndTime) = false;
                    end
                end

                % store marginal density
                binaryIndexAliveParticle = multiBernoulli(target).particleEndTime == step;
                weightsReNormalized = multiBernoulli(target).particleWeights(binaryIndexAliveParticle);
                if ~isempty(weightsReNormalized)
                    multiBernoulli(target).weightHistory{step} = weightsReNormalized/sum(weightsReNormalized);
                    multiBernoulli(target).kinematicHistory{step} = reshape(multiBernoulli(target).particlesKinematic(:,step,binaryIndexAliveParticle),4,length(weightsReNormalized));
                    multiBernoulli(target).extentHistory{step} = reshape(multiBernoulli(target).particlesExtent(:,:,step,binaryIndexAliveParticle),2,2,length(weightsReNormalized));
                end

                % truncation
                multiBernoulli(target).particleWeights(indexTruncate) = [];
                multiBernoulli(target).particleStartTime(indexTruncate) = [];
                multiBernoulli(target).particleEndTime(indexTruncate) = [];
                multiBernoulli(target).particlesKinematic(:,:,indexTruncate) = [];
                multiBernoulli(target).particlesExtent(:,:,:,indexTruncate) = [];
                sumParticleWeights = sum(multiBernoulli(target).particleWeights);
                multiBernoulli(target).particleWeights = multiBernoulli(target).particleWeights/sumParticleWeights;
                multiBernoulli(target).existence = multiBernoulli(target).existence*sumParticleWeights;

                % resampling
                multiBernoulli(target) = resampleBernoulli(multiBernoulli(target),numParticles,step);
            end
        end

        % perform update step for new targets
        targetIndex = numLegacy;
        for target = numMeasurements:-1:1
            if(any(target == newIndexes))
                targetIndex = targetIndex + 1;
                weights = zeros(numParticlesBernoulli(targetIndex),numMeasurements+1);
                weights(:,numMeasurements+1) = log(multiBernoulli(targetIndex).particleWeights);
                for measurement = 1:target
                    outputTmpDA = outputDA{measurement}(1,targetIndexes{measurement}==target);
                    if(~isinf(outputTmpDA))
                        currentWeights = likelihoodNew1{targetIndex-numLegacy}(:,measurement) * outputTmpDA;
                    else
                        currentWeights = likelihoodNew1{targetIndex-numLegacy}(:,measurement);
                    end
                    if(measurement ~= target)
                        currentWeights = currentWeights + 1;
                    end
                    currentWeights = log(currentWeights);
                    weights(:,measurement) = currentWeights;
                end

                % calculate extrinsic information for new targets (at all except last iteration) or belief (at last iteration)
                if(outer ~= numOuterIterations)
                    for measurement = 1:target
                        [weightsExtrinsicNew{targetIndex-numLegacy}(:,measurement),currentExistenceExtrinsic(targetIndex,measurement)] = getWeightsUnknownNewAll(weights,multiBernoulli(targetIndex).existence,measurement);
                    end
                else
                    multiBernoulli(targetIndex) = updateParticlesNewBernoulliAllTrajectory(multiBernoulli(targetIndex),weights);
                    if numNewBornParticles > numParticles
                        multiBernoulli(targetIndex) = resampleBernoulli(multiBernoulli(targetIndex),numParticles,step);
                        numParticlesBernoulli(targetIndex) = numParticles;
                    end
                    % store marginal density
                    multiBernoulli(targetIndex).weightHistory{step} = multiBernoulli(targetIndex).particleWeights;
                    multiBernoulli(targetIndex).kinematicHistory{step} = reshape(multiBernoulli(targetIndex).particlesKinematic(:,step,:),4,numParticles);
                    multiBernoulli(targetIndex).extentHistory{step} = reshape(multiBernoulli(targetIndex).particlesExtent(:,:,step,:),2,2,numParticles);
                end
            end
        end
    end

    % perform pruning
    isNotRedundant = [multiBernoulli.existence] >= thresholdPruning;
    multiBernoulli = multiBernoulli(isNotRedundant);

    % misdetection of PPP
    PoissonPointProcess.meanNumOfUndetectedObject = PoissonPointProcess.meanNumOfUndetectedObject*(1-probabilityEffectiveDetect);

    elapsedTime(step) = toc;

    % compute trajectory metric
    [estimatedTracks,estimatedExtents] = allTrajectoryFormation(multiBernoulli,step,parameters);
    [d_traMetric,loc_cost,miss_cost,fa_cost,switch_cost] = computeAllTrajectoryMetric(targetTracks,targetExtents,estimatedTracks,estimatedExtents,step,parameters);
    trajectoryMetric.total(step) = d_traMetric;
    trajectoryMetric.loc(step) = sum(loc_cost);
    trajectoryMetric.miss(step) = sum(miss_cost);
    trajectoryMetric.false(step) = sum(fa_cost);
    trajectoryMetric.switch(step) = sum(switch_cost);
end

fprintf('\n')

runTime = sum(elapsedTime);

% extract particles used for backward simulation
estimatedParticles = extractParticles(multiBernoulli,numSteps,parameters);

end