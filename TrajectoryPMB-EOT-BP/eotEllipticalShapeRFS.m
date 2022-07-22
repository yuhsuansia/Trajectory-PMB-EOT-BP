function [ estimatedTracks, estimatedExtents, estimatedParticles, trajectoryMetric, runTime] = eotEllipticalShapeRFS( measurementsCell,targetTracks,targetExtents,parameters )

numParticles = parameters.numParticles;
numNewBornParticles = parameters.numNewBornParticles;

meanClutter = parameters.meanClutter;
meanMeasurements = parameters.meanMeasurements;
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

detectionThreshold = parameters.detectionThreshold;
thresholdPruning = parameters.thresholdPruning;
numOuterIterations = parameters.numOuterIterations;

[numSteps, ~] = size(measurementsCell);
constantFactor = areaSize*(meanMeasurements/meanClutter);
probabilityEffectiveDetect = 1-exp(-meanMeasurements);
uniformWeight = log(1/areaSize);

multiBernoulli = repmat(struct('mark',[],'existence',[],'particleWeights',[],'particlesKinematic',[],'particlesExtent',[]),[0,1]);
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

    % estimates
    detectedTargets = 0;
    estimates{step}.wp = [];
    estimates{step}.xp = [];
    estimates{step}.ep = [];
    estimates{step}.state = zeros(4,0);
    estimates{step}.extent = zeros(2,2,0);
    estimates{step}.label = zeros(2,0);

    % load current measurements
    measurements = measurementsCell{step};
    numMeasurements = size(measurements,2);
    
    % perform prediction step
    multiBernoulli = performBernoulliPrediction(multiBernoulli,parameters);
    numLegacy = length(multiBernoulli);
    for target = 1:numLegacy
        multiBernoulli(target).existence = multiBernoulli(target).existence*(1-probabilityEffectiveDetect)/(1-multiBernoulli(target).existence*probabilityEffectiveDetect);
    end
    PoissonPointProcess.meanNumOfUndetectedObject = survivalProbability*PoissonPointProcess.meanNumOfUndetectedObject+meanBirths;
    
    % get indexes of promising new objects
    [newIndexes,measurements,likelihood1] = getPromisingNewMeasurement(multiBernoulli,measurements,parameters);
    numNew = size(newIndexes,1);
    newMultiBernoulli = repmat(struct('mark',[],'existence',[],'particleWeights',[],'particlesKinematic',[],'particlesExtent',[]),[numNew,1]);
    for target = 1:numNew
        proposalMean = measurements(:,newIndexes(target));
        proposalCovariance = 2 * totalCovariance; % stretch covariance matrix to make proposal distribution heavier-tailed then target distribution
        newMultiBernoulli(target).particlesKinematic = zeros(4,numNewBornParticles);
        newMultiBernoulli(target).particlesKinematic(1:2,:) = proposalMean + sqrtm(proposalCovariance) * randn(2,numNewBornParticles);
        newMultiBernoulli(target).particlesKinematic(3:4,:) = mvnrnd([0;0],priorVelocityCovariance,numNewBornParticles)';
        newMultiBernoulli(target).particlesExtent = iwishrndFastVector(priorExtent1,priorExtent2,numNewBornParticles);
        newMultiBernoulli(target).particleWeights = exp(normalizeLogWeights(uniformWeight-log_mvnpdf(reshape(newMultiBernoulli(target).particlesKinematic(1:2,:),2,numNewBornParticles)', proposalMean', proposalCovariance')));
        aliveProbability = PoissonPointProcess.meanNumOfUndetectedObject*(1-probabilityEffectiveDetect);
        newMultiBernoulli(target).existence = aliveProbability/(aliveProbability+1);
        newMultiBernoulli(target).mark = [step;newIndexes(target)];
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
                        likelihoodNew1{targetIndex-numLegacy}(:,measurement) = constantFactor * exp(getLogWeightsFast(measurements(:,measurement),reshape(multiBernoulli(targetIndex).particlesKinematic(1:2,:),2,numParticlesBernoulli(targetIndex)),reshape(multiBernoulli(targetIndex).particlesExtent,2,2,numParticlesBernoulli(targetIndex)) + repmat(measurementsCovariance,[1,1,numParticlesBernoulli(targetIndex)])));
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
                multiBernoulli(target) = updateParticlesNewBernoulli(multiBernoulli(target),weights);

                % perform estimation
                if(multiBernoulli(target).existence > detectionThreshold)
                    detectedTargets = detectedTargets + 1;
                    Bernoulli = multiBernoulli(target);
                    estimates{step}.state(:,detectedTargets) = sum(Bernoulli.particlesKinematic.*Bernoulli.particleWeights',2);
                    particleWeightReshape = ones(1,1,numParticlesBernoulli(target));
                    particleWeightReshape(1,1,:) = Bernoulli.particleWeights;
                    estimates{step}.extent(:,:,detectedTargets) = sum(Bernoulli.particlesExtent.*particleWeightReshape,3);
                    estimates{step}.label(:,detectedTargets) = Bernoulli.mark;
                    estimates{step}.wp{detectedTargets} = Bernoulli.particleWeights;
                    estimates{step}.xp{detectedTargets} = Bernoulli.particlesKinematic;
                    estimates{step}.ep{detectedTargets} = Bernoulli.particlesExtent;
                end
                
                % resampling
                multiBernoulli(target) = resampleBernoulliTargets(multiBernoulli(target),numParticles);
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
                    multiBernoulli(targetIndex) = updateParticlesNewBernoulli(multiBernoulli(targetIndex),weights);
                    if numNewBornParticles > numParticles
                        multiBernoulli(targetIndex) = resampleBernoulliTargets(multiBernoulli(targetIndex),numParticles);
                        numParticlesBernoulli(targetIndex) = numParticles;
                    end
                end

                % perform estimation
                if(multiBernoulli(targetIndex).existence > detectionThreshold)
                    detectedTargets = detectedTargets + 1;
                    Bernoulli = multiBernoulli(targetIndex);
                    estimates{step}.state(:,detectedTargets) = sum(Bernoulli.particlesKinematic.*Bernoulli.particleWeights',2);
                    particleWeightReshape = ones(1,1,numParticlesBernoulli(targetIndex));
                    particleWeightReshape(1,1,:) = Bernoulli.particleWeights;
                    estimates{step}.extent(:,:,detectedTargets) = sum(Bernoulli.particlesExtent.*particleWeightReshape,3);
                    estimates{step}.label(:,detectedTargets) = Bernoulli.mark;
                    estimates{step}.wp{detectedTargets} = Bernoulli.particleWeights;
                    estimates{step}.xp{detectedTargets} = Bernoulli.particlesKinematic;
                    estimates{step}.ep{detectedTargets} = Bernoulli.particlesExtent;
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
    [estimatedTracks,estimatedExtents] = trackFormation(estimates, parameters);
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
[~,~,estimatedParticles] = trackFormation(estimates, parameters);

end