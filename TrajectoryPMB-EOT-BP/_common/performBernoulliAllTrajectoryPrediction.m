function multiBernoulli = performBernoulliAllTrajectoryPrediction( multiBernoulli, scanTime, step, parameters )
numTargets = length(multiBernoulli);
drivingNoiseVariance = parameters.accelerationDeviation^2;
survivalProbability = parameters.survivalProbability;
degreeFreedomPrediction = parameters.degreeFreedomPrediction;

[A, Q] = getTransitionMatricesNew(scanTime);

for target = 1:numTargets
    numParticles = length(multiBernoulli(target).particleWeights);

    binaryIndexAliveParticles = multiBernoulli(target).particleEndTime == step-1;
    IndexAliveParticles = find(binaryIndexAliveParticles)';

    multiBernoulli(target).particleWeights = [multiBernoulli(target).particleWeights;multiBernoulli(target).particleWeights(binaryIndexAliveParticles)*(1-survivalProbability)];
    multiBernoulli(target).particleStartTime = multiBernoulli(target).particleStartTime([1:numParticles IndexAliveParticles]);
    multiBernoulli(target).particleEndTime = multiBernoulli(target).particleEndTime([1:numParticles IndexAliveParticles]);

    multiBernoulli(target).particlesKinematic = multiBernoulli(target).particlesKinematic(:,:,[1:numParticles IndexAliveParticles]);
    multiBernoulli(target).particlesExtent = multiBernoulli(target).particlesExtent(:,:,:,[1:numParticles IndexAliveParticles]);

    if ~isempty(IndexAliveParticles)
        multiBernoulli(target).particleWeights(IndexAliveParticles) = multiBernoulli(target).particleWeights(IndexAliveParticles)*survivalProbability;
        multiBernoulli(target).particleEndTime(IndexAliveParticles) = multiBernoulli(target).particleEndTime(IndexAliveParticles)+1;

        multiBernoulli(target).particlesKinematic(:,step,IndexAliveParticles) = mvnrnd((A*reshape(multiBernoulli(target).particlesKinematic(:,step-1,IndexAliveParticles),4,length(IndexAliveParticles)))',drivingNoiseVariance*Q)';
        multiBernoulli(target).particlesExtent(:,:,step,IndexAliveParticles) = wishrndFastVector(reshape(multiBernoulli(target).particlesExtent(:,:,step-1,IndexAliveParticles),2,2,length(IndexAliveParticles))/degreeFreedomPrediction,degreeFreedomPrediction,length(IndexAliveParticles));
    end
end

end