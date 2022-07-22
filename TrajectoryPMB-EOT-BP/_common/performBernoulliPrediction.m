% Florian Meyer, 2020

function multiBernoulli = performBernoulliPrediction( multiBernoulli, parameters )
numTargets = length(multiBernoulli);
drivingNoiseVariance = parameters.accelerationDeviation^2;
survivalProbability = parameters.survivalProbability;
degreeFreedomPrediction = parameters.degreeFreedomPrediction;
scanTime = parameters.scanTime;

[A, Q] = getTransitionMatricesNew(scanTime);

for target = 1:numTargets
    numParticles = length(multiBernoulli(target).particleWeights);
    multiBernoulli(target).existence = survivalProbability*multiBernoulli(target).existence;
    multiBernoulli(target).particlesKinematic = mvnrnd((A*multiBernoulli(target).particlesKinematic)',drivingNoiseVariance*Q)';
    multiBernoulli(target).particlesExtent = wishrndFastVector(multiBernoulli(target).particlesExtent/degreeFreedomPrediction,degreeFreedomPrediction,numParticles);
end

end