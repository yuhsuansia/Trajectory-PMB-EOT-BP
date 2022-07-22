function [newIndexes,measurements,likelihood1] = getPromisingNewMeasurement(multiBernoulli,measurements,parameters)

% find ``free'' measurements by updating legacy targets only
numMeasurements = size(measurements,2);
numTargets = length(multiBernoulli);
measurementsCovariance = parameters.measurementVariance * eye(2);
meanClutter = parameters.meanClutter;
meanMeasurements = parameters.meanMeasurements;
surveillanceRegion = parameters.surveillanceRegion;
areaSize =  (surveillanceRegion(2,1)-surveillanceRegion(1,1)) * (surveillanceRegion(2,2)-surveillanceRegion(1,2));
constantFactor = areaSize*(meanMeasurements/meanClutter);
currentEpsilonExtrinsic = repmat([multiBernoulli.existence]',[1,numMeasurements]);

numParticlesBernoulli = arrayfun(@(x) length(x.particleWeights),multiBernoulli);
inputDA = zeros(numTargets,numMeasurements);

likelihood1 = cell(numTargets,1);

for target = 1:numTargets
    lengthIndex = numParticlesBernoulli(target);
    likelihood1{target} = zeros(lengthIndex,numMeasurements);
    tmpMean = reshape(multiBernoulli(target).particlesKinematic(1:2,:),2,lengthIndex);
    tmpCov = reshape(multiBernoulli(target).particlesExtent,2,2,lengthIndex) + repmat(measurementsCovariance,[1,1,lengthIndex]);
    for measurement = numMeasurements:-1:1
        likelihood1{target}(:,measurement) = constantFactor * exp(getLogWeightsFast(measurements(:,measurement),tmpMean,tmpCov));
        inputDA(target,measurement) = currentEpsilonExtrinsic(target,measurement) * multiBernoulli(target).particleWeights'*likelihood1{target}(:,measurement);
    end
end

probabilitiesNew = zeros(numMeasurements,numTargets+1);
for measurement = 1:numMeasurements
    probabilitiesNew(measurement,:) = [inputDA(:,measurement);1]/(sum(inputDA(:,measurement))+1);
end
probabilitiesNew = probabilitiesNew(:,end);

% find central measurements in clusters of unused measurements and reorder measurement
[newIndexes,indexesReordered] = getCentralReordered(measurements,probabilitiesNew,measurementsCovariance,parameters);
measurements = measurements(:,indexesReordered);
likelihood1 = cellfun(@(x) x(:,indexesReordered),likelihood1,'UniformOutput',false);

end