function [clutteredMeasurements] = generateClutteredMeasurements(targetTracks, targetExtents, parameters)
meanClutter = parameters.meanClutter;
measurementVariance = parameters.measurementVariance;
meanMeasurements = parameters.meanMeasurements;
surveillanceRegion = parameters.surveillanceRegion;

numSteps = size(targetTracks,2);
numTargets = size(targetTracks,3);

clutteredMeasurements = cell(numSteps,1);
for step = 1:numSteps
    measurements = nan(2,0);
    numMeasurements = 0;
    for target = 1:numTargets
        if(isnan(targetTracks(1,step,target)))
            continue
        end
        numMeasurementsTmp = poissrnd(meanMeasurements);
        measurementsTmp1 = repmat(targetTracks(1:2,step,target),[1,numMeasurementsTmp]) + mvnrnd([0 0],targetExtents(:,:,step,target),numMeasurementsTmp)'; 
        
        measurementsTmp1 = measurementsTmp1 + sqrt(measurementVariance)*randn(2,numMeasurementsTmp);
        measurements = cat(2,measurements,measurementsTmp1);
        numMeasurements = numMeasurements + numMeasurementsTmp;
    end  
    
    numFalseAlarms = poissrnd(meanClutter);
    if(isempty(numFalseAlarms))
        numFalseAlarms = 0;
    end
    
    falseAlarms = zeros(2,numFalseAlarms);
    falseAlarms(1,:) = (surveillanceRegion(2,1) - surveillanceRegion(1,1))*rand(1,numFalseAlarms) + surveillanceRegion(1,1);
    falseAlarms(2,:) = (surveillanceRegion(2,2) - surveillanceRegion(1,2))*rand(1,numFalseAlarms) + surveillanceRegion(1,2);
    
    measurements = [falseAlarms, measurements];
    measurements = measurements(:,randperm(numFalseAlarms+numMeasurements));
    
    clutteredMeasurements{step} = measurements;
end

end