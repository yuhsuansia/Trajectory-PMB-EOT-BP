% Florian Meyer, 2020

function [centralIndexes,indexesReordered] = getCentralReordered(measurements,probabilitiesNew,measurementsCovariance,parameters)
threshold = parameters.freeThreshold;
clusterThreshold = parameters.clusterThreshold;
meanExtentBirth = parameters.priorExtent1/(parameters.priorExtent2-3);
measurementsCovariance = measurementsCovariance + meanExtentBirth;
minClusterElements = parameters.minClusterElements;

% determine which measurements are not associated to a legacy target (free) with high probability

allIndexesNumeric = (1:size(measurements,2))';

freeIndexes = probabilitiesNew >= threshold;
assignedIndexes = probabilitiesNew < threshold;

measurementsFree = measurements(:,freeIndexes);

freeIndexesNumeric = allIndexesNumeric(freeIndexes);
assignedIndexesNumeric = allIndexesNumeric(assignedIndexes);


% get clusters of free measurements using distance partitioning
clusters = getClusters(measurementsFree,measurementsCovariance,clusterThreshold)';

numElements = sum(clusters > 0,1);
[numElements,indexes] = sort(numElements,'descend');
clusters = clusters(:,indexes);

notUsedIndexes = clusters(:,numElements < minClusterElements);
notUsedIndexes = nonzeros(notUsedIndexes(:));
notUsedIndexesNumeric = freeIndexesNumeric(notUsedIndexes);
numNotUsed = size(notUsedIndexesNumeric,1);

clusters(:,numElements < minClusterElements) = [];
indexesNumericNew = zeros(0,1);
numClusters = size(clusters,2);
centralIndexes = zeros(numClusters,1);
for cluster = 1:numClusters
    
    indexes = nonzeros(clusters(:,cluster));
    currentMeasurements = measurementsFree(:,indexes);
    
    currentIndexesNumeric = freeIndexesNumeric(indexes);
    numMeasurements = size(indexes,1);
    
    if(numel(indexes)>1)
        
        distanceMatrix = zeros(numMeasurements,numMeasurements);
        for measurement1 = 1:numMeasurements
            for measurement2 = (measurement1+1):numMeasurements
                distVector = currentMeasurements(:,measurement1)-currentMeasurements(:,measurement2);
                
                distanceMatrix(measurement1,measurement2) = sqrt(distVector'/measurementsCovariance*distVector);
                distanceMatrix(measurement2,measurement1) = distanceMatrix(measurement1,measurement2);
            end
        end
        
        % order measurements within cluster according to sum of all distances to other measurements in cluster
        distanceVector = sum(distanceMatrix,2);
        [~,indexes] = sort(distanceVector,'descend');
        currentIndexesNumeric = currentIndexesNumeric(indexes);

    end
    indexesNumericNew = [currentIndexesNumeric;indexesNumericNew];
    
    % central measurement is last
    centralIndexes(1:cluster) = centralIndexes(1:cluster) + numMeasurements;
end

% put measurements of legacy targets first; then clusters of free measurements in a cluster by cluster manner; central measurement in each cluster is first
indexesReordered = [notUsedIndexesNumeric;indexesNumericNew;assignedIndexesNumeric]; 

centralIndexes = centralIndexes + numNotUsed;
centralIndexes = sort(centralIndexes,'descend');
end



function [clusters] = getClusters(measurements,measurementsCovariance,thresholdProbability)
numMeasurements = size(measurements,2);

if(~numMeasurements)
    clusters = [];
    return;
end

thresholdDistance = chi2inv(thresholdProbability,2);

% get distances among measurements
distanceVector = zeros((numMeasurements*(numMeasurements-1)/2+1),1);
distanceMatrix = zeros(numMeasurements,numMeasurements);
entry = 1;
for measurement1 = 1:numMeasurements
    for measurement2 = (measurement1+1):numMeasurements
        distVector = measurements(:,measurement1)-measurements(:,measurement2);
        
        entry = entry + 1;
        distanceVector(entry) = sqrt(distVector'/measurementsCovariance*distVector);
        
        distanceMatrix(measurement1,measurement2) = distanceVector(entry);
        distanceMatrix(measurement2,measurement1) = distanceVector(entry);
    end
end

% extract distance vectors used for partitioning
distanceVector = sort(distanceVector);
distanceVector(distanceVector>thresholdDistance) = [];
distance = distanceVector(end);


% get clusters
clusterNumbers = zeros(numMeasurements,1);
clusterId = 1;
for measurement = 1:numMeasurements
    if(clusterNumbers(measurement)==0)
        clusterNumbers(measurement) = clusterId;
        clusterNumbers = findNeighbors(measurement,clusterNumbers,clusterId,distanceMatrix,distance);
        clusterId = clusterId + 1;
    end
end
numClusters = clusterId - 1;

maxElements = sum(clusterNumbers==mode(clusterNumbers));
clusters = zeros(0,maxElements);
index = 0;
for cluster = 1:numClusters
    associationTmp = find(clusterNumbers==cluster)';
    numElements = numel(associationTmp);
    if(numElements <= maxElements)
        index = index + 1;
        clusters(index,:) = [zeros(1,maxElements-numElements),associationTmp];
    end
end

end



function [cellNumbers] = findNeighbors(index,cellNumbers,cellId,distanceMatrix,distanceThreshold)
numMeasurements = size(distanceMatrix,2);

for measurement = 1:numMeasurements
    if(measurement ~= index && distanceMatrix(measurement,index) < distanceThreshold && cellNumbers(measurement) == 0)
        cellNumbers(measurement) = cellId;
        cellNumbers = findNeighbors(index,cellNumbers,cellId,distanceMatrix,distanceThreshold);
    end
end

end