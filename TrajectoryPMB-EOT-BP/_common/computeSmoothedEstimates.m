function [smoothedTracks,smoothedExtents] = computeSmoothedEstimates(estimatedParticles,parameters)

[numSteps,numTargets] = size(estimatedParticles);
numIter = parameters.backwardSimulation;
[A,Q] = getTransitionMatricesNew(parameters.scanTime);

smoothedTracks = nan(4,numSteps,numTargets);
smoothedExtents = nan(2,2,numSteps,numTargets);

for target = 1:numTargets
    x_ffbsi = nan(4,numIter,numSteps);
    e_ffbsi = nan(2,2,numIter,numSteps);

    trajectoryStartTime = find(cellfun(@(x) ~isempty(x), estimatedParticles(:,target)),1,'first');
    trajectoryEndTime = find(cellfun(@(x) ~isempty(x), estimatedParticles(:,target)),1,'last');
    
    bins = [0 cumsum(estimatedParticles{trajectoryEndTime,target}.wp)'];
    [~,b] = histc(rand(numIter,1),bins);
    x_ffbsi(:,:,trajectoryEndTime) = estimatedParticles{trajectoryEndTime,target}.xp(:,b);
    e_ffbsi(:,:,:,trajectoryEndTime) = estimatedParticles{trajectoryEndTime,target}.ep(:,:,b);

    for timeStep = trajectoryEndTime-1:-1:trajectoryStartTime
        for iter = 1:numIter
            if isempty(estimatedParticles{timeStep,target})
                x_ffbsi(:,iter,timeStep) = A\x_ffbsi(:,iter,timeStep+1);
                e_ffbsi(:,:,iter,timeStep) = e_ffbsi(:,:,iter,timeStep+1);
            else
                w_tilde = log(estimatedParticles{timeStep,target}.wp)'+log_mvnpdf((A*estimatedParticles{timeStep,target}.xp)',x_ffbsi(:,iter,timeStep+1)',Q)';
                w_tilde = exp(normalizeLogWeights(w_tilde));
                bins = [0 cumsum(w_tilde)];
                b = find(bins>rand,1)-1;
                x_ffbsi(:,iter,timeStep) = estimatedParticles{timeStep,target}.xp(:,b);
                e_ffbsi(:,:,iter,timeStep) = estimatedParticles{timeStep,target}.ep(:,:,b);
            end
        end
    end

    for timeStep = trajectoryStartTime:trajectoryEndTime
        smoothedTracks(:,timeStep,target) = mean(x_ffbsi(:,:,timeStep),2);
        smoothedExtents(:,:,timeStep,target) = mean(e_ffbsi(:,:,:,timeStep),3);
    end
end

end