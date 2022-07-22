% Florian Meyer, 2020
% Modified by Yuxuan for further vectorization
function logWeights = getLogWeightsFastMatrix(measurements,currentParticlesKinematic,currentParticlesExtent)
numParticles = size(currentParticlesExtent,3);

allDeterminantes = permute(currentParticlesExtent(1,1,:).*currentParticlesExtent(2,2,:) - currentParticlesExtent(1,2,:).^2,[3,1,2]);
allFactors = log(1./(2*pi*sqrt(allDeterminantes)));

tmp1 = 1./repmat(allDeterminantes,[1,2]);
squeeze_extent11 = squeeze(currentParticlesExtent(1,1,:));
squeeze_extent12 = squeeze(currentParticlesExtent(1,2,:));
squeeze_extent21 = squeeze(currentParticlesExtent(2,1,:));
squeeze_extent22 = squeeze(currentParticlesExtent(2,2,:));

numMeasurements = size(measurements,2);
logWeights = zeros(numParticles,numMeasurements);

for j = 1:numMeasurements
    measurementsReptition = repmat(measurements(:,j),[1,numParticles]);

    part2 = ( measurementsReptition - currentParticlesKinematic(1:2,:) )';

    % direct calculation of innovation vector times inverted covariance matrix
    tmp2 = tmp1 .* (measurementsReptition' - currentParticlesKinematic(1:2,:)');
    part1(:,1) = tmp2(:,1) .* squeeze_extent22 - tmp2(:,2) .* squeeze_extent21;
    part1(:,2) = - tmp2(:,1) .* squeeze_extent12 + tmp2(:,2) .* squeeze_extent11;

    logWeights(:,j) = allFactors + ( -1/2*(part1(:,1).*part2(:,1) + part1(:,2).*part2(:,2)) );
end
end