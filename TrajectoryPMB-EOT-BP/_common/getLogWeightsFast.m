% Florian Meyer, 2020

function logWeights = getLogWeightsFast(measurement,currentParticlesKinematic,currentParticlesExtent)
numParticles = size(currentParticlesExtent,3);

allDeterminantes = permute(currentParticlesExtent(1,1,:).*currentParticlesExtent(2,2,:) - currentParticlesExtent(1,2,:).^2,[3,1,2]);
allFactors = log(1./(2*pi*sqrt(allDeterminantes)));

measurementsReptition = repmat(measurement,[1,numParticles]);

part2 = ( measurementsReptition - currentParticlesKinematic(1:2,:) )';

% direct calculation of innovation vector times inverted covariance matrix
tmp = 1./repmat(allDeterminantes,[1,2]) .* (measurementsReptition' - currentParticlesKinematic(1:2,:)');
part1(:,1) = tmp(:,1) .* squeeze(currentParticlesExtent(2,2,:)) - tmp(:,2) .* squeeze(currentParticlesExtent(2,1,:));
part1(:,2) = - tmp(:,1) .* squeeze(currentParticlesExtent(1,2,:)) + tmp(:,2) .* squeeze(currentParticlesExtent(1,1,:));

logWeights = allFactors + ( -1/2*(part1(:,1).*part2(:,1) + part1(:,2).*part2(:,2)) );
end