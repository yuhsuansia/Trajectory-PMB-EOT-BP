% Florian Meyer, 2020

function [a] = wishrndFastVector(parameter1,parameter2,numParticles)
d = zeros(size(parameter1));
d(1,1,:) = sqrt(parameter1(1,1,:));
d(2,1,:) = parameter1(2,1,:) ./ d(1,1,:);
d(2,2,:) = sqrt(parameter1(2,2,:) - d(2,1,:).^2);
d = permute(d,[2,1,3]);

if(size(parameter1,3) == 1)
    d = repmat(d,[1,1,numParticles]);
end

r = 2.*randg((parameter2 - [zeros(1,numParticles);ones(1,numParticles)])./2);
x = sqrt(r);
x = [x;randn(1,numParticles)];
x = permute(x,[1,3,2]);

T = zeros(2,2,numParticles);
T(1,1,:) = x(1,1,:) .* d(1,1,:) + x(3,1,:) .* d(2,1,:);
T(1,2,:) = x(1,1,:) .* d(1,2,:) + x(3,1,:) .* d(2,2,:);
T(2,2,:) = x(2,1,:) .* d(2,2,:);

a = zeros(2,2,numParticles);
a(1,1,:) = T(1,1,:) .* T(1,1,:) + T(2,1,:) .* T(2,1,:);
a(1,2,:) = T(1,1,:) .* T(1,2,:) + T(2,1,:) .* T(2,2,:);
a(2,1,:) = a(1,2,:);
a(2,2,:) = T(1,2,:) .* T(1,2,:) + T(2,2,:) .* T(2,2,:);

end