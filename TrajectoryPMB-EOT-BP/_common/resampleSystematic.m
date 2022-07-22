% Florian Meyer, 2019

function indexes = resampleSystematic(weights,numParticles)
indexes = zeros(numParticles,1);
cumWeights = cumsum(weights);

grid = zeros(1,numParticles+1);
grid(1:numParticles) = linspace(0,1-1/numParticles,numParticles) + rand/numParticles;
grid(numParticles+1) = 1;

i = 1;
j = 1;
while( i <= numParticles )
    if( grid(i) < cumWeights(j) )
        indexes(i) = j;
        i = i + 1;
    else
        j = j + 1;
    end
end

end


