function Bernoulli = resampleBernoulli(Bernoulli,numParticles,step)

binaryIndexAliveParticle = Bernoulli.particleEndTime == step;
indexAliveParticle = find(binaryIndexAliveParticle);
indexNotAliveParticle = find(~binaryIndexAliveParticle);

weightAliveParticle = Bernoulli.particleWeights(indexAliveParticle);
numParticles = min(length(indexAliveParticle),numParticles);

indexes = resampleSystematic(weightAliveParticle/sum(weightAliveParticle),numParticles);

if ~isempty(indexAliveParticle)
    Bernoulli.particleWeights = [sum(weightAliveParticle)/numParticles*ones(numParticles,1);Bernoulli.particleWeights(indexNotAliveParticle)];
    Bernoulli.particleStartTime = Bernoulli.particleStartTime([indexAliveParticle(indexes);indexNotAliveParticle]);
    Bernoulli.particleEndTime = Bernoulli.particleEndTime([indexAliveParticle(indexes);indexNotAliveParticle]);
    Bernoulli.particlesKinematic = Bernoulli.particlesKinematic(:,:,[indexAliveParticle(indexes);indexNotAliveParticle]);
    Bernoulli.particlesExtent = Bernoulli.particlesExtent(:,:,:,[indexAliveParticle(indexes);indexNotAliveParticle]);
end

end