function Bernoulli = resampleBernoulliTargets(Bernoulli,numParticles)

indexes = resampleSystematic(Bernoulli.particleWeights,numParticles);
Bernoulli.particleWeights = 1/numParticles*ones(numParticles,1);
Bernoulli.particlesKinematic = Bernoulli.particlesKinematic(:,indexes);
Bernoulli.particlesExtent = Bernoulli.particlesExtent(:,:,indexes);

end